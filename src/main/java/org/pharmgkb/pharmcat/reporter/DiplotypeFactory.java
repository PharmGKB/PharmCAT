package org.pharmgkb.pharmcat.reporter;

import java.lang.invoke.MethodHandles;
import java.util.ArrayList;
import java.util.Collection;
import java.util.List;
import java.util.Set;
import java.util.SortedSet;
import java.util.TreeSet;
import com.google.common.base.Preconditions;
import com.google.common.collect.ImmutableList;
import com.google.common.collect.ImmutableSet;
import com.google.common.collect.Lists;
import org.apache.commons.lang3.StringUtils;
import org.pharmgkb.pharmcat.Env;
import org.pharmgkb.pharmcat.definition.model.NamedAllele;
import org.pharmgkb.pharmcat.haplotype.model.BaseMatch;
import org.pharmgkb.pharmcat.haplotype.model.CombinationMatch;
import org.pharmgkb.pharmcat.haplotype.model.DiplotypeMatch;
import org.pharmgkb.pharmcat.haplotype.model.GeneCall;
import org.pharmgkb.pharmcat.phenotype.model.GenePhenotype;
import org.pharmgkb.pharmcat.reporter.model.DataSource;
import org.pharmgkb.pharmcat.reporter.model.OutsideCall;
import org.pharmgkb.pharmcat.reporter.model.result.Diplotype;
import org.pharmgkb.pharmcat.reporter.model.result.GeneReport;
import org.pharmgkb.pharmcat.reporter.model.result.Haplotype;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import static org.pharmgkb.pharmcat.reporter.model.result.GeneReport.isSinglePloidy;


/**
 * Factory class for spitting out {@link Diplotype} objects from known information about the gene
 * <p>
 * This is designed to take information from either the Matcher or from outside calls and produce consistent results
 *
 * @author Ryan Whaley
 */
public class DiplotypeFactory {
  private static final Logger sf_logger = LoggerFactory.getLogger(MethodHandles.lookup().lookupClass());
  private static final Set<String> PHENOTYPE_ONLY = ImmutableSet.of("HLA-A", "HLA-B");

  private final String m_gene;
  private final Env m_env;

  /**
   * Public constructor.
   * <p>
   * Initialize the factory with all the necessary information to make haplotype and diplotype calls
   *
   * @param gene the gene symbol for the diplotypes this will call
   */
  public DiplotypeFactory(String gene, Env env) {
    m_gene = gene;
    m_env = env;
  }


  public List<Diplotype> makeComponentDiplotypes(GeneCall geneCall, DataSource source) {
    List<Diplotype> dips = new ArrayList<>();
    SortedSet<NamedAllele> haps = new TreeSet<>();
    for (BaseMatch bm : geneCall.getHaplotypes()) {
      if (bm instanceof CombinationMatch cm) {
        haps.addAll(cm.getComponentHaplotypes());
      } else {
        haps.add(bm.getHaplotype());
      }
    }
    for (NamedAllele na : haps) {
      dips.add(new Diplotype(m_gene, m_env.makeHaplotype(m_gene, na.getName(), source)));
    }
    return dips;
  }


  public List<Diplotype> makeDiplotypes(@SuppressWarnings("rawtypes") Collection matches, DataSource source) {
    if (matches.size() > 0) {
      List<Diplotype> diplotypes = new ArrayList<>();
      for (Object obj : matches) {
        if (obj instanceof BaseMatch bm) {
          diplotypes.add(makeDiplotype(bm, source));
        } else if (obj instanceof DiplotypeMatch dm) {
          diplotypes.add(makeDiplotype(dm, source));
        } else if (obj instanceof String name) {
          diplotypes.add(makeDiplotype(name, source));
        } else {
          throw new IllegalStateException("Unknown type: " + obj.getClass());
        }
      }
      return diplotypes;
    } else {
      return ImmutableList.of(makeUnknownDiplotype(m_gene, m_env, source));
    }
  }

  private Diplotype makeDiplotype(BaseMatch baseMatch, DataSource source) {
    Diplotype diplotype = new Diplotype(m_gene, m_env.makeHaplotype(m_gene, baseMatch.getName(), source));
    fillDiplotype(diplotype, m_env, source);
    return diplotype;
    }

  private Diplotype makeDiplotype(DiplotypeMatch diplotypeMatch, DataSource source) {
    BaseMatch h1 = diplotypeMatch.getHaplotype1();
    BaseMatch h2 = diplotypeMatch.getHaplotype2();

    Diplotype diplotype = new Diplotype(m_gene, m_env.makeHaplotype(m_gene, h1.getName(), source),
        m_env.makeHaplotype(m_gene, h2.getName(), source));
    fillDiplotype(diplotype, m_env, source);

    diplotype.setCombination(h1 instanceof CombinationMatch || h2 instanceof CombinationMatch);
    return diplotype;
  }


  /**
   * Make diplotype objects based on {@link OutsideCall} objects
   */
  public List<Diplotype> makeDiplotypes(OutsideCall outsideCall, DataSource source) {
    Preconditions.checkNotNull(outsideCall);
    Preconditions.checkArgument(outsideCall.getGene().equals(m_gene));

    List<Diplotype> diplotypes = new ArrayList<>();
    if (outsideCall.getDiplotypes().size() > 0) {
      diplotypes.addAll(makeDiplotypes(outsideCall.getDiplotypes(), source));
      for (Diplotype diplotype : diplotypes) {
        // assign outside phenotype if it exists and warn on mismatch
        if (outsideCall.getPhenotype() != null && diplotype.getPhenotypes().size() > 0 && !diplotype.getPhenotypes().contains(outsideCall.getPhenotype())) {
          sf_logger.warn("Outside call phenotype does not match known phenotype: {} != {}", outsideCall.getPhenotype(), String.join("; ", diplotype.getPhenotypes()));
          diplotype.setPhenotypes(ImmutableList.of(outsideCall.getPhenotype()));
        }
        // assign outside activity score if it exists and warn on mismatch
        if (outsideCall.getActivityScore() != null && diplotype.getActivityScore() != null && !diplotype.getActivityScore().equals(outsideCall.getActivityScore())) {
          sf_logger.warn("Outside call activity score does not match known activity score: {} != {}", outsideCall.getActivityScore(), diplotype.getActivityScore());
          diplotype.setActivityScore(outsideCall.getActivityScore());
        }
      }
    }
    else {
      Diplotype outsideDiplotype = new Diplotype(outsideCall);
      GenePhenotype gp = m_env.getPhenotype(m_gene, source);
      if (gp != null && gp.isMatchedByActivityScore()
          && outsideDiplotype.isUnknownPhenotype()
          && outsideDiplotype.hasActivityScore()) {
        outsideDiplotype.setPhenotypes(ImmutableList.of(gp.getPhenotypeForActivity(outsideDiplotype)));
      }
      fillDiplotype(outsideDiplotype, m_env, source);
      diplotypes.add(outsideDiplotype);
    }
    return diplotypes;
  }


  public static String[] splitDiplotype(String gene, String diplotypeText) {
    if (diplotypeText.contains("/")) {
      if (isSinglePloidy(gene) && !GeneReport.isXChromo(gene)) {
        throw new BadOutsideCallException("Cannot specify two genotypes [" + diplotypeText + "] for single chromosome gene " +
            gene);
      }
      String[] alleles = diplotypeText.split("/");
      if (alleles.length != 2) {
        throw new BadOutsideCallException("Diplotype for " + gene + " has " + alleles.length + " alleles: (" +
            diplotypeText + ")");
      }
      return alleles;
    } else {
      if (!isSinglePloidy(gene) && !isPhenotypeOnly(gene)) {
        throw new BadOutsideCallException("Expected two genotypes separated by a '/' but saw [" + diplotypeText + "] for " +
            gene);
      }
      return new String[] {diplotypeText};
    }
  }

  public static List<String> splitHaplotype(String haplotypeText) {
    if (haplotypeText.contains(CombinationMatch.COMBINATION_JOINER)) {
      return CombinationMatch.COMBINATION_NAME_SPLITTER.splitToList(haplotypeText);
    }
    return Lists.newArrayList(haplotypeText);
  }


  /**
   * Make a faux Diplotype based on a string in the form "*1/*20".
   *
   * @param diplotypeText a string in the form "*1/*20"
   * @return a Diplotype containing the specified haplotypes
   */
  public Diplotype makeDiplotype(String diplotypeText, DataSource source) {
    Preconditions.checkArgument(StringUtils.isNotBlank(diplotypeText));

    String[] alleles = splitDiplotype(m_gene, diplotypeText);
    Haplotype hap1 = m_env.makeHaplotype(m_gene, alleles[0], source);
    Haplotype hap2 = alleles.length == 2 ? m_env.makeHaplotype(m_gene, alleles[1], source) : null;
    Diplotype diplotype = new Diplotype(m_gene, hap1, hap2);
    fillDiplotype(diplotype, m_env, source);
    return diplotype;
  }


  private static boolean isPhenotypeOnly(String gene) {
    return PHENOTYPE_ONLY.contains(gene);
  }


  public static Diplotype makeUnknownDiplotype(String gene, Env env, DataSource source) {
    Haplotype haplotype = new Haplotype(gene, Haplotype.UNKNOWN);
    Diplotype diplotype;
    if (isSinglePloidy(gene)) {
      diplotype = new Diplotype(gene, haplotype);
    } else {
      diplotype = new Diplotype(gene, haplotype, haplotype);
    }
    fillDiplotype(diplotype, env, source);
    return diplotype;
  }

  /**
   * Fill in the phenotype and activity score depending on what information is already in the diplotype
   * @param diplotype the diplotype to add information to
   * @param env Global PharmCAT environment
   * @param source which data source is supplying the phenotype/activity information
   */
  public static void fillDiplotype(Diplotype diplotype, Env env, DataSource source) {
    if (diplotype.getGene().startsWith("HLA")) {
      if (!diplotype.isUnknownAlleles() && diplotype.isUnknownPhenotype()) {
        diplotype.setPhenotypes(makeHlaPhenotype(diplotype));
        diplotype.setLookupKeys(makeHlaPhenotype(diplotype));
      }
      else if (diplotype.isUnknownAlleles() && !diplotype.isUnknownPhenotype()) {
        diplotype.setLookupKeys(diplotype.getPhenotypes());
      }
      return;
    }
    GenePhenotype gp = env.getPhenotype(diplotype.getGene(), source);
    if (gp == null) {
      return;
    }
    if (gp.isMatchedByActivityScore() && diplotype.getActivityScore() == null) {
      diplotype.setActivityScore(gp.getActivityForDiplotype(diplotype));
    }
    if (diplotype.getPhenotypes().size() == 0) {
      if (diplotype.isUnknownAlleles()) {
        diplotype.addPhenotype(gp.getPhenotypeForActivity(diplotype));
      }
      else {
        diplotype.addPhenotype(gp.getPhenotypeForDiplotype(diplotype));
      }
    }
    diplotype.addLookupKey(gp.getLookupKeyForDiplotype(diplotype));
    gp.assignActivity(diplotype.getAllele1());
    gp.assignActivity(diplotype.getAllele2());
  }

  /**
   * Make a list of phenotype strings based on what diplotypes are specified.
   *
   * <em>Note:</em> This only applies when diplotypes are available in the {@link Diplotype} object. If only phenotype is supplied
   * then this will return an empty list.
   * @param diplotype the diplotype to analyze
   * @return a List of all applicable phenotype strings for this diplotype
   */
  private static List<String> makeHlaPhenotype(Diplotype diplotype) {
    if (!diplotype.getGene().equals("HLA-A") && !diplotype.getGene().equals("HLA-B")) {
      throw new RuntimeException("Gene not supported for HLA phenotype calling: " + diplotype.getGene());
    }
    if (diplotype.isUnknownAlleles()) {
      return new ArrayList<>();
    }
    if (diplotype.getGene().equals("HLA-A")) {
      return ImmutableList.of(diplotype.containsAllele("*31:01"));
    } else {
      List<String> phenotypes = new ArrayList<>();
      phenotypes.add(diplotype.containsAllele("*15:02"));
      phenotypes.add(diplotype.containsAllele("*57:01"));
      phenotypes.add(diplotype.containsAllele("*58:01"));
      return phenotypes;
    }
  }
}
