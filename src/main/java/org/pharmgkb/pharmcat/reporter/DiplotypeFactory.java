package org.pharmgkb.pharmcat.reporter;

import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.SortedSet;
import java.util.TreeSet;
import java.util.stream.Collectors;
import com.google.common.base.Preconditions;
import com.google.common.collect.ImmutableList;
import com.google.common.collect.ImmutableSet;
import com.google.common.collect.Lists;
import org.apache.commons.lang3.StringUtils;
import org.pharmgkb.pharmcat.Env;
import org.pharmgkb.pharmcat.definition.model.GenePhenotype;
import org.pharmgkb.pharmcat.definition.model.NamedAllele;
import org.pharmgkb.pharmcat.haplotype.model.BaseMatch;
import org.pharmgkb.pharmcat.haplotype.model.CombinationMatch;
import org.pharmgkb.pharmcat.haplotype.model.DiplotypeMatch;
import org.pharmgkb.pharmcat.haplotype.model.GeneCall;
import org.pharmgkb.pharmcat.reporter.model.DataSource;
import org.pharmgkb.pharmcat.reporter.model.OutsideCall;
import org.pharmgkb.pharmcat.reporter.model.result.Diplotype;
import org.pharmgkb.pharmcat.reporter.model.result.GeneReport;
import org.pharmgkb.pharmcat.reporter.model.result.Haplotype;
import org.pharmgkb.pharmcat.reporter.model.result.Observation;
import org.pharmgkb.pharmcat.util.HaplotypeActivityComparator;

import static org.pharmgkb.pharmcat.reporter.model.result.GeneReport.isSinglePloidy;


/**
 * Factory class for spitting out {@link Diplotype} objects from known information about the gene
 * <p>
 * This is designed to take information from either the Matcher or from outside calls and produce consistent results
 *
 * @author Ryan Whaley
 */
public class DiplotypeFactory {
  public static final String UNASSIGNED_FUNCTION = "Unassigned function";
  private static final Set<String> PHENOTYPE_ONLY = ImmutableSet.of("HLA-A", "HLA-B");
  private final Map<DataSource, Map<String, Haplotype>> m_haplotypeCache = new HashMap<>();

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
      dips.add(new Diplotype(m_gene, makeHaplotype(na.getName(), source)));
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
    Diplotype diplotype = new Diplotype(m_gene, makeHaplotype(baseMatch.getName(), source));
    fillDiplotype(diplotype, m_env, source);
    return diplotype;
    }

  private Diplotype makeDiplotype(DiplotypeMatch diplotypeMatch, DataSource source) {
    BaseMatch h1 = diplotypeMatch.getHaplotype1();
    BaseMatch h2 = diplotypeMatch.getHaplotype2();

    Diplotype diplotype = new Diplotype(m_gene, makeHaplotype(h1.getName(), source),
        makeHaplotype(h2.getName(), source));
    fillDiplotype(diplotype, m_env, source);

    diplotype.setCombination(h1 instanceof CombinationMatch || h2 instanceof CombinationMatch);
    return diplotype;
  }



  public List<Diplotype> makeLeastFunctionDiplotypes(@SuppressWarnings("rawtypes") Collection matches,
      DataSource source, boolean hasTrueDiplotype) {

    if (matches.size() > 0) {
      List<Diplotype> diplotypes = new ArrayList<>();
      if (hasTrueDiplotype) {
        for (Object obj : matches) {
          List<String> hapNames = new ArrayList<>();
          if (obj instanceof BaseMatch bm) {
            hapNames.addAll(getHaps(bm));
          } else if (obj instanceof DiplotypeMatch dm) {
            hapNames.addAll(getHaps(dm.getHaplotype1()));
            hapNames.addAll(getHaps(dm.getHaplotype2()));
          } else if (obj instanceof String text) {
            String[] haplotypes = DiplotypeFactory.splitDiplotype(m_gene, text);
            for (String hap : haplotypes) {
              hapNames.addAll(DiplotypeFactory.splitHaplotype(hap));
            }
          }
          diplotypes.add(makeLeastFunctionDiplotype(hapNames, source));
        }
      } else {
        List<String> hapNames = new ArrayList<>();
        for (Object obj : matches) {
          hapNames.addAll(getHaps((BaseMatch)obj));
        }
        diplotypes.add(makeLeastFunctionDiplotype(hapNames, source));
      }
      return diplotypes;
    } else {
      return ImmutableList.of(makeUnknownDiplotype(m_gene, m_env, source));
    }
  }

  private Diplotype makeLeastFunctionDiplotype(List<String> hapNames, DataSource source) {
    List<Haplotype> haplotypes = hapNames.stream()
        .map(h -> makeHaplotype(h, source))
        .sorted(HaplotypeActivityComparator.getComparator())
        .toList();
    Haplotype hap1 = haplotypes.get(0);
    Haplotype hap2 = null;
    if (haplotypes.size() > 1) {
      hap2 = haplotypes.get(1);
    }
    Diplotype diplotype = new Diplotype(m_gene, hap1, hap2);
    fillDiplotype(diplotype, m_env, source);
    if (haplotypes.size() > 2) {
      diplotype.setObserved(Observation.INFERRED);
    }
    return diplotype;
  }



  /**
   * Make diplotype objects based on {@link OutsideCall} objects
   */
  public List<Diplotype> makeDiplotypes(OutsideCall outsideCall, DataSource source) {
    Preconditions.checkNotNull(outsideCall);
    Preconditions.checkArgument(outsideCall.getGene().equals(m_gene));

    if (outsideCall.getDiplotypes().size() > 0) {
      return makeDiplotypes(outsideCall.getDiplotypes(), source);
    } else {
      // phenotype-based diplotype
      return ImmutableList.of(new Diplotype(m_gene, outsideCall.getPhenotype()));
    }
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
    Haplotype hap1 = makeHaplotype(alleles[0], source);
    Haplotype hap2 = alleles.length == 2 ? makeHaplotype(alleles[1], source) : null;
    Diplotype diplotype = new Diplotype(m_gene, hap1, hap2);
    fillDiplotype(diplotype, m_env, source);
    return diplotype;
  }

  private boolean isPhenotypeOnly() {
    return isPhenotypeOnly(m_gene);
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

  private static void fillDiplotype(Diplotype diplotype, Env env, DataSource source) {
    if (diplotype.getGene().startsWith("HLA")) {
      diplotype.setPhenotypes(makeHlaPhenotype(diplotype));
      diplotype.setLookupKeys(makeHlaPhenotype(diplotype));
      return;
    }
    GenePhenotype gp = env.getPhenotype(diplotype.getGene(), source);
    if (gp == null) {
      return;
    }
    diplotype.addPhenotype(gp.getPhenotypeForDiplotype(diplotype));
    diplotype.addLookupKey(gp.getLookupKeyForDiplotype(diplotype));
    gp.assignActivity(diplotype.getAllele1());
    gp.assignActivity(diplotype.getAllele2());
    diplotype.calculateActivityScore();
  }

  private static List<String> makeHlaPhenotype(Diplotype diplotype) {
    if (!diplotype.getGene().equals("HLA-A") && !diplotype.getGene().equals("HLA-B")) {
      throw new RuntimeException("Gene not supported for HLA phenotype calling: " + diplotype.getGene());
    }
    if (diplotype.isUnknown()) {
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


  /**
   * Make or retrieve a cached Haplotype object that corresponds to the given allele name.
   * @param name an allele name (e.g. *20)
   * @return a Haplotype object for that allele (new or cached)
   */
  public Haplotype makeHaplotype(String name, DataSource source) {
    // return cache value if possible
    if (m_haplotypeCache.containsKey(source) && m_haplotypeCache.get(source).containsKey(name)) {
      return m_haplotypeCache.get(source).get(name);
    }

    Haplotype haplotype = new Haplotype(m_gene, name);
    GenePhenotype gp = m_env.getPhenotype(m_gene, source);
    if (gp != null) {
      haplotype.setFunction(gp.findHaplotypeFunction(name).orElse(DiplotypeFactory.UNASSIGNED_FUNCTION));
      haplotype.setActivityValue(gp.findHaplotypeActivity(name).orElse(null));
  }
    haplotype.setReference(name.equals(m_env.getReferenceAllele(m_gene)));

    m_haplotypeCache.computeIfAbsent(source, (s) -> new HashMap<>())
        .put(name, haplotype);
    return haplotype;
  }


  private static List<String> getHaps(BaseMatch baseMatch) {
    if (baseMatch instanceof CombinationMatch cm) {
      return cm.getComponentHaplotypes().stream()
          .map(NamedAllele::getName)
          .collect(Collectors.toList());
    } else {
      return Lists.newArrayList(baseMatch.getHaplotype().getName());
    }
  }
}
