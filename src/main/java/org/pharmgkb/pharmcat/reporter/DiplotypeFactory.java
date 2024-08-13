package org.pharmgkb.pharmcat.reporter;

import java.util.ArrayList;
import java.util.Collection;
import java.util.List;
import java.util.Set;
import java.util.SortedSet;
import java.util.TreeSet;
import com.google.common.collect.ImmutableList;
import com.google.common.collect.ImmutableSet;
import com.google.common.collect.Lists;
import org.pharmgkb.pharmcat.Env;
import org.pharmgkb.pharmcat.definition.model.NamedAllele;
import org.pharmgkb.pharmcat.haplotype.model.BaseMatch;
import org.pharmgkb.pharmcat.haplotype.model.CombinationMatch;
import org.pharmgkb.pharmcat.haplotype.model.DiplotypeMatch;
import org.pharmgkb.pharmcat.haplotype.model.GeneCall;
import org.pharmgkb.pharmcat.haplotype.model.HaplotypeMatch;
import org.pharmgkb.pharmcat.reporter.model.DataSource;
import org.pharmgkb.pharmcat.reporter.model.result.Diplotype;
import org.pharmgkb.pharmcat.reporter.model.result.GeneReport;
import org.pharmgkb.pharmcat.reporter.model.result.Haplotype;

import static org.pharmgkb.pharmcat.Constants.isSinglePloidy;


/**
 * Factory class for spitting out {@link Diplotype} objects from known information about the gene
 * <p>
 * This is designed to take information from either the Matcher or from outside calls and produce consistent results
 *
 * @author Ryan Whaley
 */
public class DiplotypeFactory {
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
      dips.add(new Diplotype(m_gene, na.getName(), null, m_env, source));
    }
    return dips;
  }


  public List<Diplotype> makeDiplotypes(Collection<DiplotypeMatch> matches, DataSource source) {
    if (matches.isEmpty()) {
      return ImmutableList.of(makeUnknownDiplotype(m_gene, m_env, source));
    }
    return matches.stream()
        .map((dm) -> {
          BaseMatch h1 = dm.getHaplotype1();
          BaseMatch h2 = dm.getHaplotype2();
          String h2Name = null;
          if (h2 != null) {
            h2Name = h2.getName();
          }

          Diplotype diplotype = new Diplotype(m_gene, h1.getName(), h2Name, m_env, source);
          diplotype.setCombination(h1 instanceof CombinationMatch || h2 instanceof CombinationMatch);
          return diplotype;
        })
        .toList();
  }

  public List<Diplotype> makeDiplotypesFromHaplotypeMatches(Collection<HaplotypeMatch> matches, DataSource source) {
    if (matches.isEmpty()) {
      return ImmutableList.of(makeUnknownDiplotype(m_gene, m_env, source));
    }
    return matches.stream()
        .map((hm) -> new Diplotype(m_gene, hm.getName(), null, m_env, source))
        .toList();
  }



  public static String[] splitDiplotype(String gene, String diplotypeText) {
    if (diplotypeText.contains("/")) {
      if (isSinglePloidy(gene) && !GeneReport.isXChromo(gene)) {
        throw new BadOutsideCallException("Cannot specify two genotypes [" + diplotypeText +
            "] for single chromosome gene " + gene);
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
    if (CombinationMatch.isCombinationName(haplotypeText))  {
      return CombinationMatch.COMBINATION_NAME_SPLITTER.splitToList(CombinationMatch.extractCombinationName(haplotypeText));
    }
    return Lists.newArrayList(haplotypeText);
  }


  private static boolean isPhenotypeOnly(String gene) {
    return PHENOTYPE_ONLY.contains(gene);
  }


  public static Diplotype makeUnknownDiplotype(String gene, Env env, DataSource source) {
    Haplotype haplotype = new Haplotype(gene, Haplotype.UNKNOWN);
    Diplotype diplotype;
    if (isSinglePloidy(gene)) {
      diplotype = new Diplotype(gene, haplotype, null, env, source);
    } else {
      diplotype = new Diplotype(gene, haplotype, haplotype, env, source);
    }
    return diplotype;
  }

  /**
   * Make a list of phenotype strings based on what diplotypes are specified.
   * <p>
   * <em>Note:</em> This only applies when diplotypes are available in the {@link Diplotype} object.
   * If only phenotype is supplied, then this will return an empty list.
   *
   * @param diplotype the diplotype to analyze
   * @return a List of all applicable phenotype strings for this diplotype
   */
  public static List<String> makeHlaPhenotype(Diplotype diplotype) {
    if (!diplotype.getGene().equals("HLA-A") && !diplotype.getGene().equals("HLA-B")) {
      throw new RuntimeException("Gene not supported for HLA phenotype calling: " + diplotype.getGene());
    }
    if (diplotype.isUnknownAlleles()) {
      return new ArrayList<>();
    }
    if (diplotype.getGene().equals("HLA-A")) {
      return ImmutableList.of(hlaPhenotype(diplotype, "*31:01"));
    } else {
      List<String> phenotypes = new ArrayList<>();
      phenotypes.add(hlaPhenotype(diplotype, "*15:02"));
      phenotypes.add(hlaPhenotype(diplotype, "*57:01"));
      phenotypes.add(hlaPhenotype(diplotype, "*58:01"));
      return phenotypes;
    }
  }

  private static String hlaPhenotype(Diplotype diplotype, String allele) {
    if (nameContainsAllele(diplotype.getAllele1(), allele) || nameContainsAllele(diplotype.getAllele2(), allele)) {
      return allele + " positive";
    } else {
      return allele + " negative";
    }
  }

  private static boolean nameContainsAllele(Haplotype hap, String allele) {
    return hap != null && hap.getName().contains(allele);
  }
}
