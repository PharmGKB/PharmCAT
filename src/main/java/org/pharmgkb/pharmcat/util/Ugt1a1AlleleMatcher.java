package org.pharmgkb.pharmcat.util;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Set;
import javax.annotation.Nonnull;
import com.google.common.base.Preconditions;
import com.google.common.collect.ImmutableList;
import com.google.common.collect.ImmutableSet;
import org.pharmgkb.pharmcat.reporter.model.VariantReport;
import org.pharmgkb.pharmcat.reporter.model.result.GeneReport;


/**
 * This class will handle calling UGT1A1 diplotypes based on single variant position genotypes. This should only be used
 * in the case of a GeneCall that is either unphased or phased with multiple calls.
 *
 * @author Ryan Whaley
 */
public class Ugt1a1AlleleMatcher {

  private static final String sf_gtDelimiter = "[|/]";
  private static final List<String> sf_groupA = ImmutableList.of("*28", "*80");
  private static final List<String> sf_groupB = ImmutableList.of("*6", "*27", "*37");
  private static final List<String> sf_groupAB = new ArrayList<>();
  static {
    sf_groupAB.addAll(sf_groupA);
    sf_groupAB.addAll(sf_groupB);
  }

  /**
   * Returns true if this matcher is appropriate to use for the given {@link GeneReport}
   * @param gene a GeneReport to check
   * @return true if this matcher should be used, false otherwise
   */
  public static boolean shouldBeUsedOn(@Nonnull GeneReport gene) {
    return gene.getGene().equals("UGT1A1") && (!gene.isPhased() || gene.getMatcherDiplotypes().size() > 1);
  }

  /**
   * Makes diplotype calls (i.e. *1/*80) based on variant report information in the given {@link GeneReport}
   * @param report a Gene Reprot
   * @return a Set of diplotype calls (i.e. "*1/*80")
   */
  public static Set<String> makeLookupCalls(GeneReport report) {
    Preconditions.checkNotNull(report);
    Preconditions.checkArgument(report.getGene().equals("UGT1A1"), "Can only be used on UGT1A1");

    boolean geneMissing = report.getVariantReports().stream().allMatch(VariantReport::isMissing);
    if (geneMissing) {
      return ImmutableSet.of();
    }

    List<String> haplotypes = matchHaplotypes(report);

    long countHomoAB = sf_groupAB.stream()
        .filter(a -> haplotypes.stream().filter(h -> h.equals(a)).count() >= 2)
        .count();
    long countHetA = sf_groupA.stream()
        .filter(a -> haplotypes.stream().filter(h -> h.equals(a)).count() >= 1)
        .count();
    long countHetB = sf_groupB.stream()
        .filter(a -> haplotypes.stream().filter(h -> h.equals(a)).count() >= 1)
        .count();

    if (countHomoAB >= 1 || (countHetA >= 1 && countHetB >= 1) || countHetB >= 2) {
      return ImmutableSet.of("*80/*80");
    }
    else if (countHetA == 0 && countHetB == 0) {
      return ImmutableSet.of("*1/*1");
    }
    else {
      return ImmutableSet.of("*1/*80");
    }
  }

  /**
   * Generates a list of found Haplotype names in this {@link GeneReport} based on {@link VariantReport} data.
   * 
   * Basically, give a list of all found alleles based on specific positions. Most alleles are straight-forward but a 
   * special case exists for *80. It can be called based on two different positions with one (233760233) taking priority 
   * over the other (233759924).
   * 
   * @param report a {@link GeneReport} for the UGT1A1 gene
   * @return a List of String names for alleles found in this {@link GeneReport} (each allele can occur more than once)
   */
  private static List<String> matchHaplotypes(@Nonnull GeneReport report) {

    List<String> haplotypes = new ArrayList<>();

    VariantReport pos60233 = report.getVariantReports().stream()
        .filter(v -> v.getPosition() == 233760233).findFirst().orElseThrow(RuntimeException::new);
    VariantReport pos59924 = report.getVariantReports().stream()
        .filter(v -> v.getPosition() == 233759924).findFirst().orElseThrow(RuntimeException::new);
    
    if (!pos60233.isMissing()) {
      Arrays.stream(pos60233.getCall().split(sf_gtDelimiter))
          .filter(a -> a.startsWith("CATAT"))
          .forEach(a -> haplotypes.add("*80"));
    } else if (!pos59924.isMissing() && pos60233.isMissing()) {
      Arrays.stream(pos59924.getCall().split(sf_gtDelimiter))
          .filter(a -> a.equals("T"))
          .forEach(a -> haplotypes.add("*80"));
    }

    for (VariantReport variant : report.getVariantReports()) {
      if (variant.getCall() == null) {
        continue;
      }

      String[] alleles = variant.getCall().split(sf_gtDelimiter);

      if (variant.getPosition() == 233757013) {
        Arrays.stream(alleles)
            .filter(a -> a.equals("G"))
            .forEach(a -> haplotypes.add("*60"));
      }

      if (variant.getPosition() == 233760498) {
        Arrays.stream(alleles)
            .filter(a -> a.equals("A"))
            .forEach(a -> haplotypes.add("*6"));
      }

      if (variant.getPosition() == 233760973) {
        Arrays.stream(alleles)
            .filter(a -> a.equals("A"))
            .forEach(a -> haplotypes.add("*27"));
      }

      if (variant.getPosition() == 233760233) {
        Arrays.stream(alleles)
            .filter(a -> a.equals("CATAT"))
            .forEach(a -> haplotypes.add("*28"));
        Arrays.stream(alleles)
            .filter(a -> a.equals("C"))
            .forEach(a -> haplotypes.add("*36"));
        Arrays.stream(alleles)
            .filter(a -> a.equals("CATATAT"))
            .forEach(a -> haplotypes.add("*37"));
      }
    }
    return haplotypes;
  }
}
