package org.pharmgkb.pharmcat.reporter.model.result;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.List;
import java.util.Set;
import java.util.SortedSet;
import java.util.TreeSet;
import javax.annotation.Nonnull;
import com.google.common.collect.ImmutableList;
import org.pharmgkb.common.comparator.HaplotypeNameComparator;
import org.pharmgkb.pharmcat.haplotype.model.GeneCall;
import org.pharmgkb.pharmcat.haplotype.model.Variant;


/**
 * Gene report data specifically for the UGT1A1 gene. This gene requires specific, non-standard rules for calling it so
 * override all the necessary methods.
 *
 * @author Ryan Whaley
 */
public class GeneReportUgt1a1 extends GeneReport {

  private static final List<String> sf_groupA = ImmutableList.of("*28", "*80");
  private static final List<String> sf_groupB = ImmutableList.of("*6", "*27", "*37");
  private static final List<String> sf_groupAB = new ArrayList<>();
  static {
    sf_groupAB.addAll(sf_groupA);
    sf_groupAB.addAll(sf_groupB);
  }

  private List<String> m_haplotypes = new ArrayList<>();
  private Set<String> m_displayFunctions = new TreeSet<>();

  /**
   * public constructor
   */
  public GeneReportUgt1a1(@Nonnull String geneSymbol) {
    super(geneSymbol);
  }

  /**
   * Sets internal data for
   * @param call a {@link GeneCall} that has been made by the NamedAlleleMatcher
   */
  @Override
  public void setCallData(@Nonnull GeneCall call) {
    super.setCallData(call);

    for (Variant variant : call.getVariants()) {
      String[] alleles = variant.getVcfCall().split("[|/]");

      if (variant.getPosition() == 233757013) {
        Arrays.stream(alleles)
            .filter(a -> a.equals("G"))
            .forEach(a -> m_haplotypes.add("*60"));
      }

      if (variant.getPosition() == 233759924) {
        Arrays.stream(alleles)
            .filter(a -> a.equals("T"))
            .forEach(a -> m_haplotypes.add("*80"));
      }

      if (variant.getPosition() == 233760498) {
        Arrays.stream(alleles)
            .filter(a -> a.equals("A"))
            .forEach(a -> m_haplotypes.add("*6"));
      }

      if (variant.getPosition() == 233760973) {
        Arrays.stream(alleles)
            .filter(a -> a.equals("A"))
            .forEach(a -> m_haplotypes.add("*27"));
      }

      if (variant.getPosition() == 233760233) {
        Arrays.stream(alleles)
            .filter(a -> a.equals("CATAT"))
            .forEach(a -> m_haplotypes.add("*28"));
        Arrays.stream(alleles)
            .filter(a -> a.equals("C"))
            .forEach(a -> m_haplotypes.add("*36"));
        Arrays.stream(alleles)
            .filter(a -> a.equals("CATATAT"))
            .forEach(a -> m_haplotypes.add("*37"));
      }
    }
  }

  public List<String> makeDiplotypeCalls() {
    long countHomoAB = sf_groupAB.stream()
        .filter(a -> m_haplotypes.stream().filter(h -> h.equals(a)).count() >= 2)
        .count();
    long countHetA = sf_groupA.stream()
        .filter(a -> m_haplotypes.stream().filter(h -> h.equals(a)).count() >= 1)
        .count();
    long countHetB = sf_groupB.stream()
        .filter(a -> m_haplotypes.stream().filter(h -> h.equals(a)).count() >= 1)
        .count();

    if (countHomoAB >= 2 || (countHetA >= 1 && countHetB >= 1) || countHetB >= 2) {
      return ImmutableList.of("*80/*80");
    }
    else if (countHetA == 0 && countHetB == 0) {
      return ImmutableList.of("*1/*1");
    }
    else {
      return ImmutableList.of("*1/*80");
    }
  }

  @Override
  public Collection<String> printDisplayCalls() {
    Collection<String> alleles = getDistinctAlleles();
    List<String> displayCalls = new ArrayList<>();

    for (String allele : alleles) {
      long count = m_haplotypes.stream().filter(h -> h.equals(allele)).count();
      if (count >= 2) {
        displayCalls.add(allele + " (homozygous)");
      } else {
        displayCalls.add(allele + " (heterozygous)");
      }
    }

    return displayCalls;
  }

  public void addDisplayFunction(String function) {
    if (function != null) {
      m_displayFunctions.add(function);
    }
  }

  @Override
  public Collection<String> printDisplayFunctions() {
    return m_displayFunctions;
  }

  public Collection<String> getDistinctAlleles() {
    SortedSet<String> alleles = new TreeSet<>(HaplotypeNameComparator.getComparator());
    alleles.addAll(m_haplotypes);
    return alleles;
  }
}
