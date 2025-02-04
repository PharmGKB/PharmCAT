package org.pharmgkb.pharmcat.reporter.caller;

import java.util.Arrays;
import java.util.Comparator;
import java.util.List;
import java.util.Objects;
import java.util.SortedSet;
import java.util.TreeSet;
import com.google.common.base.Preconditions;
import org.apache.commons.lang3.StringUtils;
import org.checkerframework.checker.nullness.qual.Nullable;
import org.pharmgkb.pharmcat.Env;
import org.pharmgkb.pharmcat.UnexpectedStateException;
import org.pharmgkb.pharmcat.reporter.model.DataSource;
import org.pharmgkb.pharmcat.reporter.model.VariantReport;
import org.pharmgkb.pharmcat.reporter.model.result.Diplotype;
import org.pharmgkb.pharmcat.reporter.model.result.GeneReport;


/**
 * This class calls SLCO1B1 using only rs4149056. This is only to be used when the matcher can't make a call.
 *
 * @author Ryan Whaley
 */
public class Slco1b1CustomCaller {
  private static final String GENE = "SLCO1B1";
  private static final String sf_callPosition = "rs4149056";
  private static final String sf_vcfCallSplitter = "[|/]";


  public static boolean isSlco1b1(String gene) {
    return GENE.equals(gene);
  }

  public static boolean isSlco1b1(GeneReport geneReport) {
    return geneReport.getGene().equals(GENE);
  }


  /**
   * Gets the diplotype calls that should be used to look up guideline annotations.
   */
  public static SortedSet<Diplotype> inferDiplotypes(GeneReport report, Env env,
      DataSource source) {
    Preconditions.checkArgument(report.getGene().equals(GENE), "Can only be used on SLCO1B1");

    if (report.isOutsideCall() || !report.getSourceDiplotypes().stream().allMatch(Diplotype::isUnknown)) {
      return report.getSourceDiplotypes();
    }

    List<VariantReport> variants = report.getVariantReports().stream()
        .filter(v -> sf_callPosition.equals(v.getDbSnpId()))
        .toList();
    if (variants.isEmpty()) {
      return report.getSourceDiplotypes();
    }
    if (variants.size() > 1) {
      throw new UnexpectedStateException("More than one report found for " + sf_callPosition);
    }

    VariantReport variant = variants.get(0);
    String[] haps = convertToHaplotypes(variant);
    if (haps == null) {
      return report.getSourceDiplotypes();
    }

    Diplotype diplotype = new Diplotype(GENE, haps[0], haps[1], env, source, 0);
    diplotype.setVariant(variant);
    diplotype.setInferred(true);

    String[] alleles = Objects.requireNonNull(splitVariant(variant));
    diplotype.setInferredSourceDiplotype(new Diplotype(GENE,
        sf_callPosition + " " + alleles[0], sf_callPosition + " " + alleles[1], env, source, 0));

    SortedSet<Diplotype> dips = new TreeSet<>();
    dips.add(diplotype);

    // add variant info to source diplotypes
    report.getSourceDiplotypes().first().setVariant(variant);
    return dips;
  }

  private static @Nullable String[] splitVariant(VariantReport variant) {
    String[] alleles = variant.getCall().split(sf_vcfCallSplitter);
    if (alleles.length != 2) {
      return null;
    }
    // sort descending (T before C so *1 before *5)
    Arrays.sort(alleles, Comparator.reverseOrder());
    return alleles;
  }

  /**
   * Make a diplotype string based off of the calls for the given variant.
   */
  private static @Nullable String[] convertToHaplotypes(VariantReport variant) {
    if (StringUtils.isBlank(variant.getCall())) {
      return null;
    }
    String[] alleles = splitVariant(variant);
    String[] haps = new String[2];
    haps[0] = alleleToHap(alleles[0]);
    if (haps[0] == null) {
      return null;
    }
    haps[1] = alleleToHap(alleles[1]);
    if (haps[1] == null) {
      return null;
    }
    return haps;
  }

  private static String alleleToHap(String allele) {
    return switch (allele) {
      case "T" -> "*1";
      case "C" -> "*5";
      default -> null;
    };
  }
}
