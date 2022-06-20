package org.pharmgkb.pharmcat.util;

import java.util.Arrays;
import java.util.List;
import java.util.Optional;
import java.util.stream.Collectors;
import com.google.common.base.Preconditions;
import com.google.common.collect.ImmutableSet;
import org.apache.commons.lang3.StringUtils;
import org.pharmgkb.pharmcat.ParseException;
import org.pharmgkb.pharmcat.UnexpectedStateException;
import org.pharmgkb.pharmcat.reporter.DiplotypeFactory;
import org.pharmgkb.pharmcat.reporter.model.VariantReport;
import org.pharmgkb.pharmcat.reporter.model.result.Diplotype;
import org.pharmgkb.pharmcat.reporter.model.result.GeneReport;


/**
 * This class calls SLCO1B1 using only rs4149056. This is only to be used when the matcher can't make a call
 *
 * @author Ryan Whaley
 */
public class Slco1b1AlleleMatcher {
  
  private static final String SLCO1B1 = "SLCO1B1";
  private static final String sf_callPosition = "rs4149056";
  private static final String sf_vcfCallSplitter = "[|/]";

  /**
   * Returns true if this matcher is appropriate to use for the given {@link GeneReport}
   * @param report a GeneReport to check
   * @return true if this matcher should be used, false otherwise
   */
  public static boolean shouldBeUsedOn(GeneReport report) {
    return report != null
        && report.getGene().equals(SLCO1B1)
        && (report.getMatcherDiplotypes().isEmpty() || report.getMatcherDiplotypes().stream().allMatch(Diplotype::isUnknown));
  }

  /**
   * Makes the diplotype call that should be used for lookup of guideline annotation groups
   * @param report an SLCO1B1 {@link GeneReport}
   * @param diplotypeFactory The factory class responsible for constructing diplotypes
   * @return an Optional Diplotype result, can be empty if the necessary position is missing
   */
  public static Optional<Diplotype> makeLookupCalls(GeneReport report, DiplotypeFactory diplotypeFactory) {
    Preconditions.checkNotNull(report);
    Preconditions.checkArgument(report.getGene().equals(SLCO1B1), "Can only be used on SLCO1B1");
    
    List<VariantReport> variants = report.getVariantReports().stream()
        .filter(v -> v.getDbSnpId() != null && v.getDbSnpId().equals(sf_callPosition)).toList();
    if (variants.size() == 0) return Optional.empty();
    if (variants.size() > 1) throw new UnexpectedStateException("More than one report found for " + sf_callPosition);

    VariantReport variant = variants.get(0);
    if (StringUtils.isBlank(variant.getCall())) return Optional.empty();
    
    String diplotypeText = makeDiplotype(variant);
    Diplotype diplotype = diplotypeFactory.makeDiplotypes(ImmutableSet.of(diplotypeText)).get(0);
    
    diplotype.setVariant(variant);
    
    return Optional.of(diplotype);
  }

  /**
   * Make a diplotype string based off of the calls for the given variant
   * @param variant a {@link VariantReport} for rs4149056
   * @return a String diplotype in the form "*1/*5"
   */
  private static String makeDiplotype(VariantReport variant) {
    String[] alleles = variant.getCall().split(sf_vcfCallSplitter);

    if (Arrays.equals(alleles, new String[]{"T","T"})) {
      return "*1/*1";
    }
    else if (Arrays.equals(alleles, new String[]{"C","C"})) {
      return "*5/*5";
    }
    else if ((Arrays.equals(alleles, new String[]{"T","C"})) || (Arrays.equals(alleles, new String[]{"C","T"}))) {
      return "*1/*5";
    }
    else {
      throw new ParseException("Unexpected genotype for " + variant);
    }

  }
}
