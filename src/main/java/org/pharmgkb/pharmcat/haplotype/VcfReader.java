package org.pharmgkb.pharmcat.haplotype;

import java.io.BufferedReader;
import java.io.IOException;
import java.lang.invoke.MethodHandles;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.List;
import java.util.SortedMap;
import java.util.TreeMap;
import java.util.regex.Pattern;
import javax.annotation.Nonnull;
import javax.annotation.Nullable;
import javax.annotation.concurrent.ThreadSafe;
import com.google.common.base.Preconditions;
import com.google.common.collect.ImmutableSet;
import org.pharmgkb.common.comparator.ChromosomePositionComparator;
import org.pharmgkb.parser.vcf.VcfLineParser;
import org.pharmgkb.parser.vcf.VcfParser;
import org.pharmgkb.parser.vcf.model.VcfMetadata;
import org.pharmgkb.parser.vcf.model.VcfPosition;
import org.pharmgkb.parser.vcf.model.VcfSample;
import org.pharmgkb.pharmcat.ParseException;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;


/**
 * This class reads VCF files and pulls the sample's alleles for positions of interest (i.e. is necessary to make a
 * haplotype call).
 *
 * @author Mark Woon
 */
@ThreadSafe
public class VcfReader {
  private static final Logger sf_logger = LoggerFactory.getLogger(MethodHandles.lookup().lookupClass());
  private static final Pattern sf_gtDelimiter = Pattern.compile("[\\|/]");
  private static final Pattern sf_noCallPattern = Pattern.compile("^[\\.\\|/]+$");
  private static final Pattern sf_allelePattern = Pattern.compile("^[AaCcGgTt]+$");
  private ImmutableSet<String> m_locationsOfInterest;


  /**
   * Constructor.
   *
   * @param locationsOfInterest set of chr:positions to pull alleles for
   */
  public VcfReader(ImmutableSet<String> locationsOfInterest) {
    m_locationsOfInterest = locationsOfInterest;
  }


  /**
   * Read VCF file and pull the sample's alleles for positions of interest.
   *
   * @return map of {@code <chr:position, SampleAllele>}
   */
  public SortedMap<String, SampleAllele> read(Path vcfFile) throws IOException {

    Preconditions.checkNotNull(vcfFile);
    Preconditions.checkArgument(Files.isRegularFile(vcfFile), "Not a file");
    Preconditions.checkArgument(Files.isReadable(vcfFile), "Not readable");
    Preconditions.checkArgument(vcfFile.toString().endsWith(".vcf"));

    // <chr:position, allele>
    try (BufferedReader reader = Files.newBufferedReader(vcfFile)) {
      // read VCF file
      SimpleLineParser lineParser = new SimpleLineParser();
      new VcfParser.Builder()
          .fromReader(reader)
          .parseWith(lineParser)
          .build().parse();
      return lineParser.getAlleleMap();
    }
  }


  private class SimpleLineParser implements VcfLineParser {
    // <chr:position, allele>
    private SortedMap<String, SampleAllele> m_alleleMap = new TreeMap<>(ChromosomePositionComparator.getComparator());


    public SortedMap<String, SampleAllele> getAlleleMap() {
      return m_alleleMap;
    }


    @Override
    public void parseLine(VcfMetadata metadata, VcfPosition position, List<VcfSample> sampleData) {

      String chrPos = position.getChromosome() + ":" + position.getPosition();

      if (sampleData.isEmpty()) {
        sf_logger.warn("Missing sample data on {}", chrPos);
        return;
      }

      if (!m_locationsOfInterest.contains(chrPos)) {
        sf_logger.warn("Ignoring {}", chrPos);
        return;
      }
      if (m_alleleMap.containsKey(chrPos)) {
        sf_logger.warn("Ignoring duplicate entry for {}", chrPos);
        return;
      }

      if (sampleData.size() > 1) {
        sf_logger.warn("Multiple samples found, only using first");
      }

      String gt = sampleData.get(0).getProperty("GT");
      if (gt == null) {
        sf_logger.warn("No genotype for {}", chrPos);
        return;
      }
      if (sf_noCallPattern.matcher(gt).matches()) {
        sf_logger.info("No call ({}) at {}", gt, chrPos);
        return;
      }

      int[] alleleIdxs = sf_gtDelimiter.splitAsStream(gt)
          .mapToInt(Integer::parseInt)
          .toArray();
      // normalize alleles to use same syntax as haplotype definition
      List<String> alleles = new ArrayList<>();
      if (position.getAltBases().size() == 0) {
        String gt1 = position.getAllele(0);
        validateAlleles(chrPos, gt1, null);

        String g[] = normalizeAlleles(gt1, null);
        alleles.add(g[0]);

      } else {
        for (String gt2 : position.getAltBases()) {
          String gt1 = position.getAllele(0);
          validateAlleles(chrPos, gt1, gt2);

          String g[] = normalizeAlleles(gt1, gt2);
          String g1 = g[0];
          String g2 = g[1];
          if (alleles.size() == 0) {
            alleles.add(g1);
            alleles.add(g2);
          } else {
            if (!alleles.get(0).equals(g1)) {
              throw new IllegalArgumentException("Getting different reference allele for " + chrPos);
            }
            alleles.add(g2);
          }
        }
      }

      String a1 = alleles.get(alleleIdxs[0]);
      String a2 = null;
      if (alleleIdxs.length > 1) {
        a2 = alleles.get(alleleIdxs[1]);
      } else {
        sf_logger.warn("Only a single allele for {}", chrPos);
      }

      // genotype divided by "|" if phased and "/" if unphased
      boolean isPhased = true;
      if (gt.contains("/") && a2 != null && !a1.equalsIgnoreCase(a2)) {
        isPhased = false;
      }

      List<String> vcfAlleles = new ArrayList<>();
      vcfAlleles.add(position.getRef());
      vcfAlleles.addAll(position.getAltBases());

      m_alleleMap.put(chrPos, new SampleAllele(position.getChromosome(), position.getPosition(), a1, a2, isPhased, vcfAlleles));
    }


    /**
     * Normalize alleles from VCF to match syntax from allele definitions
     */
    private String[] normalizeAlleles(@Nonnull String refAllele, @Nullable String varAllele) {

      refAllele = refAllele.toUpperCase();

      if (varAllele != null) {
        varAllele = varAllele.toUpperCase();
        return new String[] { refAllele, varAllele };
      } else {
        return new String[] { refAllele };
      }
    }
  }


  /**
   * Validate GT input per VCF 4.2 specification.
   */
  private static void validateAlleles(@Nonnull String chrPos, @Nonnull String gt1, @Nullable String gt2) {

    StringBuilder problems = new StringBuilder();
    if (gt1.startsWith("<")) {
      problems.append("Don't know how to handle ref structural variant '")
          .append(gt1)
          .append("'");
    } else if (gt1.toUpperCase().contains("N")) {
      problems.append("Don't know how to handle ambiguous allele in ref '")
          .append(gt1)
          .append("'");
    } else if (!sf_allelePattern.matcher(gt1).matches()) {
      problems.append("Unsupported bases in ref '")
          .append(gt1)
          .append("'");
    }

    if (gt2 != null) {
      if (gt2.startsWith("<")) {
        if (problems.length() > 0) {
          problems.append(System.lineSeparator());
        }
        problems.append("Don't know how to handle alt structural variant '")
            .append(gt2)
            .append("'");
      } else {
        if (!sf_allelePattern.matcher(gt1).matches()) {
          if (problems.length() > 0) {
            problems.append(System.lineSeparator());
          }
          if (gt2.toUpperCase().contains("N")) {
            problems.append("Don't know how to handle ambiguous allele in alt '")
                .append(gt2)
                .append("'");
          } else if (gt2.contains("*")) {
            problems.append("Don't know how to handle missing allele in alt '")
                .append(gt2)
                .append("'");
          } else {
            problems.append("Unsupported bases in alt '")
                .append(gt1)
                .append("'");
          }
        }
      }
    }

    if (problems.length() > 0) {
      throw new ParseException("Problem at " + chrPos + ":" + System.lineSeparator() + problems.toString());
    }
  }
}
