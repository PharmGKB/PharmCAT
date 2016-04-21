package org.pharmgkb.pharmcat.haplotype;

import java.io.BufferedReader;
import java.io.IOException;
import java.lang.invoke.MethodHandles;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.List;
import java.util.Set;
import java.util.SortedMap;
import java.util.TreeMap;
import java.util.regex.Pattern;
import javax.annotation.Nonnull;
import com.google.common.base.Preconditions;
import org.pharmgkb.common.comparator.ChromosomePositionComparator;
import org.pharmgkb.parser.vcf.VcfLineParser;
import org.pharmgkb.parser.vcf.VcfParser;
import org.pharmgkb.parser.vcf.model.VcfMetadata;
import org.pharmgkb.parser.vcf.model.VcfPosition;
import org.pharmgkb.parser.vcf.model.VcfSample;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;


/**
 * @author Mark Woon
 */
public class VcfReader {
  private static final Logger sf_logger = LoggerFactory.getLogger(MethodHandles.lookup().lookupClass());
  private static final Pattern sf_gtDelimiter = Pattern.compile("[\\|/]");
  private Set<String> m_locationsOfInterest;


  /**
   * Constructor.
   *
   * @param locationsOfInterest set of chr:positions to pull alleles for
   */
  public VcfReader(Set<String> locationsOfInterest) {
    m_locationsOfInterest = locationsOfInterest;
  }


  public SortedMap<String, SampleAllele>  read(Path vcfFile) throws IOException {

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
      if (gt.equals(".") || gt.equals("./.")) {
        return;
      }

      int[] alleleIdxs = sf_gtDelimiter.splitAsStream(gt)
          .mapToInt(Integer::parseInt)
          .toArray();
      // normalize alleles to use same syntax as haplotype definition
      List<String> alleles = new ArrayList<>();
      if (position.getAltBases().size() == 0) {
        alleles.add(position.getAllele(0).toUpperCase());
        sf_logger.warn("Unexpected state: no alt alleles for {}",  chrPos);
      } else {
        for (int y = 0; y < position.getAltBases().size(); y += 1) {
          String g[] = normalizeAlleles(position.getAllele(0).toUpperCase(), position.getAltBases().get(y).toUpperCase());
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
      m_alleleMap.put(chrPos, new SampleAllele(position.getChromosome(), position.getPosition(), a1, a2, isPhased));
    }

    private String[] normalizeAlleles(@Nonnull String refAllele, @Nonnull String varAllele) {
      return new String[] { refAllele, varAllele };
    }
  }
}
