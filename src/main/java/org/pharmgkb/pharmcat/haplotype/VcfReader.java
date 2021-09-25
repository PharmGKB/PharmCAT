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
import com.google.common.base.Preconditions;
import com.google.common.collect.ImmutableMap;
import com.google.common.collect.SortedSetMultimap;
import com.google.common.collect.TreeMultimap;
import org.checkerframework.checker.nullness.qual.Nullable;
import org.pharmgkb.common.comparator.ChromosomePositionComparator;
import org.pharmgkb.parser.vcf.VcfLineParser;
import org.pharmgkb.parser.vcf.VcfParser;
import org.pharmgkb.parser.vcf.model.ContigMetadata;
import org.pharmgkb.parser.vcf.model.VcfMetadata;
import org.pharmgkb.parser.vcf.model.VcfPosition;
import org.pharmgkb.parser.vcf.model.VcfSample;
import org.pharmgkb.pharmcat.ParseException;
import org.pharmgkb.pharmcat.definition.model.VariantLocus;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;


/**
 * This class reads VCF files and pulls the sample's alleles for positions of interest (i.e. is necessary to make a
 * haplotype call).
 *
 * @author Mark Woon
 */
public class VcfReader implements VcfLineParser {
  private static final Logger sf_logger = LoggerFactory.getLogger(MethodHandles.lookup().lookupClass());
  private static final Pattern sf_gtDelimiter = Pattern.compile("[|/]");
  private static final Pattern sf_noCallPattern = Pattern.compile("^[.|/]+$");
  private static final Pattern sf_allelePattern = Pattern.compile("^[AaCcGgTt]+$");
  private final ImmutableMap<String, VariantLocus> m_locationsOfInterest;
  private VcfMetadata m_vcfMetadata;
  private String m_genomeBuild;
  // <chr:position, allele>
  private final SortedMap<String, SampleAllele> m_alleleMap = new TreeMap<>(ChromosomePositionComparator.getComparator());
  // <chr:position, warning>
  private final SortedSetMultimap<String, String> m_warnings = TreeMultimap.create();


  /**
   * Constructor.
   * Reads in VCF file and pull the sample's alleles for positions of interest.
   *
   * @param locationsOfInterest set of chr:positions to pull alleles for
   */
  public VcfReader(ImmutableMap<String, VariantLocus> locationsOfInterest, Path vcfFile) throws IOException {
    m_locationsOfInterest = locationsOfInterest;
    read(vcfFile);
  }

  /**
   * Constructor.  Primarily for testing.
   * Reads in VCF file and pull the sample's alleles for <em>ALL</em> positions.
   */
  VcfReader(Path vcfFile) throws IOException {
    m_locationsOfInterest = null;
    read(vcfFile);
  }


  public VcfMetadata getVcfMetadata() {
    return m_vcfMetadata;
  }


  /**
   * Gets the genome build the VCF file is using.
   * This is pulled from the VCF contig assembly metadata.
   */
  public @Nullable String getGenomeBuild() {
    return m_genomeBuild;
  }


  /**
   * Gets sample data.
   *
   * @return map of {@code <chr:position, SampleAllele>}
   */
  public SortedMap<String, SampleAllele> getAlleleMap() {
    return m_alleleMap;
  }


  /**
   * Gets warnings from reading data, keyed to chromosomal position.
   */
  public SortedSetMultimap<String, String> getWarnings() {
    return m_warnings;
  }


  /**
   * Read VCF file.
   */
  private void read(Path vcfFile) throws IOException {

    Preconditions.checkNotNull(vcfFile);
    Preconditions.checkArgument(Files.isRegularFile(vcfFile), "%s is not a file", vcfFile);
    Preconditions.checkArgument(Files.isReadable(vcfFile), "%s is not readable", vcfFile);
    Preconditions.checkArgument(vcfFile.toString().endsWith(".vcf"), "%s is not a VCF file", vcfFile);

    // <chr:position, allele>
    try (BufferedReader reader = Files.newBufferedReader(vcfFile)) {
      // read VCF file
      VcfParser vcfParser = new VcfParser.Builder()
          .fromReader(reader)
          .parseWith(this)
          .build();
      m_vcfMetadata = vcfParser.parseMetadata();
      for (ContigMetadata cm : m_vcfMetadata.getContigs().values()) {
        if (cm.getAssembly() != null) {
          if (m_genomeBuild == null) {
            m_genomeBuild = cm.getAssembly();
          } else if (!m_genomeBuild.equals(cm.getAssembly())) {
            throw new IllegalStateException("VCF file uses different assemblies (" + m_genomeBuild + " vs " +
                cm.getAssembly() + " for contig)");
          }
        }
      }
      vcfParser.parse();
    }
  }

  private void addWarning(String chrPos, String msg) {
    m_warnings.put(chrPos, msg);
    sf_logger.warn(msg);
  }


  @Override
  public void parseLine(VcfMetadata metadata, VcfPosition position, List<VcfSample> sampleData) {

    String chrPos = position.getChromosome() + ":" + position.getPosition();

    if (sampleData.isEmpty()) {
      sf_logger.warn("Missing sample data on {}", chrPos);
      return;
    }

    VariantLocus varLoc = null;
    if (m_locationsOfInterest != null) {
      varLoc = m_locationsOfInterest.get(chrPos);
      if (varLoc == null) {
        sf_logger.warn("Ignoring {}", chrPos);
        return;
      }
    }
    if (m_alleleMap.containsKey(chrPos)) {
      addWarning(chrPos, "Duplicate entry: first valid position wins");
      return;
    }

    if (sampleData.size() > 1) {
      addWarning(chrPos, "Multiple samples found, only using first entry");
    }

    String gt = sampleData.get(0).getProperty("GT");
    if (gt == null) {
      addWarning(chrPos, "Ignoring: no genotype");
      return;
    }
    if (sf_noCallPattern.matcher(gt).matches()) {
      addWarning(chrPos, "Ignoring: no call (" + gt + ")");
      return;
    }

    // normalize alleles to use same syntax as haplotype definition
    List<String> alleles = new ArrayList<>();
    if (position.getAltBases().size() == 0) {
      String gt1 = position.getAllele(0);
      validateAlleles(chrPos, gt1, null);

      String[] g = normalizeAlleles(gt1, null);
      alleles.add(g[0]);

    } else {
      for (String gt2 : position.getAltBases()) {
        String gt1 = position.getAllele(0);
        validateAlleles(chrPos, gt1, gt2);

        String[] g = normalizeAlleles(gt1, gt2);
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

    int[] alleleIdxs = sf_gtDelimiter.splitAsStream(gt)
        .filter(a -> !a.equals("."))
        .mapToInt(Integer::parseInt)
        .toArray();
    String a1 = alleles.get(alleleIdxs[0]);
    String a2 = ".";
    if (alleleIdxs.length > 1) {
      a2 = alleles.get(alleleIdxs[1]);
    } else {
      if (alleleIdxs[0] == 0) {
        addWarning(chrPos, "Ignoring: only a single allele found.  Since it's reference, treating this as a missing position.");
        return;
      }
      // missing allele is guaranteed to produce a no-call if assumeReferenceInDefinitions is true
      addWarning(chrPos, "Only a single allele found.");
    }

    // genotype divided by "|" if phased and "/" if unphased
    // however, we'll treat homozygous as phased
    boolean isPhased = true;
    //noinspection RedundantIfStatement
    if (gt.contains("/") && a2 != null && !a1.equalsIgnoreCase(a2)) {
      isPhased = false;
    }

    List<String> vcfAlleles = new ArrayList<>();
    vcfAlleles.add(position.getRef());
    vcfAlleles.addAll(position.getAltBases());

    SampleAllele sampleAllele = new SampleAllele(position.getChromosome(), position.getPosition(), a1, a2, isPhased,
        vcfAlleles);
    m_alleleMap.put(chrPos, sampleAllele);
  }


  /**
   * Normalize alleles from VCF to match syntax from allele definitions
   */
  private String[] normalizeAlleles(String refAllele, @Nullable String varAllele) {
    Preconditions.checkNotNull(refAllele);

    refAllele = refAllele.toUpperCase();

    if (varAllele != null) {
      varAllele = varAllele.toUpperCase();
      return new String[] { refAllele, varAllele };
    } else {
      return new String[] { refAllele };
    }
  }


  /**
   * Validate GT input per VCF 4.2 specification.
   */
  private static void validateAlleles(String chrPos, String gt1, @Nullable String gt2) {

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
