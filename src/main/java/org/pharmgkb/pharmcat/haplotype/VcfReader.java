package org.pharmgkb.pharmcat.haplotype;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.lang.invoke.MethodHandles;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.SortedMap;
import java.util.TreeMap;
import java.util.regex.Pattern;
import java.util.zip.GZIPInputStream;
import com.google.common.base.Preconditions;
import com.google.common.collect.ImmutableMap;
import com.google.common.collect.SortedSetMultimap;
import com.google.common.collect.TreeMultimap;
import org.checkerframework.checker.nullness.qual.Nullable;
import org.pharmgkb.common.comparator.ChromosomePositionComparator;
import org.pharmgkb.parser.vcf.VcfFormatException;
import org.pharmgkb.parser.vcf.VcfLineParser;
import org.pharmgkb.parser.vcf.VcfParser;
import org.pharmgkb.parser.vcf.model.ContigMetadata;
import org.pharmgkb.parser.vcf.model.FormatMetadata;
import org.pharmgkb.parser.vcf.model.VcfMetadata;
import org.pharmgkb.parser.vcf.model.VcfPosition;
import org.pharmgkb.parser.vcf.model.VcfSample;
import org.pharmgkb.pharmcat.VcfFile;
import org.pharmgkb.pharmcat.definition.DefinitionReader;
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
  public static final String MSG_AD_FORMAT_MISSING = "AD format is not defined.  Assuming AD field is valid.";
  private static final Logger sf_logger = LoggerFactory.getLogger(MethodHandles.lookup().lookupClass());
  private static final Pattern sf_gtDelimiter = Pattern.compile("[|/]");
  private static final Pattern sf_noCallPattern = Pattern.compile("^[.|/]+$");
  private static final Pattern sf_allelePattern = Pattern.compile("^[AaCcGgTt]+$");
  private static final String sf_filterCodeRef = "PCATxREF";
  private static final String sf_filterCodeAlt = "PCATxALT";
  private static final String sf_filterCodeIndel = "PCATxINDEL";
  private final ImmutableMap<String, VariantLocus> m_locationsOfInterest;
  private final ImmutableMap<String, String> m_locationsByGene;
  private final boolean m_findCombinations;
  private final boolean m_useSpecificSample;
  private String m_sampleId;
  private int m_sampleIdx = -1;
  private VcfMetadata m_vcfMetadata;
  private boolean m_adFormatDefined;
  private boolean m_useAdFormat = true;
  private String m_genomeBuild;
  // <chr:position, allele>
  private final SortedMap<String, SampleAllele> m_alleleMap = new TreeMap<>(ChromosomePositionComparator.getComparator());
  // <chr:position, warning>
  private final SortedSetMultimap<String, String> m_warnings = TreeMultimap.create();


  /**
   * Constructor.
   * Reads in VCF file and pull the sample's alleles for positions of interest.
   */
  public VcfReader(DefinitionReader definitionReader, BufferedReader vcfReader, @Nullable String sampleId,
      boolean findCombinations) throws IOException {
    m_locationsOfInterest = definitionReader.getLocationsOfInterest();
    m_locationsByGene = definitionReader.getLocationsByGene();
    m_sampleId = sampleId;
    m_useSpecificSample = m_sampleId != null;
    m_findCombinations = findCombinations;
    read(vcfReader);
  }


  /**
   * Constructor.  Primarily used for testing.
   * This will read in the VCF file and pull the (first) sample's alleles for positions of interest.
   */
  public VcfReader(DefinitionReader definitionReader, Path vcfFile) throws IOException {
    m_locationsOfInterest = definitionReader.getLocationsOfInterest();
    m_locationsByGene = definitionReader.getLocationsByGene();
    m_sampleId = null;
    m_useSpecificSample = false;
    m_findCombinations = false;
    read(vcfFile);
  }

  /**
   * Constructor.  Primarily for testing.
   * This will read in the VCF file and pull the sample's alleles for <em>ALL</em> positions.
   */
  VcfReader(Path vcfFile) throws IOException {
    m_locationsOfInterest = null;
    m_locationsByGene = null;
    m_sampleId = null;
    m_useSpecificSample = false;
    m_findCombinations = false;
    read(vcfFile);
  }


  public @Nullable String getSampleId() {
    return m_sampleId;
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
    Preconditions.checkArgument(VcfFile.isVcfFile(vcfFile), "%s is not a VCF file", vcfFile);

    try (BufferedReader reader = openVcfFile(vcfFile)) {
      read(reader);
    }
  }

  /**
   * Read VCF data via reader.
   */
  private void read(BufferedReader reader) throws IOException {
    // read VCF file
    try (VcfParser vcfParser = new VcfParser.Builder()
        .fromReader(reader)
        .parseWith(this)
        .build()) {
      m_vcfMetadata = vcfParser.parseMetadata();
      FormatMetadata adFormat = m_vcfMetadata.getFormats().get("AD");
      if (adFormat != null) {
        m_adFormatDefined = true;
        String number = adFormat.getPropertiesRaw().get("Number");
        if (!"R".equals(number)) {
          m_useAdFormat = false;
          if (!".".equals(number)) {
            addWarning("VCF", "INFO header for AD has unexpected number (" + number +
                "). Expecting 'R'. Treating number as '.' and ignoring AD field.");
          }
        }
      }
      if (m_useSpecificSample) {
        for (int x = 0; x < m_vcfMetadata.getNumSamples(); x += 1) {
          if (m_sampleId.equals(m_vcfMetadata.getSampleName(x))) {
            m_sampleIdx = x;
            break;
          }
        }
        if (m_sampleIdx == -1) {
          throw new IllegalStateException("Cannot find sample '" + m_sampleId + "'");
        }
      } else {
        m_sampleId = m_vcfMetadata.getSampleName(0);
        m_sampleIdx = 0;
      }
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

    String gt = sampleData.get(m_sampleIdx).getProperty("GT");
    VariantLocus varLoc;
    Set<String> undocumentedVariations = new HashSet<>();
    if (m_locationsOfInterest != null) {
      varLoc = m_locationsOfInterest.get(chrPos);
      if (varLoc == null) {
        sf_logger.warn("Ignoring {}", chrPos);
        return;
      }
      if (!position.getRef().equals(varLoc.getRef())) {
        addWarning(chrPos, "Discarded genotype at this position because REF in VCF (" + position.getRef() +
            ") does not match expected reference (" + varLoc.getRef() + ")");
        return;
      } else if (position.getFilters().contains(sf_filterCodeRef)) {
        addWarning(chrPos, "PharmCAT preprocessor detected REF mismatch (filter " + sf_filterCodeRef +
            ") but this does not match current data.  Was the VCF preprocessed with a different version of PharmCAT?");
      }

      Set<Integer> gtCalls = new HashSet<>();
      if (gt != null && !sf_noCallPattern.matcher(gt).matches()) {
        Arrays.stream(sf_gtDelimiter.split(gt))
            .filter(a -> !a.equals("."))
            .map(Integer::parseInt)
            .forEach(gtCalls::add);
      }

      // strip out structural alts (e.g. <*>)
      List<String> altBases = position.getAltBases().stream()
          .filter(a -> !a.startsWith("<"))
          .toList();
      boolean expectMultiBase = varLoc.getAlts().stream().anyMatch(a -> a.length() > 1);
      boolean hasMultiBase = altBases.stream().anyMatch(a -> a.length() > 1);
      if (expectMultiBase && !hasMultiBase) {
        if (altBases.isEmpty()) {
          addWarning(chrPos, "Genotype at this position has no ALT allele and an indel or repeat is expected. " +
              "PharmCAT cannot validate this position");
        } else {
          addWarning(chrPos, "Genotype at this position has SNPs ( " +
              String.join("/", position.getAltBases()) + ") but PharmCAT expects indel or repeat (" +
              String.join("/", varLoc.getAlts()) + ")");
        }
      } else {
        for (int x = 0; x < altBases.size(); x += 1) {
          if (gtCalls.contains(x + 1)) {
            String alt = altBases.get(x);
            if (!varLoc.getAlts().contains(alt)) {
              undocumentedVariations.add(alt);
            }
          }
        }
        if (!undocumentedVariations.isEmpty()) {
          StringBuilder msgBuilder = new StringBuilder()
              .append("The genetic variation at this position does not match what is in the allele definition " +
                  "(expected ")
              .append(String.join("/", varLoc.getAlts()))
              .append(", found ")
              .append(String.join("/", undocumentedVariations))
              .append(" in VCF)");
          if (treatAsReference(chrPos)) {
            msgBuilder.append(".  Undocumented variations will be replaced with reference.");
          }
          addWarning(chrPos, msgBuilder.toString());
        } else if (position.getFilters().contains(sf_filterCodeAlt) && sampleData.size() == 1) {
          addWarning(chrPos, "PharmCAT preprocessor detected ALT mismatch (filter " + sf_filterCodeAlt +
              ") but this does not match current data (expected " + String.join("/", varLoc.getAlts()) +
              " and got " + String.join("/", position.getAltBases()) +
              ").  Was the VCF preprocessed with a different version of PharmCAT?");
        }
      }

    } else {
      // for some reason we don't have locations of interest, so pass on warnings from preprocessor
      if (position.getFilters().contains(sf_filterCodeRef)) {
        addWarning(chrPos, "Discarded genotype at this position because REF in VCF (" + position.getRef() +
            ") does not match expected reference");
        return;
      }
      if (position.getFilters().contains(sf_filterCodeAlt)) {
        addWarning(chrPos, "The genetic variation at this position does not match what is in the allele definition");
      }
      if (position.getFilters().contains(sf_filterCodeIndel)) {
        addWarning(chrPos, "Genotype at this position uses unexpected format for INDEL");
      }
    }
    if (m_alleleMap.containsKey(chrPos)) {
      addWarning(chrPos, "Duplicate entry: first valid position wins");
      return;
    }

    if (sampleData.size() > 1 && !m_useSpecificSample) {
      addWarning(chrPos, "Multiple samples found, only using first entry.");
    }

    if (gt == null) {
      addWarning(chrPos, "Ignoring: no genotype");
      return;
    }
    if (sf_noCallPattern.matcher(gt).matches()) {
      addWarning(chrPos, "Ignoring: no call (" + gt + ")");
      return;
    }
    // already filtered out no-calls, guaranteed to have at least 1 GT here
    String[] gtArray = sf_gtDelimiter.split(gt);
    int[] alleleIndices = Arrays.stream(gtArray)
        .filter(a -> !a.equals("."))
        .mapToInt(Integer::parseInt)
        .toArray();
    boolean isHaploid = gtArray.length == 1;

    // normalize alleles to use the same syntax as haplotype definition
    List<String> alleles = new ArrayList<>();
    if (position.getAltBases().isEmpty()) {
      String gt1 = position.getAllele(0);
      if (!validateAlleles(chrPos, gt1, null, false)) {
        return;
      }

      String[] g = normalizeAlleles(gt1, null);
      alleles.add(g[0]);

    } else {
      String gt1 = position.getAllele(0);
      List<Integer> selected = Arrays.stream(alleleIndices).boxed().toList();
      for (int x = 1; x <= position.getAltBases().size(); x += 1) {
        String gt2 = position.getAllele(x);
        boolean isSelected = selected.contains(x);
        if (!validateAlleles(chrPos, gt1, gt2, isSelected)) {
          return;
        }

        String[] g = normalizeAlleles(gt1, gt2);
        String g1 = g[0];
        String g2 = g[1];
        if (alleles.isEmpty()) {
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


    if (m_useAdFormat) {
      // reference: https://gatk.broadinstitute.org/hc/en-us/articles/360035531692-VCF-Variant-Call-Format
      String allelicDepth = sampleData.get(0).getProperty("AD");
      if (allelicDepth != null) {
        if (!m_adFormatDefined) {
          addWarning("VCF", MSG_AD_FORMAT_MISSING);
        }
        if (!".".equals(allelicDepth)) {
          // try to catch reference overlap style VCF where
          // VCF always specifies heterozygous GT and uses AD field to determine actual alleles
          // see https://github.com/PharmGKB/PharmCAT/issues/90 for example
          try {
            List<Integer> depths = Arrays.stream(allelicDepth.split(","))
                .filter(a -> !a.equals("."))
                .map(Integer::parseInt)
                .toList();
            Map<Integer, Integer> genotype = new HashMap<>();
            Arrays.stream(alleleIndices)
                .forEach(g -> genotype.merge(g, 1, Integer::sum));
            // GT is het, but AD is not (only one side has any reads)
            if (genotype.size() != 1 && depths.stream().filter(d -> d > 0).count() == 1) {
              addWarning(chrPos, "Discarding genotype at this position because GT field indicates heterozygous (" +
                  gt + ") but AD field indicates homozygous (" + allelicDepth + ")");
              return;
            }
          } catch (NumberFormatException ex) {
            addWarning(chrPos, "Invalid allelic depth (AD) field: " + allelicDepth);
          }
        }
      }
    }

    for (int alleleIdx : alleleIndices) {
      if (alleleIdx >= alleles.size()) {
        throw new VcfFormatException("Invalid GT allele value (" + alleleIdx + ") for " + chrPos +
            " (only " + (alleles.size() - 1) + " ALT allele" + (alleles.size() > 1 ? "s" : "") + " specified)");
      }
    }

    String a1 = alleles.get(alleleIndices[0]);
    String a2 = isHaploid ? null : ".";
    if (alleleIndices.length > 1) {
      a2 = alleles.get(alleleIndices[1]);
      if (alleleIndices.length > 2) {
        addWarning(chrPos, alleleIndices.length + " genotypes found.  Only using first 2 genotypes.");
      }
    } else if (!position.getChromosome().equals("chrM") &&
        !position.getChromosome().equals("chrY") &&
        !position.getChromosome().equals("chrX")) {
      if (alleleIndices[0] == 0) {
        // treating "./0" like any other missing position
        addWarning(chrPos, "Ignoring: only a single genotype found.  " +
            "Since it's reference, treating this as a missing position.");
        return;
      }
      // going to accept "./1" because there's a variation, and we want to pass this on to the reporter
      addWarning(chrPos, "Only a single genotype found.");
    }

    // genotype divided by "|" if phased and "/" if unphased
    boolean isPhased = true;
    boolean isEffectivelyPhased = true;
    if (gt.contains("/")) {
      isPhased = false;
      // we'll also treat homozygous as phased
       if (a2 != null && !a1.equalsIgnoreCase(a2)) {
         isEffectivelyPhased = false;
       }
    }

    List<String> vcfAlleles = new ArrayList<>();
    vcfAlleles.add(position.getRef());
    vcfAlleles.addAll(position.getAltBases());

    boolean treatUndocumentedAsReference = m_locationsByGene != null && treatAsReference(chrPos);
    SampleAllele sampleAllele = new SampleAllele(position.getChromosome(), position.getPosition(), a1, a2, isPhased,
        isEffectivelyPhased, vcfAlleles, undocumentedVariations, treatUndocumentedAsReference);
    m_alleleMap.put(chrPos, sampleAllele);
  }


  private boolean treatAsReference(String chrPos) {
    String gene = m_locationsByGene.get(chrPos);
    if (NamedAlleleMatcher.TREAT_UNDOCUMENTED_VARIATIONS_AS_REFERENCE.contains(gene)) {
      if ("DPYD".equals(gene)) {
        // DPYD ignores find-combinations mode
        return true;
      }
      return !m_findCombinations;
    }
    return false;
  }


  /**
   * Normalize alleles from VCF to match syntax from allele definitions.
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
  private boolean validateAlleles(String chrPos, String gt1, @Nullable String gt2, boolean isG2Selected) {

    boolean isValid = true;
    if (gt1.startsWith("<")) {
      addWarning(chrPos, "Discarded genotype at this position because REF uses structural variation '" + gt1 + "'");
      isValid = false;
    } else if (gt1.toUpperCase().contains("N")) {
      addWarning(chrPos, "Discarded genotype at this position because REF uses ambiguous allele in '" + gt1 + "'");
      isValid = false;
    } else if (!sf_allelePattern.matcher(gt1).matches()) {
      addWarning(chrPos, "Discarded genotype at this position because REF uses unknown base in '" + gt1 + "'");
      isValid = false;
    }

    if (gt2 != null) {
      String prefix = isG2Selected ? "Discarded genotype at this position because " : "";
      if (gt2.startsWith("<")) {
        if (isG2Selected) {
          isValid = false;
        }
        addWarning(chrPos, prefix + "ALT uses structural variation '" + gt2 + "'");
      } else if (gt2.toUpperCase().contains("N")) {
        if (isG2Selected) {
          isValid = false;
        }
        addWarning(chrPos, prefix + "ALT uses ambiguous allele in '" + gt2 + "'");

      } else if (gt2.contains("*")) {
        if (isG2Selected) {
          isValid = false;
        }
        addWarning(chrPos, prefix + "ALT uses missing allele in '" + gt2 + "'");
      } else if (!sf_allelePattern.matcher(gt2).matches()) {
        if (isG2Selected) {
          isValid = false;
        }
        addWarning(chrPos, prefix + "ALT uses unknown base in '" + gt2 + "'");
      }
    }

    return isValid;
  }


  static BufferedReader openVcfFile(Path vcfFile) throws  IOException {
    String filename = vcfFile.toString();
    if (filename.endsWith(".vcf.bgz") || filename.endsWith(".vcf.gz")) {
      return new BufferedReader(new InputStreamReader(new GZIPInputStream(Files.newInputStream(vcfFile))));
    }
    return Files.newBufferedReader(vcfFile);
  }
}
