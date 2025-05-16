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
import com.google.common.collect.ImmutableSet;
import com.google.common.collect.SortedSetMultimap;
import com.google.common.collect.TreeMultimap;
import org.jspecify.annotations.Nullable;
import org.pharmgkb.common.comparator.ChromosomePositionComparator;
import org.pharmgkb.parser.vcf.VcfFormatException;
import org.pharmgkb.parser.vcf.VcfLineParser;
import org.pharmgkb.parser.vcf.VcfParser;
import org.pharmgkb.parser.vcf.model.ContigMetadata;
import org.pharmgkb.parser.vcf.model.FormatMetadata;
import org.pharmgkb.parser.vcf.model.VcfMetadata;
import org.pharmgkb.parser.vcf.model.VcfPosition;
import org.pharmgkb.parser.vcf.model.VcfSample;
import org.pharmgkb.pharmcat.ParseException;
import org.pharmgkb.pharmcat.VcfFile;
import org.pharmgkb.pharmcat.definition.DefinitionReader;
import org.pharmgkb.pharmcat.definition.model.VariantLocus;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import static org.pharmgkb.pharmcat.Constants.isLowestFunctionGene;


/**
 * This class reads VCF files and pulls the sample's alleles for positions of interest (i.e. is necessary to make a
 * haplotype call).
 *
 * @author Mark Woon
 */
public class VcfReader implements VcfLineParser {
  public static final String MSG_AD_FORMAT_MISSING = "AD format is not defined.  Assuming AD field is valid.";
  public static final Pattern GT_DELIMITER = Pattern.compile("[|/]");
  private static final Logger sf_logger = LoggerFactory.getLogger(MethodHandles.lookup().lookupClass());
  private static final Set<String> sf_haploidChromosomes = ImmutableSet.of("chrY", "chrM");
  private static final Pattern sf_allelePattern = Pattern.compile("^[AaCcGgTt]+$");
  private static final String sf_filterCodeRef = "PCATxREF";
  private static final String sf_filterCodeAlt = "PCATxALT";
  private static final String sf_filterCodeIndel = "PCATxINDEL";
  private final @Nullable ImmutableMap<String, VariantLocus> m_locationsOfInterest;
  private final @Nullable ImmutableMap<String, String> m_locationsByGene;
  private final boolean m_findCombinations;
  private @Nullable String m_sampleId;
  private int m_sampleIdx = -1;
  private @Nullable VcfMetadata m_vcfMetadata;
  private boolean m_adFormatDefined;
  private boolean m_useAdFormat = true;
  private @Nullable String m_genomeBuild;
  // <chr:position, allele>
  private final SortedMap<String, SampleAllele> m_alleleMap = new TreeMap<>(ChromosomePositionComparator.getComparator());
  // <chr:position, warning>
  private final SortedSetMultimap<String, String> m_warnings = TreeMultimap.create();
  private final Set<String> m_discardedPositions = new HashSet<>();


  /**
   * Constructor.
   * Reads in a VCF file and pulls the sample's alleles at positions of interest.
   *
   * @throws ParseException if there are no samples in the VCF file
   */
  public VcfReader(DefinitionReader definitionReader, BufferedReader vcfReader, @Nullable String sampleId,
      boolean findCombinations) throws IOException, ParseException {
    m_locationsOfInterest = definitionReader.getLocationsOfInterest();
    m_locationsByGene = definitionReader.getLocationsByGene();
    m_sampleId = sampleId;
    m_findCombinations = findCombinations;
    read(vcfReader);
  }


  /**
   * Constructor.  Primarily used for testing.
   * This will read in the VCF file and pull the (first) sample's alleles at positions of interest.
   *
   * @throws ParseException if there are no samples in the VCF file
   */
  public VcfReader(DefinitionReader definitionReader, Path vcfFile) throws IOException {
    m_locationsOfInterest = definitionReader.getLocationsOfInterest();
    m_locationsByGene = definitionReader.getLocationsByGene();
    m_sampleId = null;
    m_findCombinations = false;
    read(vcfFile);
  }

  /**
   * Constructor.  Primarily for testing.
   * This will read in the VCF file and pull the (first) sample's alleles at <em>ALL</em> positions.
   *
   * @throws ParseException if there are no samples in the VCF file
   */
  VcfReader(Path vcfFile) throws IOException {
    this(vcfFile, null);
  }

  /**
   * Constructor.  Primarily for testing.
   * This will read in the VCF file and pull the sample's alleles at <em>ALL</em> positions.
   *
   * @throws ParseException if there are no samples in the VCF file
   */
  VcfReader(Path vcfFile, @Nullable String sampleId) throws IOException {
    m_locationsOfInterest = null;
    m_locationsByGene = null;
    m_sampleId = sampleId;
    m_findCombinations = false;
    read(vcfFile);
  }


  /**
   * Gets the Sample ID of the data to read.
   */
  public @Nullable String getSampleId() {
    // this should never be null after read() is called
    return m_sampleId;
  }


  public @Nullable VcfMetadata getVcfMetadata() {
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
   * Reads the VCF file.
   *
   * @throws ParseException if there are no samples in the VCF file
   */
  private void read(Path vcfFile) throws IOException, ParseException {
    Preconditions.checkNotNull(vcfFile);
    Preconditions.checkArgument(VcfFile.isVcfFile(vcfFile), "%s is not a VCF file", vcfFile);

    try (BufferedReader reader = openVcfFile(vcfFile)) {
      read(reader);
    }
  }

  /**
   * Read VCF data via reader.
   *
   * @throws ParseException if there are no samples in the VCF file
   */
  private void read(BufferedReader reader) throws IOException, ParseException {
    // read VCF
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
      if (m_sampleId != null) {
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
    addWarning(chrPos, msg, msg);
  }

  private void addWarning(String chrPos, String reportMsg, String clMsg) {
    // saved to report
    m_warnings.put(chrPos, reportMsg);
    // prints to the command line
    sf_logger.warn(clMsg);
  }


  /**
   * Parses a single line from VCF.
   *
   * @throws ParseException if there are no samples in the VCF file
   */
  @Override
  public void parseLine(VcfMetadata metadata, VcfPosition position, List<VcfSample> sampleData) throws ParseException {

    if (sampleData.isEmpty()) {
      throw new ParseException("VCF does not contain sample data");
    }

    final String chrPos = position.getChromosome() + ":" + position.getPosition();
    if (m_alleleMap.containsKey(chrPos)) {
      // The preprocessor will split out PGx variants and non-PGx variants into 2 separate lines and add a filter code
      // on the second entry. Don't bother printing warning in this case.
      if (!position.getFilters().contains(sf_filterCodeRef) &&
          !position.getFilters().contains(sf_filterCodeAlt) &&
          !position.getFilters().contains(sf_filterCodeIndel)) {
        addWarning(chrPos, "Duplicate entry found in VCF; first valid entry trumps others.",
            "Duplicate entry: first valid position wins");
      }
      return;
    }

    // only care about locations of interest
    VariantLocus varLoc = null;
    if (m_locationsOfInterest != null) {
      varLoc = m_locationsOfInterest.get(chrPos);
      if (varLoc == null) {
        sf_logger.warn("Ignoring {}", chrPos);
        return;
      }
      if (!position.getRef().equals(varLoc.getRef())) {
        addWarning(chrPos, "Discarded genotype at this position because REF in VCF (" + position.getRef() +
            ") does not match expected reference (" + varLoc.getRef() + ")");
        m_discardedPositions.add(chrPos);
        return;
      } else if (position.getFilters().contains(sf_filterCodeRef)) {
        addWarning(chrPos, "PharmCAT preprocessor detected REF mismatch (filter " + sf_filterCodeRef +
            ") but this does not match current data.  Was the VCF preprocessed with a different version of PharmCAT?");
      }
    } else {
      // for some reason we don't have locations of interest, so pass on warnings from preprocessor
      if (position.getFilters().contains(sf_filterCodeRef)) {
        addWarning(chrPos, "Discarded genotype at this position because REF in VCF (" + position.getRef() +
            ") does not match expected reference");
        m_discardedPositions.add(chrPos);
        return;
      }
      if (position.getFilters().contains(sf_filterCodeAlt)) {
        addWarning(chrPos, "The genetic variation at this position does not match what is in the allele definition");
      }
      if (position.getFilters().contains(sf_filterCodeIndel)) {
        addWarning(chrPos, "Genotype at this position uses unexpected format for INDEL");
      }
    }


    String gt = sampleData.get(m_sampleIdx).getProperty("GT");
    if (gt == null) {
      addWarning(chrPos, "Ignoring: no genotype");
      m_discardedPositions.add(chrPos);
      return;
    }
    String[] gtArray = GT_DELIMITER.split(gt);
    List<Integer> gtNonMissing = Arrays.stream(gtArray)
        .filter(g -> !g.equals("."))
        .map(Integer::parseInt)
        .toList();
    if (gtNonMissing.isEmpty()) {
      addWarning(chrPos, "Ignoring: no call (" + gt + ")");
      m_discardedPositions.add(chrPos);
      return;
    }

    if (sf_haploidChromosomes.contains(position.getChromosome())) {
      // expect a single allele
      if (gtNonMissing.size() > 1) {
        addWarning(chrPos, gtNonMissing.size() + " genotypes found (GT=" + gt + getGeneForWarning(chrPos) +
            ") for haploid chromosome. Will only use first non-missing genotype.");
      }
    } else {
      // diploid chromosome
      if (gtNonMissing.size() > 2) {
        addWarning(chrPos, gtNonMissing.size() + " genotypes found (GT=" + gt + getGeneForWarning(chrPos) +
            "). Will only use first two genotypes.");
      } else if (gtNonMissing.size() == 1) {
        if (!position.getChromosome().equals("chrX")) {
          if (gtNonMissing.get(0) == 0) {
            // treating "./0" or "0/." like any other missing position
            addWarning(chrPos, "Ignoring: only a single genotype found (GT=" + gt +
                ").  Since it's reference, treating this as a missing position.");
            m_discardedPositions.add(chrPos);
            return;
          }
          addWarning(chrPos, gtNonMissing.size() + " genotype found (GT=" + gt + getGeneForWarning(chrPos) +
              "), expecting 2.");
        }
      }
    }

    Set<String> undocumentedVariations = new HashSet<>();
    Set<Integer> gtUnique = new HashSet<>(gtNonMissing);
    if (varLoc != null) {
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
          addWarning(chrPos, "Genotype at this position has SNP" +
              (position.getAltBases().size() > 1 ? "s" : "") + " (" +
              String.join("/", position.getAltBases()) + ") but PharmCAT expects indel or repeat (" +
              String.join("/", varLoc.getAlts()) + ")");
        }
      } else {
        for (int x = 0; x < altBases.size(); x += 1) {
          if (gtUnique.contains(x + 1)) {
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
          if (treatUndocumentedAsReference(chrPos)) {
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
    }

    if (sampleData.size() > 1 && m_sampleId == null) {
      // only warn once
      if (m_alleleMap.isEmpty()) {
        addWarning(chrPos, "Multiple samples found, only using first entry.");
      }
    }

    // normalize alleles to use the same syntax as haplotype definition
    List<String> alleles = new ArrayList<>();
    if (position.getAltBases().isEmpty()) {
      String gt1 = position.getAllele(0);
      if (!validateAlleles(chrPos, gt1, null, false)) {
        m_discardedPositions.add(chrPos);
        return;
      }

      String[] g = normalizeAlleles(gt1, null);
      alleles.add(g[0]);

    } else {
      String gt1 = position.getAllele(0);
      for (int x = 1; x <= position.getAltBases().size(); x += 1) {
        String gt2 = position.getAllele(x);
        boolean isSelected = gtNonMissing.contains(x);
        if (!validateAlleles(chrPos, gt1, gt2, isSelected)) {
          m_discardedPositions.add(chrPos);
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
      String allelicDepth = sampleData.get(m_sampleIdx).getProperty("AD");
      if (allelicDepth != null) {
        if (!m_adFormatDefined) {
          addWarning("VCF", MSG_AD_FORMAT_MISSING);
        }
        if (!".".equals(allelicDepth)) {
          // try to catch reference overlap style VCF where
          // VCF always specifies heterozygous GT and uses AD field to determine actual alleles
          // see https://github.com/PharmGKB/PharmCAT/issues/90 for details
          try {
            List<Integer> depths = Arrays.stream(allelicDepth.split(","))
                .filter(a -> !a.equals("."))
                .map(Integer::parseInt)
                .toList();
            Map<Integer, Integer> genotype = new HashMap<>();
            gtNonMissing.forEach(g -> genotype.merge(g, 1, Integer::sum));
            // GT is het, but AD is not (only one side has any reads)
            if (genotype.size() != 1 && depths.stream().filter(d -> d > 0).count() == 1) {
              addWarning(chrPos, "Discarding genotype at this position because GT field indicates heterozygous (" +
                  gt + ") but AD field indicates homozygous (" + allelicDepth + ")");
              m_discardedPositions.add(chrPos);
              return;
            }
          } catch (NumberFormatException ex) {
            addWarning(chrPos, "Invalid allelic depth (AD) field: " + allelicDepth);
          }
        }
      }
    }

    for (int alleleIdx : gtNonMissing) {
      if (alleleIdx >= alleles.size()) {
        throw new VcfFormatException("Invalid GT allele value (" + alleleIdx + ") for " + chrPos +
            " (only " + (alleles.size() - 1) + " ALT allele" + (alleles.size() > 1 ? "s" : "") + " specified)");
      }
    }

    String a1 = null;
    if (!gtArray[0].equals(".")) {
      a1 = alleles.get(Integer.parseInt(gtArray[0]));
    }
    String a2 = null;
    if (!sf_haploidChromosomes.contains(position.getChromosome()) || a1 == null) {
      if (gtArray.length > 1 && !gtArray[1].equals(".")) {
        a2 = alleles.get(Integer.parseInt(gtArray[1]));
      }
    }

    // genotype divided by "|" if phased and "/" if unphased
    boolean isPhased = true;
    boolean isEffectivelyPhased = true;
    Integer phaseSet = null;
    if (gt.contains("/")) {
      isPhased = false;
      if (gtUnique.size() > 1) {
        // we'll also treat homozygous as phased
        isEffectivelyPhased = false;
       }
    } else {
      String p = sampleData.get(m_sampleIdx).getProperty("PS");
      if (p != null && !p.equals(".")) {
        phaseSet = Integer.valueOf(p);
        if (gtUnique.size() > 1) {
          isEffectivelyPhased = false;
        }
      }
    }

    List<String> vcfAlleles = new ArrayList<>();
    vcfAlleles.add(position.getRef());
    vcfAlleles.addAll(position.getAltBases());

    boolean treatUndocumentedAsReference = m_locationsByGene != null && treatUndocumentedAsReference(chrPos);
    SampleAllele sampleAllele = new SampleAllele(position.getChromosome(), position.getPosition(), a1, a2, isPhased,
        isEffectivelyPhased, phaseSet, vcfAlleles, gt, undocumentedVariations, treatUndocumentedAsReference);
    m_alleleMap.put(chrPos, sampleAllele);
    if (m_discardedPositions.contains(chrPos)) {
      addWarning(chrPos, "Duplicate entry found in VCF; this entry trumps previous invalid entry.",
          "Duplicate entry: first valid position wins");
    }
  }

  /**
   * Gets gene {@code chrPos} is in surrounded by parentheses for use in warnings.
   * <p>
   * Example: {@code " (TPMT)"}
   */
  private String getGeneForWarning(String chrPos) {
    if (m_locationsByGene == null) {
      return "";
    }
    return ", for " + m_locationsByGene.get(chrPos);
  }

  private boolean treatUndocumentedAsReference(String chrPos) {
    if (m_locationsByGene == null) {
      return false;
    }
    String gene = m_locationsByGene.get(chrPos);
    if (gene != null && NamedAlleleMatcher.TREAT_UNDOCUMENTED_VARIATIONS_AS_REFERENCE.contains(gene)) {
      if (isLowestFunctionGene(gene)) {
        // lowest-function genes ignore find-combinations mode
        return true;
      }
      // if finding combinations, don't treat as reference so that the exact change is reflected
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
