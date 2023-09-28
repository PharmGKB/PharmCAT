package org.pharmgkb.pharmcat.definition.model;

import java.io.IOException;
import java.util.*;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import java.util.stream.Collectors;
import com.google.gson.annotations.Expose;
import com.google.gson.annotations.SerializedName;
import org.pharmgkb.pharmcat.ParseException;
import org.pharmgkb.pharmcat.haplotype.Iupac;
import org.pharmgkb.pharmcat.reporter.model.DataSource;
import org.pharmgkb.pharmcat.util.DataManager;
import org.pharmgkb.pharmcat.util.VcfHelper;


/**
 * This class represents one complete allele translation set for a gene.
 *
 * @author Ryan Whaley
 */
public class DefinitionFile {
  public static final String FORMAT_VERSION = "2";
  // metadata
  @Expose
  @SerializedName("formatVersion")
  private String m_formatVersion = FORMAT_VERSION;
  @Expose
  @SerializedName(value = "source")
  private DataSource m_source;
  @Expose
  @SerializedName(value = "version", alternate = {"cpicVersion"})
  private String m_version;
  @Expose
  @SerializedName("modificationDate")
  private Date m_modificationDate;
  @Expose
  @SerializedName("gene")
  private String m_geneSymbol;
  @Expose
  @SerializedName("orientation")
  private String m_orientation;
  @Expose
  @SerializedName("chromosome")
  private String m_chromosome;
  @Expose
  @SerializedName("genomeBuild")
  private String m_genomeBuild;
  @Expose
  @SerializedName("refSeqChromosomeId")
  private String m_refSeqChromosome;
  // definitions
  @Expose
  @SerializedName("variants")
  private VariantLocus[] m_variants;
  @Expose
  @SerializedName("namedAlleles")
  private SortedSet<NamedAllele> m_namedAlleles;


  /**
   * The format version of the definition file.
   */
  public String getFormatVersion() {
    return m_formatVersion;
  }

  /**
   * The source that this definition was pulled from.
   */
  public DataSource getSource() {
    return m_source;
  }

  public void setSource(DataSource source) {
    m_source = source;
  }


  /**
   * The version of this definition.
   */
  public String getVersion() {
    return m_version;
  }


  /**
   * The date this file should be considered to be last modified (should be manually set by curators).
   */
  public Date getModificationDate() {
    return m_modificationDate;
  }


  /**
   * The symbol of the gene these alleles are on.
   */
  public String getGeneSymbol() {
    return m_geneSymbol;
  }

  /**
   * The orientation of the gene relative to the chromosome
   */
  public String getOrientation() {
    return m_orientation;
  }

  /**
   * The name of the chromosome this translation is on.
   */
  public String getChromosome() {
    return m_chromosome;
  }

  /**
   * The human genome assembly (build) the positions in this translation are from (e.g. b38 or b37)
   */
  public String getGenomeBuild() {
    return m_genomeBuild;
  }


  /**
   * The RefSeq identifier for the chromosome this translation is on (should agree with build).
   */
  public String getRefSeqChromosome() {
    return m_refSeqChromosome;
  }


  /**
   * The {@link VariantLocus} objects used to define {@link NamedAllele}s in this translation
   */
  public VariantLocus[] getVariants() {
    return m_variants;
  }


  /**
   * All the named alleles defined in this translation
   */
  public SortedSet<NamedAllele> getNamedAlleles() {
    return m_namedAlleles;
  }


  @Override
  public String toString() {
    return "Allele definition for " + m_geneSymbol;
  }


  @Override
  public boolean equals(Object o) {
    if (this == o) {
      return true;
    }
    if (!(o instanceof DefinitionFile that)) {
      return false;
    }
    return Objects.equals(m_formatVersion, that.getFormatVersion()) &&
        Objects.equals(m_modificationDate, that.getModificationDate()) &&
        Objects.equals(m_geneSymbol, that.getGeneSymbol()) &&
        Objects.equals(m_orientation, that.getOrientation()) &&
        Objects.equals(m_chromosome, that.getChromosome()) &&
        Objects.equals(m_genomeBuild, that.getGenomeBuild()) &&
        Objects.equals(m_refSeqChromosome, that.getRefSeqChromosome()) &&
        Arrays.equals(m_variants, that.getVariants()) &&
        Objects.equals(m_namedAlleles, that.getNamedAlleles());
  }

  @Override
  public int hashCode() {
    return Objects.hash(m_formatVersion, m_modificationDate, m_geneSymbol, m_orientation,
        m_chromosome, m_genomeBuild, m_refSeqChromosome, Arrays.hashCode(m_variants), m_namedAlleles);
  }


  /**
   * Remove ignored allele specified in {@link DefinitionExemption}.
   * Should only be called during initial generation of this {@link DefinitionFile} by {@link DataManager}.
   */
  public void removeIgnoredNamedAlleles(DefinitionExemption exemption) {
    System.out.println("  Removing " + exemption.getIgnoredAlleles());
    m_namedAlleles = m_namedAlleles.stream()
        .filter((na) -> !exemption.shouldIgnoreAllele(na.getName()))
        .collect(Collectors.toCollection(TreeSet::new));
  }


  /**
   * Remove ignored positions specified in {@link DefinitionExemption}.
   * Should only be called during initial generation of this {@link DefinitionFile} by {@link DataManager}.
   */
  public void removeIgnoredPositions(DefinitionExemption exemption) {
    // find ignored positions
    Set<Integer> ignoredPositions = new HashSet<>();
    List<VariantLocus> newVariants = new ArrayList<>();
    for (int x = 0; x < m_variants.length; x += 1) {
      if (exemption.shouldIgnorePosition(m_variants[x])) {
        System.out.println("  Removing position " + x + " (" + m_variants[x] + ")");
        ignoredPositions.add(x);
      } else {
        newVariants.add(m_variants[x]);
      }
    }
    if (exemption.getIgnoredPositions().size() != ignoredPositions.size()) {
      throw new IllegalStateException("Should have " + exemption.getIgnoredPositions().size() + " ignored positions, " +
          "but only found " + ignoredPositions.size());
    }
    // update variants
    m_variants = newVariants.toArray(new VariantLocus[0]);

    SortedSet<NamedAllele> updatedNamedAlleles = new TreeSet<>();
    for (NamedAllele namedAllele : m_namedAlleles) {
      String[] cpicAlleles = new String[namedAllele.getCpicAlleles().length - ignoredPositions.size()];
      if (m_variants.length != cpicAlleles.length) {
        throw new IllegalStateException("Number of variants (" + m_variants.length + ") and number of CPIC alleles (" +
            cpicAlleles.length + ") don't match up for " + namedAllele.getName());
      }
      for (int x = 0, y = 0; x < namedAllele.getCpicAlleles().length; x += 1) {
        if (ignoredPositions.contains(x)) {
          continue;
        }
        cpicAlleles[y] = namedAllele.getCpicAlleles()[x];
        y += 1;
      }
      // if there's nothing left that differs from reference allele then don't include this named allele in output
      if (!Arrays.stream(cpicAlleles).allMatch(Objects::isNull)) {
        updatedNamedAlleles.add(new NamedAllele(namedAllele.getId(), namedAllele.getName(), null, cpicAlleles,
            namedAllele.isReference()));
      }
    }
    m_namedAlleles = updatedNamedAlleles;
  }


  /**
   * Translate variants from CPIC to VCF (i.e. {@code cpicAlleles} to {@code alleles}).
   */
  public void doVcfTranslation(VcfHelper vcfHelper) throws IOException {

    NamedAllele referenceNamedAllele = null;
    for (NamedAllele na : m_namedAlleles) {
      na.initializeForImport(m_variants);
      if (na.isReference()) {
        referenceNamedAllele = na;
      }
    }
    if (referenceNamedAllele == null) {
      throw new IllegalStateException(m_geneSymbol + " does not have reference named allele");
    }

    for (VariantLocus vl : m_variants) {
      translateVariantLocus(referenceNamedAllele, vl, vcfHelper);
    }

    VariantLocus[] sortedVariants = new VariantLocus[m_variants.length];
    System.arraycopy(m_variants, 0, sortedVariants, 0, m_variants.length);
    Arrays.sort(sortedVariants);
    boolean mustResort = !Arrays.equals(sortedVariants, m_variants);

    SortedSet<NamedAllele> updatedNamedAlleles = new TreeSet<>();
    for (NamedAllele na : m_namedAlleles) {
      String[] fixedAlleles = new String[na.getCpicAlleles().length];
      for (int x = 0; x < fixedAlleles.length; x += 1) {
        String cpicAllele = na.getCpicAlleles()[x];
        if (cpicAllele != null) {
          VariantLocus vl = m_variants[x];
          fixedAlleles[x] = vl.getCpicToVcfAlleleMap().get(cpicAllele);
          if (fixedAlleles[x] == null) {
            if (Iupac.lookup(cpicAllele).isAmbiguity()) {
              fixedAlleles[x] = cpicAllele;
            } else {
              throw new IllegalStateException("Don't know how to translate CPIC allele '" + cpicAllele + "' @ " +
                  " position " + x + " (" + vl + ") for " + m_geneSymbol + " " + na.getName() +
                  "; expecting " + vl.getRef() + " / " + String.join(" / ", vl.getAlts()));
            }
          }
        }
      }
      NamedAllele updated;
      if (mustResort) {
        updated = reorderHaplotypeAlleles(na, m_variants, sortedVariants, fixedAlleles);
      } else {
        updated = new NamedAllele(na.getId(), na.getName(), fixedAlleles, na.getCpicAlleles(), na.isReference());
      }
      updatedNamedAlleles.add(updated);
    }
    m_namedAlleles = Collections.unmodifiableSortedSet(updatedNamedAlleles);
    m_variants = sortedVariants;
  }

  private static final Pattern sf_hgvsRepeatPattern = Pattern.compile("g\\.[\\d_]+([ACGT]+\\[\\d+])$");
  private static final Pattern sf_hgvsInsPattern = Pattern.compile("g\\.[\\d_]+(del[ACGT]*)?(ins[ACGT]+)$");
  private static final Pattern sf_hgvsDelPattern = Pattern.compile("g\\.[\\d_]+del[ACGT]*$");

  private void translateVariantLocus(NamedAllele referenceNamedAllele, VariantLocus vl, VcfHelper vcfHelper)
      throws IOException {
    String errorLocation = m_geneSymbol + " @ " + vl.getCpicPosition();

    String refAllele = Objects.requireNonNull(referenceNamedAllele.getCpicAllele(vl));
    if (Iupac.isWobble(refAllele)) {
      String hgvs = VariantLocus.HGVS_NAME_SPLITTER.splitToList(vl.getChromosomeHgvsName()).get(0);
      VcfHelper.VcfData vcf = vcfHelper.hgvsToVcf(m_refSeqChromosome + ":" + hgvs);
      System.out.println(m_geneSymbol + " reference (" + referenceNamedAllele.getName() +
          ") @ " + vl.getCpicPosition() + " is ambiguous (" + refAllele + "):  using " + vcf.ref + " for VCF");
      refAllele = vcf.ref;
    }
    SortedSet<String> altAlleles = new TreeSet<>();
    boolean isSnp = true;
    List<String> repeats = new ArrayList<>();
    List<String> nonRepeats = new ArrayList<>();
    for (String allele : vl.getCpicAlleles()) {
      if (allele.length() > 1) {
        isSnp = false;
      }
      if (allele.contains("(") || allele.contains(")")) {
        if (!allele.contains("(") || !allele.contains(")")) {
          throw new IllegalStateException(errorLocation + ": allele has mismatched parentheses - " + allele);
        }
        repeats.add(allele);
      } else {
        nonRepeats.add(allele);
      }
      if (allele.contains("[") || allele.contains("]")) {
        throw new IllegalStateException(errorLocation + ": allele uses square brackets - " + allele);
      }
      if (!allele.equals(refAllele) && !Iupac.isWobble(allele)) {
        altAlleles.add(allele);
      }
    }
    if (repeats.size() > 0 && repeats.size() != vl.getCpicAlleles().size()) {
      boolean haveSingle = false;
      if (nonRepeats.size() == 1) {
        String repeatedSequence = repeats.get(0);
        repeatedSequence = repeatedSequence.substring(0, repeatedSequence.indexOf("("));
        haveSingle = nonRepeats.get(0).equals(repeatedSequence);
      }
      if (!haveSingle) {
        throw new IllegalStateException(errorLocation + ": has " + repeats.size() + " repeat alleles but " +
            vl.getCpicAlleles().size() + " total alleles (" + vl.getCpicAlleles() + ")");
      }
    }

    List<String> hgvsNames = VariantLocus.HGVS_NAME_SPLITTER.splitToList(vl.getChromosomeHgvsName());

    if (!isSnp && repeats.size() == 0 && altAlleles.size() != 1) {
      // in/dels - must have HGVS to represent each change
      throw new IllegalStateException(errorLocation + ": has " + altAlleles.size() + " alt alleles; max is 1");
    }


    Map<String, String> vcfMap = new HashMap<>();
    List<String> missingAlts = new ArrayList<>(altAlleles);
    long vcfPosition = -1;

    if (isSnp) {
      for (String h : hgvsNames) {
        String hgvs = m_refSeqChromosome + ":" + h;
        VcfHelper.VcfData vcf;
        try {
          vcf = vcfHelper.hgvsToVcf(hgvs);
        } catch (ParseException ex) {
          throw new ParseException(errorLocation + " / " + hgvs + ": " + ex.getMessage(), ex);
        }

        if (vcfPosition == -1) {
          vcfPosition = vcf.pos;
        } else if (vcfPosition != vcf.pos) {
          throw new ParseException(errorLocation + ": SNP position mismatch (" + vcfPosition + " vs. " +
              vcf.pos + " for " + hgvs + ")");
        }

        if (!refAllele.equals(vcf.ref)) {
          throw new ParseException(errorLocation + ": VCF's reference allele does not match (" +
              refAllele + " vs. " + vcf.ref + " for " + hgvs + ")");
        }
        if (!missingAlts.remove(vcf.alt)) {
          throw new ParseException(errorLocation + ": VCF's alt allele does not match (expecting " +
              altAlleles + ",  got " + vcf.alt + " for " + hgvs + ")");
        }

        vcfMap.put(refAllele, vcf.ref);
        vcfMap.put(vcf.alt, vcf.alt);
      }
      // warnings
      if (vcfPosition != vl.getCpicPosition()) {
        System.out.println(errorLocation + ": pos/vcf mismatch " + vl.getCpicPosition() +
            " vs. " + vcfPosition);
      }

    } else if (repeats.size() > 0) {
      Map<String, VcfHelper.VcfData> firstPass = new HashMap<>();
      for (String h : hgvsNames) {
        String repeatAlt;
        // treat dups as a form of repeat
        if (h.endsWith("dup")) {
          String repeatedSequence = repeats.get(0);
          repeatedSequence = repeatedSequence.substring(0, repeatedSequence.indexOf("("));
          repeatAlt = repeatedSequence + "(2)";
        } else {
          Matcher m = sf_hgvsRepeatPattern.matcher(h);
          if (!m.matches()) {
            throw new IllegalStateException(errorLocation + ": Invalid HGVS repeat (" + h + ")");
          }
          repeatAlt = m.group(1).replaceAll("\\[", "(").replaceAll("]", ")");
        }
        if (repeatAlt.equals(refAllele)) {
          continue;
        }
        if (!missingAlts.remove(repeatAlt)) {
          throw new IllegalStateException(errorLocation + ": Repeat alt allele does not match (expecting " +
              altAlleles + ",  got " + repeatAlt + ")");
        }

        String hgvs = m_refSeqChromosome + ":" + h;
        VcfHelper.VcfData vcf = vcfHelper.hgvsToVcf(hgvs);
        firstPass.put(repeatAlt, vcf);
      }
      VcfHelper.VcfData vcf = VcfHelper.normalizeRepeats(m_chromosome, firstPass.values());

      vcfPosition = vcf.pos;
      vcfMap.put(refAllele, vcf.ref);
      SortedSet<String> repAlts = new TreeSet<>(vcf.getAlts());
      if (altAlleles.size() != repAlts.size()) {
        throw new IllegalStateException(errorLocation + ": Expected " + altAlleles.size() +
            " repeats but VCF normalization produced " + repAlts.size());
      }
      Iterator<String> repIt = repAlts.iterator();
      for (String alt : altAlleles) {
        vcfMap.put(alt, repIt.next());
      }

    } else {
      // in/del/dup
      for (String h : hgvsNames) {
        String hgvs = m_refSeqChromosome + ":" + h;
        VcfHelper.VcfData vcf = vcfHelper.hgvsToVcf(hgvs);

        String alt;
        Matcher m = sf_hgvsDelPattern.matcher(h);
        if (m.matches()) {
          alt = "del" + refAllele;
        } else if (h.endsWith("dup")) {
          alt = refAllele + refAllele;
        } else if (h.contains("ins")) {
          m = sf_hgvsInsPattern.matcher(h);
          if (!m.matches()) {
            throw new IllegalStateException(errorLocation + ": unsupported ins or delins - " + h);
          }
          alt = m.group(2);
        } else {
          throw new IllegalStateException(errorLocation + ": Unsupported HGVS - " + h);
        }

        vcfPosition = validateVcfPosition(vcfPosition, vcf, errorLocation);
        validateVcfRef(vcfMap, refAllele, vcf, errorLocation);

        if (!missingAlts.remove(alt)) {
          throw new IllegalStateException(errorLocation + ": Alt allele does not match (expecting " + altAlleles +
              ",  got " + alt + ")");
        }
        vcfMap.put(alt, vcf.alt);
      }
    }

    if (missingAlts.size() > 0) {
      if (altAlleles.size() == 1) {
        throw new IllegalStateException(errorLocation + ": Missing alts " + missingAlts);
      } else {
        if (!vcfMap.entrySet().stream().allMatch((e) -> e.getKey().equals(e.getValue()))) {
          // CPIC alleles needs to be translated
          throw new IllegalStateException(errorLocation + ": Don't know how to translate " + missingAlts);
        } else {
          // no translation, use as is
          for (String alt : missingAlts) {
            vcfMap.put(alt, alt);
          }
        }
      }
    }

    vl.setPosition(vcfPosition);
    vl.setCpicToVcfAlleleMap(vcfMap);
    String vcfRef = vcfMap.get(refAllele);
    vl.setRef(vcfRef);
    altAlleles.forEach((a) -> {
      String vcfAlt = vcfMap.get(a);
      vl.addAlt(vcfAlt);
    });
  }

  private long validateVcfPosition(long vcfPosition, VcfHelper.VcfData vcf, String errorLocation) {
    if (vcfPosition != -1 && vcfPosition != vcf.pos) {
      throw new IllegalStateException(errorLocation + ": VCF position mismatch (" + vcfPosition + " vs. " + vcf.pos +
          ")");
    }
    return vcf.pos;
  }

  private void validateVcfRef(Map<String, String> vcfMap, String refAllele, VcfHelper.VcfData vcf,
      String errorLocation) {

    if (vcfMap.containsKey(refAllele)) {
      if (!vcf.alt.equals(vcfMap.get(refAllele))) {
        throw new IllegalStateException(errorLocation + ": VCF ref mismatch (" + vcfMap.get(refAllele) + " vs. " +
            vcf.ref);
      }
    } else {
      vcfMap.put(refAllele, vcf.ref);
    }
  }

  /**
   * Makes sure allele names do not have a "/" in them.
   */
  public void validateAlleleNames() throws IllegalStateException {
    for (NamedAllele na : m_namedAlleles) {
      if (na.getName().contains("/")) {
        throw new IllegalStateException(m_source + " has a " + m_geneSymbol +
            " NamedAllele has invalid name with a '/': " + na.getName());
      }
    }
  }


  /**
   * Get updated {@link NamedAllele} with re-ordered alleles based on re-sorted positions.
   */
  private NamedAllele reorderHaplotypeAlleles(NamedAllele hap, VariantLocus[] oldPositions,
      VariantLocus[] newPositions, String[] fixedAlleles) {

    List<VariantLocus> oldPos = Arrays.stream(oldPositions).toList();
    // resort alleles, cpicAlleles
    String[] alleles = new String[newPositions.length];
    String[] cpicAlleles = new String[newPositions.length];
    for (int x = 0; x < newPositions.length; x += 1) {
      if (oldPositions[x] == newPositions[x]) {
        alleles[x] = fixedAlleles[x];
        cpicAlleles[x] = hap.getCpicAllele(x);
      } else {
        int oldIdx = oldPos.indexOf(newPositions[x]);
        alleles[x] = fixedAlleles[oldIdx];
        cpicAlleles[x] = hap.getCpicAllele(oldIdx);
      }
    }

    return new NamedAllele(hap.getId(), hap.getName(), alleles, cpicAlleles, hap.getMissingPositions(),
        hap.isReference(), hap.getNumCombinations(), hap.getNumPartials());
  }
}
