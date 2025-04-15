package org.pharmgkb.pharmcat.haplotype;

import java.util.HashSet;
import java.util.List;
import java.util.Set;
import com.google.common.base.Preconditions;
import com.google.gson.annotations.Expose;
import com.google.gson.annotations.SerializedName;
import org.apache.commons.lang3.ObjectUtils;
import org.jspecify.annotations.Nullable;
import org.pharmgkb.common.comparator.ChromosomeNameComparator;


/**
 * This class represents a sample's allele at a specific location.
 * <p>
 * Its comparator only takes its chromosomal position into account.
 *
 * @author Mark Woon
 */
public class SampleAllele implements Comparable<SampleAllele> {
  @Expose
  @SerializedName("chr")
  private final String m_chromosome;
  @Expose
  @SerializedName("position")
  private final int m_position;
  @Expose
  @SerializedName("allele1")
  private final @Nullable String m_allele1;
  @Expose
  @SerializedName("allele2")
  private final @Nullable String m_allele2;
  @Expose
  @SerializedName("computedAllele1")
  private final @Nullable String m_computedAllele1;
  @Expose
  @SerializedName("computedAllele2")
  private final @Nullable String m_computedAllele2;
  @Expose
  @SerializedName("phased")
  private final boolean m_isPhased;
  @Expose
  @SerializedName("effectivelyPhased")
  private final boolean m_isEffectivelyPhased;
  @Expose
  @SerializedName("vcfAlleles")
  private final List<String> m_vcfAlleles;
  @Expose
  @SerializedName("gt")
  private final String m_gt;
  @Expose
  @SerializedName("vcfCall")
  private final String m_vcfCall;
  @Expose
  @SerializedName("undocumentedVariations")
  private final Set<String> m_undocumentedVariations = new HashSet<>();
  @Expose
  @SerializedName("treatUndocumentedVariationsAsReference")
  private boolean m_treatUndocumentedVariationsAsReference;
  @Expose
  @SerializedName("phaseSet")
  private @Nullable Integer m_phaseSet;


  public SampleAllele(String chromosome, long position, @Nullable String a1, @Nullable String a2, boolean isPhased,
      boolean isEffectivelyPhased, @Nullable Integer phaseSet, List<String> vcfAlleles, String gt,
      @Nullable Set<String> undocumentedVariations, boolean treatUndocumentedAsReference) {
    Preconditions.checkNotNull(vcfAlleles);
    m_chromosome = chromosome;
    m_position = (int)position;

    // must set undocumentedVariations and vcfAlleles before getting computed alleles
    if (undocumentedVariations != null) {
      m_undocumentedVariations.addAll(undocumentedVariations);
    }
    m_vcfAlleles = vcfAlleles;
    m_gt = gt;

    StringBuilder callBuilder = new StringBuilder();
    if (a1 != null) {
      m_allele1 = a1.toUpperCase();
      m_computedAllele1 = computeAllele(m_allele1, treatUndocumentedAsReference);
      callBuilder.append(m_allele1);
    } else {
      m_allele1 = null;
      m_computedAllele1 = ".";
      callBuilder.append(".");
    }
    if (a2 != null) {
      m_allele2 = a2.toUpperCase();
      m_computedAllele2 = computeAllele(m_allele2, treatUndocumentedAsReference);
      if (isPhased) {
        callBuilder.append("|");
      } else {
        callBuilder.append("/");
      }
      callBuilder.append(m_allele2);
    } else {
      m_allele2 = null;
      String[] gtArray = VcfReader.GT_DELIMITER.split(gt);
      if (gtArray.length > 1 && gtArray[1].equals(".")) {
        if (isPhased) {
          callBuilder.append("|");
        } else {
          callBuilder.append("/");
        }
        callBuilder.append(".");
        m_computedAllele2 = ".";
      } else {
        m_computedAllele2 = null;
      }
    }
    m_vcfCall = callBuilder.toString();
    m_isPhased = isPhased;
    m_isEffectivelyPhased = isEffectivelyPhased;
    if (m_isPhased) {
      m_phaseSet = phaseSet;
    }
  }

  private String computeAllele(String allele, boolean treatUndocumentedAsReference) {
    if (treatUndocumentedAsReference && m_undocumentedVariations.contains(allele)) {
      m_treatUndocumentedVariationsAsReference = true;
      return m_vcfAlleles.get(0);
    }
    return allele;
  }

  /**
   * This constructor should only be use in tests.
   */
  protected SampleAllele(String chromosome, long position, @Nullable String a1, @Nullable String a2,
      boolean isPhased, List<String> vcfAlleles, String gt) {
    this(chromosome, position, a1, a2, isPhased, isPhased, null, vcfAlleles, gt, null, false);
  }


  public String getChromosome() {
    return m_chromosome;
  }

  public int getPosition() {
    return m_position;
  }

  public String getChrPosition() {
    return m_chromosome + ":" + m_position;
  }

  /**
   * Gets the first allele from the VCF (per GT column).
   */
  public @Nullable String getAllele1() {
    return m_allele1;
  }

  /**
   * Gets the second allele from the VCF (per GT column).
   */
  public @Nullable String getAllele2() {
    return m_allele2;
  }


  /**
   * Gets the first allele, taking treatUndocumentedVariationsAsReference into account.
   */
  public @Nullable String getComputedAllele1() {
    return m_computedAllele1;
  }

  /**
   * Gets the second allele, taking treatUndocumentedVariationsAsReference into account.
   */
  public @Nullable String getComputedAllele2() {
    return m_computedAllele2;
  }


  /**
   * Gets whether the sample is phased at this location according to VCF.
   */
  public boolean isPhased() {
    return m_isPhased;
  }

  /**
   * Gets the PS value from the VCF.
   */
  public @Nullable Integer getPhaseSet() {
    return m_phaseSet;
  }

  /**
   * Gets whether the sample should be treated as phased at this location (i.e. really phased or homozygous).
   */
  public boolean isEffectivelyPhased() {
    return m_isEffectivelyPhased;
  }

  /**
   * Gets the list of all possible alleles from the VCF (i.e. REF and ALT).
   */
  public List<String> getVcfAlleles() {
    return m_vcfAlleles;
  }

  /**
   * Gets the GT value from the VCF.
   */
  public String getGt() {
    return m_gt;
  }

  /**
   * Gets the alleles separated by VCF phasing delimiter (i.e. "/" or "|").
   * Missing allele will be represented by ".".
   */
  public String getVcfCall() {
    return m_vcfCall;
  }

  public boolean isHomozygous() {
    if (m_computedAllele1 == null) {
      // both alleles can never be null
      Preconditions.checkState(m_computedAllele2 != null);
      return false;
    }
    return m_computedAllele1.equals(m_computedAllele2);
  }


  /**
   * Gets undocumented variations for sample at this position.
   */
  public Set<String> getUndocumentedVariations() {
    return m_undocumentedVariations;
  }

  public boolean isTreatUndocumentedVariationsAsReference() {
    return m_treatUndocumentedVariationsAsReference;
  }


  @Override
  public String toString() {
    return m_vcfCall + " @ " + m_chromosome + ":" + m_position;
  }

  @Override
  public int compareTo(SampleAllele o) {

    int rez = ChromosomeNameComparator.getComparator().compare(m_chromosome, o.getChromosome());
    if (rez != 0) {
      return rez;
    }
    rez = ObjectUtils.compare(m_position, o.getPosition());
    if (rez != 0) {
      return rez;
    }
    rez = ObjectUtils.compare(m_allele1, o.getAllele1());
    if (rez != 0) {
      return rez;
    }
    return ObjectUtils.compare(m_allele2, o.getAllele2());
  }
}
