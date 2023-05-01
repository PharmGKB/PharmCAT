package org.pharmgkb.pharmcat.haplotype;

import java.util.HashSet;
import java.util.List;
import java.util.Set;
import com.google.common.base.Preconditions;
import org.apache.commons.lang3.ObjectUtils;
import org.checkerframework.checker.nullness.qual.Nullable;
import org.pharmgkb.common.comparator.ChromosomeNameComparator;


/**
 * This class represents a sample's allele at a specific location.
 * <p>
 * Its comparator only takes its chromosomal position into account.
 *
 * @author Mark Woon
 */
public class SampleAllele implements Comparable<SampleAllele> {
  private final String m_chromosome;
  private final int m_position;
  private final String m_allele1;
  private final String m_allele2;
  private final String m_computedAllele1;
  private final String m_computedAllele2;
  private final boolean m_isPhased;
  private final boolean m_isEffectivelyPhased;
  private final List<String> m_vcfAlleles;
  private final Set<String> m_undocumentedVariations = new HashSet<>();
  private boolean m_treatUndocumentedVariationsAsReference;


  public SampleAllele(String chromosome, long position, String a1, @Nullable String a2, boolean isPhased,
      boolean isEffectivelyPhased, List<String> vcfAlleles, @Nullable Set<String> undocumentedVariations,
      boolean treatUndocumentedAsReference) {
    Preconditions.checkNotNull(vcfAlleles);
    m_chromosome = chromosome;
    m_position = (int)position;

    // must set undocumentedVariations and vcfAlleles before getting computed alleles
    if (undocumentedVariations != null) {
      m_undocumentedVariations.addAll(undocumentedVariations);
    }
    m_vcfAlleles = vcfAlleles;

    m_allele1 = a1.toUpperCase();
    m_computedAllele1 = computeAllele(m_allele1, treatUndocumentedAsReference);
    if (a2 != null) {
      m_allele2 = a2.toUpperCase();
      m_computedAllele2 = computeAllele(m_allele2, treatUndocumentedAsReference);
    } else {
      m_allele2 = null;
      m_computedAllele2 = null;
    }
    m_isPhased = isPhased;
    m_isEffectivelyPhased = isEffectivelyPhased;
  }

  private String computeAllele(String allele, boolean treatUndocumentedAsReference) {
    if (treatUndocumentedAsReference && m_undocumentedVariations.contains(allele)) {
      m_treatUndocumentedVariationsAsReference = true;
      return m_vcfAlleles.get(0);
    }
    return allele;
  }

  /**
   * This constructor is primarily for use in tests only.
   */
  protected SampleAllele(String chromosome, long position, String a1, @Nullable String a2,
      boolean isPhased, List<String> vcfAlleles) {
    this(chromosome, position, a1, a2, isPhased, isPhased, vcfAlleles, null, false);
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
  public String getAllele1() {
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
  public String getComputedAllele1() {
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


  public boolean isHomozygous() {
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
    return m_allele1 + (m_allele2 == null ? "" : "/" + m_allele2) + " @ " + m_chromosome + ":" + m_position;
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
