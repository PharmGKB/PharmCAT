package org.pharmgkb.pharmcat.haplotype;

import java.util.List;
import java.util.regex.Matcher;
import javax.annotation.Nonnull;
import javax.annotation.Nullable;
import com.google.common.base.Preconditions;
import org.apache.commons.lang3.ObjectUtils;
import org.apache.commons.lang3.StringUtils;
import org.pharmgkb.common.comparator.ChromosomeNameComparator;
import org.pharmgkb.pharmcat.definition.model.VariantLocus;
import org.pharmgkb.pharmcat.definition.model.VariantType;


/**
 * This class represents a sample's allele at a specific location.
 * <p>
 * It's comparator only takes it's chromosomal position into account.
 *
 * @author Mark Woon
 */
public class SampleAllele implements Comparable<SampleAllele> {
  static boolean STRICT_CHECKING = false;
  private String m_chromosome;
  private int m_position;
  private String m_allele1;
  private String m_allele2;
  private boolean m_isPhased;
  private List<String> m_vcfAlleles;

  public SampleAllele(@Nonnull String chromosome, long position, @Nonnull String a1, @Nullable String a2,
      boolean isPhased, @Nonnull List<String> vcfAlleles) {
    m_chromosome = chromosome;
    m_position = (int)position;
    if (a1.contains("ins") || a1.contains("del")) {
      m_allele1 = a1;
    } else {
      m_allele1 = a1.toUpperCase();
    }
    if (a2 != null) {
      if (a2.contains("ins") || a2.contains("del")) {
        m_allele2 = a2;
      } else {
        m_allele2 = a2.toUpperCase();
      }
    }
    m_isPhased = isPhased;
    m_vcfAlleles = vcfAlleles;
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

  public String getAllele1() {
    return m_allele1;
  }

  public String getAllele2() {
    return m_allele2;
  }
  
  public boolean isMissing() {
    return StringUtils.isEmpty(m_allele1) && StringUtils.isEmpty(m_allele2); 
  }

  public boolean isPhased() {
    return m_isPhased;
  }

  public List<String> getVcfAlleles() {
    return m_vcfAlleles;
  }

  @Override
  public String toString() {
    return m_allele1 + "/" + m_allele2 + " @ " + m_chromosome + ":" + m_position;
  }

  @Override
  public int compareTo(@Nonnull SampleAllele o) {

    int rez = ChromosomeNameComparator.getComparator().compare(m_chromosome, o.getChromosome());
    if (rez != 0) {
      return rez;
    }
    return ObjectUtils.compare(m_position, o.getPosition());
  }


  /**
   * Interprets the alleles in this {@link SampleAllele} in terms of the given {@link VariantLocus}.
   * This will return a <strong>new</strong> {@link SampleAllele} if the {@link VariantLocus} is not a SNP, with it's
   * alleles modified to use the format used by the allele definitions.
   */
  public SampleAllele forVariant(VariantLocus variant) {

    if (variant.getType() == VariantType.SNP) {
      return this;
    }
    String a1 = m_allele1;
    String a2 = m_allele2;

    if (variant.getType() == VariantType.INS) {
      // VCF:         TC  -> TCA
      // definition:  del -> insA
      a1 = convertInsertion(m_allele1);
      a2 = convertInsertion(m_allele2);

    } else if (variant.getType() == VariantType.DEL) {
      // VCF:         TC -> T
      // definition:  C  -> delC
      a1 = convertDeletion(variant, m_allele1);
      a2 = convertDeletion(variant, m_allele2);
    } else if (variant.getType() == VariantType.REPEAT) {
      // VCF:         ATAA -> ATATATATATATATATAA
      // definition:  A(TA)6TAA  -> A(TA)7TAA
      Matcher m = VariantLocus.REPEAT_PATTERN.matcher(variant.getReferenceRepeat());
      if (!m.matches()) {
        throw new IllegalStateException("Invalid repeat format for " + variant.getChromosomeHgvsName());
      }
      String prefix = m.group(1);
      String repeat = m.group(2);
      String postfix = m.group(4);
      a1 = convertRepeat(variant, prefix, repeat, postfix, m_allele1);
      a2 = convertRepeat(variant, prefix, repeat, postfix, m_allele2);
    }
    return new SampleAllele(m_chromosome, m_position, a1, a2, m_isPhased, m_vcfAlleles);
  }


  /**
   * Convert from VCF insertion to allele definition insertion format.
   * <pre><code>
   * VCF:         TC  -> TCA
   * definition:  del -> insA
   * </code></pre>
   */
  private @Nonnull String convertInsertion(@Nonnull String allele) {

    String ref = m_vcfAlleles.get(0);
    if (allele.equals(ref)) {
      return "del";
    }

    // must be an ALT, and therefore longer than REF
    Preconditions.checkState(allele.length() > ref.length(), "Not an insertion: " + ref + " >" + allele);
    return "ins" + allele.substring(ref.length());
  }

  /**
   * Convert from VCF deletion to allele definition deletion format.
   * <pre><code>
   * VCF:         TC -> T
   * definition:  C  -> delC
   * </code></pre>
   */
  private @Nonnull String convertDeletion(@Nonnull VariantLocus variant, @Nonnull String allele) {

    String ref = m_vcfAlleles.get(0);
    if (allele.equals(ref)) {
      return allele.substring(1);
    }

    // must be an ALT, and therefore shorter than REF
    Preconditions.checkState(allele.length() < ref.length(), "Not an deletion: " + ref + " >" + allele + " @ " +
        variant.getChromosomeHgvsName());
    return "del" + ref.substring(1);
  }

  private @Nonnull String convertRepeat(@Nonnull VariantLocus variant, @Nonnull String prefix, @Nonnull String repeat,
      @Nonnull String postfix, @Nonnull String allele) {

    if (allele.contains("(")) {
      // already a repeat
      if (STRICT_CHECKING) {
        Matcher m = VariantLocus.REPEAT_PATTERN.matcher(allele);
        if (!m.matches()) {
          throw new IllegalArgumentException("Sample has " + allele + ", which is an invalid repeat format @ " +
              variant.getChromosomeHgvsName());
        }
        if (!m.group(1).startsWith(repeat) || !m.group(2).endsWith(postfix) || !m.group(4).endsWith(postfix)) {
          throw new IllegalArgumentException("Sample has " + allele + ", which doesn't match expected repeat (" +
              prefix + "(" + repeat + ")#" + postfix + " @ " + variant.getChromosomeHgvsName());
        }
      }
      return allele;
    }

    // validate
    if (!allele.startsWith(prefix) || !allele.endsWith(postfix)) {
      if (STRICT_CHECKING) {
        throw new IllegalArgumentException("Sample has " + allele + ", which doesn't match expected repeat " +
            prefix + "(" + repeat + ")#" + postfix + " @ " + variant.getChromosomeHgvsName());
      }
      return allele;
    }

    String rep = allele.substring(prefix.length(), allele.length() - postfix.length());
    int numReps =  rep.length() / repeat.length();
    if (!rep.equals(StringUtils.repeat(repeat, numReps))) {
      if (STRICT_CHECKING) {
        throw new IllegalArgumentException("Sample has " + allele + ", which doesn't match expected repeat " +
            prefix + "(" + repeat + ")#" + postfix + " @ " + variant.getChromosomeHgvsName());
      } else {
        return allele;
      }
    }

    return prefix + "(" + repeat + ")" + numReps + postfix;
  }
}
