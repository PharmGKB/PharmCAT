package org.cpic.haplotype;

import java.util.regex.Pattern;
import javax.annotation.Nonnull;
import org.apache.commons.lang3.ObjectUtils;


/**
 * This class represents a sample's allele at a specific location.
 * <p>
 * It's comparator only takes it's chromosomal position into account.
 *
 * @author Mark Woon
 */
public class SampleAllele implements Comparable<SampleAllele> {
  private static final Pattern sf_chrNumPattern = Pattern.compile("^\\d+$");
  private String m_chromosome;
  private int m_position;
  private String m_allele1;
  private String m_allele2;
  private boolean m_isPhased;

  public SampleAllele(@Nonnull String chromosome, long position, @Nonnull String a1, String a2, boolean isPhased) {
    m_chromosome = chromosome;
    m_position = (int)position;
    m_allele1 = a1;
    m_allele2 = a2;
    m_isPhased = isPhased;
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

  public boolean isPhased() {
    return m_isPhased;
  }

  @Override
  public String toString() {
    return m_allele1 + "/" + m_allele2 + " @ " + m_chromosome + ":" + m_position;
  }

  @Override
  public int compareTo(@Nonnull SampleAllele o) {

    String chr1 = m_chromosome.substring(3);
    String chr2 = o.getChromosome().substring(3);

    if (sf_chrNumPattern.matcher(chr1).matches()) {
      if (sf_chrNumPattern.matcher(chr2).matches()) {
        int rez = ObjectUtils.compare(Integer.parseInt(chr1), Integer.parseInt(chr2));
        if (rez != 0) {
          return rez;
        }
      } else {
        // numeric vs. non-numeric
        return -1;
      }
    } else {
      if (sf_chrNumPattern.matcher(chr2).matches()) {
        // non-numeric vs. numeric
        return 1;
      } else {
        int rez = ObjectUtils.compare(chr1, chr2);
        if (rez != 0) {
          return rez;
        }
      }
    }
    return ObjectUtils.compare(m_position, o.getPosition());
  }
}
