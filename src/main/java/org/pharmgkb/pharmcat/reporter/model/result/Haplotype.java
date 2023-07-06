package org.pharmgkb.pharmcat.reporter.model.result;

import com.google.gson.annotations.Expose;
import com.google.gson.annotations.SerializedName;
import org.apache.commons.lang3.ObjectUtils;
import org.apache.commons.lang3.builder.EqualsBuilder;
import org.apache.commons.lang3.builder.HashCodeBuilder;
import org.pharmgkb.pharmcat.util.HaplotypeNameComparator;

import static org.pharmgkb.pharmcat.reporter.TextConstants.isUnspecified;


/**
 * Model to represent a haplotype and all derived information
 *
 * @author Ryan Whaley
 */
public class Haplotype implements Comparable<Haplotype> {
  public static final String UNKNOWN = "Unknown";

  @Expose
  @SerializedName("gene")
  private final String m_gene;
  @Expose
  @SerializedName("name")
  private final String m_name;
  @Expose
  @SerializedName("function")
  private String m_function;
  @Expose
  @SerializedName("reference")
  private boolean m_reference = false;
  @Expose
  @SerializedName("activityValue")
  private String m_activityValue;

  /**
   * public constructor
   */
  public Haplotype(String gene, String name) {
    m_gene = gene;
    m_name = name;
  }

  /**
   * Gets the name of the gene for this haplotype
   */
  public String getGene() {
    return m_gene;
  }

  /**
   * Gets just the name of this haplotype, e.g. "*10"
   */
  public String getName() {
    return m_name;
  }

  /**
   * Gets the display name of this haplotype, e.g. "GENEX*10"
   */
  public String toString() {
    if (m_name.startsWith("*")) {
      return m_gene + m_name;
    }
    else {
      return m_gene + " " + m_name;
    }
  }

  /**
   * Gets the function for this allele.
   *
   * @return the function
   */
  public String getFunction() {
    return m_function;
  }

  public void setFunction(String function) {
    m_function = function;
  }

  /**
   * Whether this haplotype is the "reference" haplotype for its gene or not
   * @return true if this haplotype is a reference haplotype
   */
  public boolean isReference() {
    return m_reference;
  }

  public void setReference(boolean reference) {
    m_reference = reference;
  }


  /**
   * True if this haplotype is unknown
   */
  public boolean isUnknown() {
    return m_name.equals(UNKNOWN);
  }


  public String getActivityValue() {
    return m_activityValue;
  }

  public void setActivityValue(String activityValue) {
    m_activityValue = activityValue;
  }

  public String toFormattedFunction() {
    if (isUnspecified(m_activityValue)) {
      return m_function;
    } else {
      return String.format("%s (%s)", m_activityValue, m_function);
    }
  }

  /**
   * Uses the {@link HaplotypeNameComparator} class to compare based on allele name
   */
  @Override
  public int compareTo(Haplotype o) {
    int rez = ObjectUtils.compare(getGene(), o.getGene());
    if (rez != 0) {
      return rez;
    }

    return HaplotypeNameComparator.getComparator().compare(getName(), o.getName());
  }

  @Override
  public boolean equals(Object o) {
    if (o == this) {
      return true;
    }
    if (!(o instanceof Haplotype h)) {
      return false;
    }
    return new EqualsBuilder()
        .append(m_gene, h.getGene())
        .append(m_name, h.getName())
        .isEquals();
  }

  @Override
  public int hashCode() {
    return new HashCodeBuilder(31, 167)
        .append(m_gene)
        .append(m_name)
        .toHashCode();
  }
}
