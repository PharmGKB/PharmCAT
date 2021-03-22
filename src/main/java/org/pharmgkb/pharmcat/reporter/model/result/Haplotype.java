package org.pharmgkb.pharmcat.reporter.model.result;

import com.google.gson.annotations.Expose;
import com.google.gson.annotations.SerializedName;
import org.apache.commons.lang3.ObjectUtils;
import org.apache.commons.lang3.builder.EqualsBuilder;
import org.apache.commons.lang3.builder.HashCodeBuilder;
import org.pharmgkb.common.comparator.HaplotypeNameComparator;


/**
 * Model to represent a haplotype and all derived information
 *
 * @author Ryan Whaley
 */
public class Haplotype implements Comparable<Haplotype> {
  public static final String UNKNOWN = "Unknown";

  @Expose
  @SerializedName("gene")
  private final String f_gene;
  @Expose
  @SerializedName("name")
  private final String f_name;
  @Expose
  @SerializedName("function")
  private String m_function;
  @Expose
  @SerializedName("reference")
  private boolean m_reference = false;

  /**
   * public constructor
   *
   * this can apply some known transformations to the allele name in some cases. the name you get out may not match the
   * name you put in.
   */
  public Haplotype(String gene, String name) {
    f_gene = gene;
    f_name = name;
  }

  /**
   * Gets the name of the gene for this haplotype
   */
  public String getGene() {
    return f_gene;
  }

  /**
   * Gets just the name of this haplotype, e.g. "*10"
   */
  public String getName() {
    return f_name;
  }

  /**
   * Gets the display name of this haplotype, e.g. "GENEX*10"
   */
  public String toString() {
    if (f_name.startsWith("*")) {
      return f_gene + f_name;
    }
    else {
      return f_gene + " " + f_name;
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
    return f_name.equals(UNKNOWN);
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
    if (!(o instanceof Haplotype)) {
      return false;
    }
    if (o == this) {
      return true;
    }

    Haplotype h = (Haplotype)o;
    return new EqualsBuilder()
        .append(f_gene, h.getGene())
        .append(f_name, h.getName())
        .isEquals();
  }

  @Override
  public int hashCode() {
    return new HashCodeBuilder(31, 167)
        .append(f_gene)
        .append(f_name)
        .toHashCode();
  }
}
