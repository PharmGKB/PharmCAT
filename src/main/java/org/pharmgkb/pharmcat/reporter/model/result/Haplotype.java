package org.pharmgkb.pharmcat.reporter.model.result;

import java.util.Optional;
import org.apache.commons.lang3.ObjectUtils;
import org.apache.commons.lang3.builder.EqualsBuilder;
import org.apache.commons.lang3.builder.HashCodeBuilder;
import org.pharmgkb.common.comparator.HaplotypeNameComparator;
import org.pharmgkb.pharmcat.definition.IncidentalFinder;


/**
 * Model to represent a haplotype and all derived information
 *
 * @author Ryan Whaley
 */
public class Haplotype implements Comparable<Haplotype> {

  private String m_gene;
  private String m_name;
  private String m_calledFunction;
  private String m_guidelineFunction;
  private boolean m_incidental = false;
  private boolean m_reference = false;

  /**
   * public constructor
   *
   * this can apply some known transformations to the allele name in some cases. the name you get out may not match the
   * name you put in.
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
   * Gets a key for this haplotype used to lookup matching annotation groups. This integrates known per-gene
   * modifications that need to be done to allele calls before lookup.
   */
  public String printLookup() {
    switch (m_gene) {
      case "CFTR":
        if (m_name.equals("Reference") || m_incidental) {
          return "Other";
        }
        return m_name;
      case "DPYD":
        if (m_name.equals("Reference")) {
          return "Any normal function variant or no variant detected";
        }
      default:
        return m_name;
    }
  }

  public String printDisplay() {
    return m_name;
  }

  /**
   * Gets the function for this haplotype that the caller specified
   */
  private Optional<String> getCalledFunction() {
    return Optional.ofNullable(m_calledFunction);
  }

  public void setCalledFunction(String calledFunction) {
    m_calledFunction = calledFunction;
  }

  /**
   * Gets the function for this haplotype that the guideline specified
   */
  private String getGuidelineFunction() {
    return m_guidelineFunction;
  }

  public void setGuidelineFunction(String guidelineFunction) {
    m_guidelineFunction = guidelineFunction;
  }

  /**
   * Gets the preferred function for this allele.
   *
   * Defaults to called function but will fall back to guideline function if called is null
   *
   * @return the called function, or the guideline function if that's not present
   */
  public String getFunction() {

    return getCalledFunction().orElse(getGuidelineFunction());
  }

  /**
   * Flags whether this allele is an incidental finding
   * @return true if this allele is an incidental finding, false otherwise
   */
  public boolean isIncidental() {
    return m_incidental;
  }

  public void setIncidental(IncidentalFinder incidental) {
    m_incidental = incidental.isFinding(this);
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
   * Uses the {@link HaplotypeNameComparator} class to compare based on allele name
   */
  @Override
  public int compareTo(Haplotype o) {
    int rez = ObjectUtils.compare(getGene(), o.getGene());
    if (rez != 0) {
      return rez;
    }

    rez = HaplotypeNameComparator.getComparator().compare(getName(), o.getName());
    if (rez != 0) {
      return rez;
    }
    return 0;
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
