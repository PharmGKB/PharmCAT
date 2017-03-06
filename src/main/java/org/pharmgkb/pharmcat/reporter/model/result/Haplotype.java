package org.pharmgkb.pharmcat.reporter.model.result;

import java.util.Optional;
import javax.annotation.Nonnull;
import org.apache.commons.lang3.ObjectUtils;
import org.pharmgkb.common.comparator.HaplotypeNameComparator;
import org.pharmgkb.pharmcat.definition.IncidentalFinder;


/**
 * Model to represent a haplotype and all derived information
 *
 * @author whaleyr
 */
public class Haplotype implements Comparable<Haplotype> {

  private String m_gene;
  private String m_name;
  private String m_calledFunction;
  private String m_guidelineFunction;
  private boolean m_incidental = false;

  /**
   * public constructor
   *
   * this can apply some known transformations to the allele name in some cases. the name you get out may not match the
   * name you put in.
   */
  public Haplotype(@Nonnull String gene, @Nonnull String name) {
    m_gene = gene;

    switch (m_gene) {
      case "CYP2C19":
        m_name = name.replaceAll("\\*4[AB]", "*4");
        break;
      default:
        m_name = name;
    }
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
   * Uses the {@link HaplotypeNameComparator} class to compare based on allele name
   */
  @Override
  public int compareTo(@Nonnull Haplotype o) {
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
}
