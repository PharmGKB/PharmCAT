package org.pharmgkb.pharmcat.reporter.model.result;

import java.util.Collection;
import java.util.List;
import java.util.Objects;
import java.util.Set;
import java.util.TreeSet;
import java.util.stream.Collectors;
import javax.annotation.Nonnull;
import com.google.common.collect.Multimap;
import com.google.common.collect.TreeMultimap;
import org.pharmgkb.pharmcat.reporter.model.DosingGuideline;
import org.pharmgkb.pharmcat.reporter.model.Group;
import org.pharmgkb.pharmcat.reporter.model.RelatedGene;


/**
 * This class is a wrapper around the {@link DosingGuideline} class that also stores the results of the matching code
 * found elsewhere.
 *
 * @author Ryan Whaley
 */
public class GuidelineReport implements Comparable<GuidelineReport> {

  private DosingGuideline m_dosingGuideline;
  private Set<Group> m_matchingGroups;
  private Multimap<String,String> m_matchedDiplotypes = TreeMultimap.create();
  private boolean m_reportable = false;

  public GuidelineReport( DosingGuideline guideline){
    m_dosingGuideline = guideline;
  }

  /**
   * Gets the name of the guideline, a pass-through to the stored guideline.
   */
  public String getName() {
    return m_dosingGuideline.getName();
  }

  /**
   * Gets just the symbols of the related genes of the guideline. Calculated from data in the original guideline.
   */
  public Set<String> getRelatedGeneSymbols() {
    return m_dosingGuideline.getRelatedGenes().stream()
        .map(RelatedGene::getSymbol).collect(Collectors.toSet());
  }

  /**
   * Gets the summary text of the guideline, a pass-through to the stored guideline.
   */
  public String getSummaryHtml() {
    return m_dosingGuideline.getSummaryMarkdown().getHtml();
  }

  /**
   * Gets all annotation groups of the guideline, a pass-through to the stored guideline.
   */
  public List<Group> getGroups() {
    return m_dosingGuideline.getGroups();
  }

  /**
   * Gets only the matching annotation groups based on the called genotypes
   */
  public Set<Group> getMatchingGroups() {
    return m_matchingGroups;
  }

  public void addMatchingGroup(Group group) {
    if (m_matchingGroups == null) {
      m_matchingGroups = new TreeSet<>();
    }
    m_matchingGroups.add(group);
  }

  /**
   * Gets the URL for the whole annotation
   */
  public String getUrl() {
    return "https://www.pharmgkb.org/guideline/" + m_dosingGuideline.getId();
  }

  /**
   * Gets the matched diplotypes for this annotation (should be a subset of all called genotypes)
   */
  public Multimap<String,String> getMatchedDiplotypes() {
    return m_matchedDiplotypes;
  }

  public void putMatchedDiplotype(String id, String diplotype) {
    m_matchedDiplotypes.put(id, diplotype);
  }

  /**
   * Does this annotation have enough information in the called genes to report a specific annotation group?
   */
  public boolean isReportable() {
    return m_reportable;
  }

  /**
   * Sets whether this annotation has enough information in the called genes to report a specific annotation group
   */
  public void setReportable(Collection<String> calledGenes) {
    m_reportable = getRelatedGeneSymbols().stream().allMatch(calledGenes::contains);
  }

  @Override
  public int compareTo(@Nonnull GuidelineReport o) {
    int rez = Boolean.compare(isReportable(), o.isReportable());
    if (rez != 0) {
      return rez * -1;
    }
    rez = Objects.compare(getName(), o.getName(), String.CASE_INSENSITIVE_ORDER);
    if (rez != 0) {
      return rez;
    }
    return 0;
  }
}
