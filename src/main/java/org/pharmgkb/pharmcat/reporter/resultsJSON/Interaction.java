package org.pharmgkb.pharmcat.reporter.resultsJSON;

import java.util.List;
import java.util.Set;
import java.util.TreeSet;
import com.google.common.collect.Multimap;
import com.google.common.collect.TreeMultimap;
import org.pharmgkb.pharmcat.reporter.model.CPICinteraction;
import org.pharmgkb.pharmcat.reporter.model.Group;
import org.pharmgkb.pharmcat.reporter.model.RelatedChemical;


public class Interaction {

  private String name;

  private List<RelatedChemical> relatedChemicals;
  private String source;
  private String summaryHtml;
  private String textHtml;
  private Set<Group> m_matchingGroups;
  private String m_url;
  private Multimap<String,String> m_matchedDiplotypes = TreeMultimap.create();

  public Interaction( CPICinteraction inter){
    this.name = inter.getName();
    this.relatedChemicals = inter.getRelatedChemicals();
    this.source = inter.getSource();
    this.summaryHtml = inter.getSummaryHtml();
    this.textHtml = inter.getTextHtml();
    m_url = "https://www.pharmgkb.org/guideline/"+inter.getId();
  }

  public String getName() {
    return name;
  }

  public void setName(String name) {
    this.name = name;
  }

  public List<RelatedChemical> getRelatedChemicals() {
    return relatedChemicals;
  }

  public void setRelatedChemicals(List<RelatedChemical> relatedChemicals) {
    this.relatedChemicals = relatedChemicals;
  }

  public String getSource() {
    return source;
  }

  public void setSource(String source) {
    this.source = source;
  }

  public String getSummaryHtml() {
    return summaryHtml;
  }

  public void setSummaryHtml(String summaryHtml) {
    this.summaryHtml = summaryHtml;
  }

  public String getTextHtml() {
    return textHtml;
  }

  public void setTextHtml(String textHtml) {
    this.textHtml = textHtml;
  }

  public Set<Group> getMatchingGroups() {
    return m_matchingGroups;
  }

  public void addMatchingGroup(Group group) {
    if (m_matchingGroups == null) {
      m_matchingGroups = new TreeSet<>();
    }
    m_matchingGroups.add(group);
  }

  public Multimap<String,String> getMatchedDiplotypes() {
    return m_matchedDiplotypes;
  }

  public void putMatchedDiplotype(String id, String diplotype) {
    m_matchedDiplotypes.put(id, diplotype);
  }

  public String getUrl() {
    return m_url;
  }
}
