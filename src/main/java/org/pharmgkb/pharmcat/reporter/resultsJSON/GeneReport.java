package org.pharmgkb.pharmcat.reporter.resultsJSON;

import java.util.ArrayList;
import java.util.List;
import java.util.Set;
import java.util.stream.Collectors;
import org.pharmgkb.pharmcat.haplotype.model.json.GeneCall;
import org.pharmgkb.pharmcat.reporter.model.CPICException;


public class GeneReport {

  private String m_gene;
  private Set<String> m_diplotypes;
  private List<CPICException> m_exceptList = new ArrayList<>();

  public GeneReport(GeneCall call) {
    m_gene = call.getGene();
    m_diplotypes = call.getDiplotypes().stream()
        .map(m -> call.getGene() + ":" + m.getName())
        .collect(Collectors.toSet());
  }

  public Set<String> getDips(){
    return m_diplotypes;
  }

  public String getGene(){
    return m_gene;
  }

  public List<CPICException> getExceptionList(){
    return m_exceptList;
  }

  public void addException( CPICException except ){
    m_exceptList.add(except);
  }
}

