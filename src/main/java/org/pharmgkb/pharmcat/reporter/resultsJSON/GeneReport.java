package org.pharmgkb.pharmcat.reporter.resultsJSON;

import java.util.ArrayList;
import java.util.List;
import java.util.Set;
import java.util.SortedSet;
import java.util.TreeSet;
import java.util.stream.Collectors;
import org.pharmgkb.pharmcat.haplotype.model.GeneCall;
import org.pharmgkb.pharmcat.haplotype.model.Variant;
import org.pharmgkb.pharmcat.reporter.model.CPICException;


public class GeneReport {

  private String m_gene;
  private Set<String> m_diplotypes;
  private SortedSet<Variant> m_variants = new TreeSet<>();
  private List<CPICException> m_exceptList = new ArrayList<>();

  public GeneReport(GeneCall call) {
    m_gene = call.getGene();
    m_diplotypes = call.getDiplotypes().stream()
        .map(m -> call.getGene() + ":" + m.getName())
        .collect(Collectors.toSet());
    m_variants.addAll(call.getVariants());
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

  public SortedSet<Variant> getVariants() {
    return m_variants;
  }
}

