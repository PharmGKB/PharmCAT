package org.pharmgkb.pharmcat.reporter.model.result;

import java.util.ArrayList;
import java.util.List;
import java.util.Set;
import java.util.SortedSet;
import java.util.TreeSet;
import java.util.stream.Collectors;
import org.pharmgkb.pharmcat.haplotype.MatchData;
import org.pharmgkb.pharmcat.haplotype.NamedAlleleMatcher;
import org.pharmgkb.pharmcat.haplotype.model.GeneCall;
import org.pharmgkb.pharmcat.haplotype.model.Variant;
import org.pharmgkb.pharmcat.reporter.model.CPICException;


/**
 * This class is used to help collect Gene-related data for later reporting
 */
public class GeneReport {

  private String m_gene;
  private Set<String> m_diplotypes;
  private Set<String> m_uncalledHaplotypes;
  private SortedSet<Variant> m_variants = new TreeSet<>();
  private List<CPICException> m_exceptList = new ArrayList<>();
  private MatchData m_matchData;

  /**
   * public constructor
   * @param call a {@link GeneCall} that has been made by the {@link NamedAlleleMatcher}
   */
  public GeneReport(GeneCall call) {
    m_gene = call.getGene();
    m_diplotypes = call.getDiplotypes().stream()
        .map(m -> call.getGene() + ":" + m.getName())
        .collect(Collectors.toSet());
    m_variants.addAll(call.getVariants());
    m_matchData = call.getMatchData();
    m_uncalledHaplotypes = call.getUncallableHaplotypes();
  }

  /**
   * Diplotypes that have been called for this gene
   */
  public Set<String> getDips(){
    return m_diplotypes;
  }

  /**
   * The gene symbol for this gene
   */
  public String getGene(){
    return m_gene;
  }

  /**
   * The list of exceptions that apply to this gene
   */
  public List<CPICException> getExceptionList(){
    return m_exceptList;
  }

  public void addException( CPICException except ){
    m_exceptList.add(except);
  }

  /**
   * The variants that constitute this gene
   */
  public SortedSet<Variant> getVariants() {
    return m_variants;
  }

  public MatchData getMatchData() {
    return m_matchData;
  }

  public Set<String> getUncalledHaplotypes() {
    return m_uncalledHaplotypes;
  }
}

