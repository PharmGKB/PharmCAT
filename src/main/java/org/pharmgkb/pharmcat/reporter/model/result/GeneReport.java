package org.pharmgkb.pharmcat.reporter.model.result;

import java.util.ArrayList;
import java.util.List;
import java.util.Objects;
import java.util.Set;
import java.util.SortedSet;
import java.util.TreeSet;
import javax.annotation.Nonnull;
import com.google.common.collect.ImmutableSortedSet;
import org.pharmgkb.pharmcat.haplotype.MatchData;
import org.pharmgkb.pharmcat.haplotype.NamedAlleleMatcher;
import org.pharmgkb.pharmcat.haplotype.model.GeneCall;
import org.pharmgkb.pharmcat.haplotype.model.Variant;
import org.pharmgkb.pharmcat.reporter.model.AstrolabeCall;
import org.pharmgkb.pharmcat.reporter.model.DrugLink;
import org.pharmgkb.pharmcat.reporter.model.PharmcatException;


/**
 * This class is used to help collect Gene-related data for later reporting
 */
public class GeneReport implements Comparable<GeneReport> {
  private static final String UNCALLED = "uncalled";

  private String m_gene;
  private Set<String> m_diplotypes = new TreeSet<>();
  private Set<String> m_uncalledHaplotypes;
  private SortedSet<Variant> m_variants = new TreeSet<>();
  private List<PharmcatException> m_exceptList = new ArrayList<>();
  private MatchData m_matchData;
  private List<DrugLink> m_relatedDrugs = new ArrayList<>();
  private SortedSet<String> m_functions = new TreeSet<>();

  /**
   * public constructor
   */
  public GeneReport(@Nonnull String geneSymbol) {
    m_gene = geneSymbol;
    addDip(UNCALLED);
  }

  /**
   * public constructor
   * @param call a {@link GeneCall} that has been made by the {@link NamedAlleleMatcher}
   */
  public void setCallData(@Nonnull GeneCall call) {
    m_gene = call.getGene();
    if (!call.getDiplotypes().isEmpty()) {
      m_diplotypes.clear();
      call.getDiplotypes().forEach(d -> addDip(d.getName()));
    }
    m_variants.addAll(call.getVariants());
    m_matchData = call.getMatchData();
    m_uncalledHaplotypes = call.getUncallableHaplotypes();
    if (call.getDiplotypes() != null && call.getDiplotypes().size() > 0) {
      call.getDiplotypes().forEach(d -> m_functions.add(d.getFunction()));
    }
  }

  public void setAstrolabeData(@Nonnull AstrolabeCall call) {
    m_gene = call.getGene();
    if (call.getDiplotypes() != null) {
      m_diplotypes.clear();
      call.getDiplotypes().forEach(this::addDip);
    }
  }

  /**
   * Diplotypes that have been called for this gene
   */
  public Set<String> getDips(){
    return m_diplotypes;
  }

  private void addDip(String dip) {
    String cleanDip = m_gene.equals("CYP2C19") ? dip.replaceAll("\\*4[AB]", "*4") : dip;

    m_diplotypes.add(cleanDip);
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
  public List<PharmcatException> getExceptionList(){
    return m_exceptList;
  }

  public void addException( PharmcatException except ){
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

  public SortedSet<String> getFunctions() {
    if (!isCalled()) {
      return ImmutableSortedSet.of(UNCALLED);
    }
    return m_functions;
  }

  public String toString() {
    return m_gene + " Report";
  }

  public boolean isCalled() {
    return m_diplotypes != null && m_diplotypes.size() > 0 && !m_diplotypes.contains(UNCALLED);
  }

  @Override
  public int compareTo(@Nonnull GeneReport o) {
    int rez = Objects.compare(getGene(), o.getGene(), String.CASE_INSENSITIVE_ORDER);
    if (rez != 0) {
      return rez;
    }
    return 0;
  }

  public List<DrugLink> getRelatedDrugs() {
    return m_relatedDrugs;
  }

  public void addRelatedDrugs(GuidelineReport guideline) {

    guideline.getRelatedDrugs().stream()
        .map(d -> new DrugLink(d, guideline.getId(), guideline.isRxChange()))
        .forEach(m -> m_relatedDrugs.add(m));
  }
}

