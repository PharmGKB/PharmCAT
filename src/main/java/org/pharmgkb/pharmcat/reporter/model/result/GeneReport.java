package org.pharmgkb.pharmcat.reporter.model.result;

import java.util.ArrayList;
import java.util.Collection;
import java.util.List;
import java.util.Objects;
import java.util.Set;
import java.util.SortedSet;
import java.util.TreeSet;
import java.util.stream.Collectors;
import javax.annotation.Nonnull;
import com.google.common.collect.ImmutableList;
import com.google.common.collect.ImmutableSet;
import org.pharmgkb.pharmcat.haplotype.MatchData;
import org.pharmgkb.pharmcat.haplotype.NamedAlleleMatcher;
import org.pharmgkb.pharmcat.haplotype.model.GeneCall;
import org.pharmgkb.pharmcat.haplotype.model.Variant;
import org.pharmgkb.pharmcat.reporter.model.AstrolabeCall;
import org.pharmgkb.pharmcat.reporter.model.DrugLink;
import org.pharmgkb.pharmcat.reporter.model.MatchLogic;
import org.pharmgkb.pharmcat.reporter.model.PharmcatException;


/**
 * This class is used to help collect Gene-related data for later reporting
 */
public class GeneReport implements Comparable<GeneReport> {
  private static final String UNCALLED = "not called";
  public static final String NA = "N/A";

  private String m_gene;
  private Set<String> m_uncalledHaplotypes;
  private SortedSet<Variant> m_variants = new TreeSet<>();
  private List<PharmcatException> m_exceptList = new ArrayList<>();
  private MatchData m_matchData;
  private List<DrugLink> m_relatedDrugs = new ArrayList<>();
  private boolean m_astrolabeCall = false;
  private List<Diplotype> m_diplotypes = new ArrayList<>();

  /**
   * public constructor
   */
  public GeneReport(@Nonnull String geneSymbol) {
    m_gene = geneSymbol;
  }

  /**
   * Sets data in this report based on data found in a {@link GeneCall}
   * @param call a {@link GeneCall} that has been made by the {@link NamedAlleleMatcher}
   */
  public void setCallData(@Nonnull GeneCall call) {
    m_gene = call.getGene();
    m_variants.addAll(call.getVariants());
    m_matchData = call.getMatchData();
    m_uncalledHaplotypes = call.getUncallableHaplotypes();
  }

  /**
   * Sets data in this report based on data found in a {@link AstrolabeCall}
   * @param call a {@link AstrolabeCall} to pull data from
   */
  public void setAstrolabeData(@Nonnull AstrolabeCall call) {
    m_astrolabeCall = true;
    m_gene = call.getGene();
  }

  public List<Diplotype> getDiplotypes() {
    return m_diplotypes;
  }

  public void addDiplotype(Diplotype diplotype) {
    m_diplotypes.add(diplotype);
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

  public void applyExceptions(Collection<PharmcatException> exceptions) {
    if (exceptions != null) {
      exceptions.forEach(this::addException);
    }

    // add incidental allele message if present
    if (getDiplotypes().stream().anyMatch(Diplotype::isIncidental)) {
      PharmcatException pe = new PharmcatException();
      pe.setName("Incidental Finding");
      pe.setMessage("CPIC does not provide recommendations for this genotype. Variant is included in CPIC " +
          "supplemental table of variants recommended by the American College of Medical Genetics (ACMG) that should " +
          "be tested to determine carrier status as a part of population screening programs.");
      m_exceptList.add(pe);
    }
  }

  /**
   * Evaluates an exception and adds it to the report if it's relevent.
   * @param except an exception to be evaluated and possibly added
   */
  private void addException(PharmcatException except) {
    MatchLogic match = except.getMatches();

    if (match.getGene().equals(m_gene)) {
      boolean critHap = match.getHapsCalled().isEmpty()
          || match.getHapsCalled().stream().anyMatch(h -> m_diplotypes.stream().anyMatch(d -> d.hasAllele(h)));
      boolean critUnmatchedHap = match.getHapsMissing().isEmpty()
          || match.getHapsMissing().stream().allMatch(h -> m_uncalledHaplotypes.contains(h));
      boolean critDip = match.getDips().isEmpty()
          || match.getDips().stream().allMatch(d -> m_diplotypes.stream().anyMatch(e -> e.printBare().equals(d)));
      boolean critMissVariant = match.getVariantsMissing().isEmpty() ||
          match.getVariantsMissing().stream().allMatch(v -> m_variants.isEmpty() || m_variants.stream().noneMatch(a -> a.getRsid().equals(v)));

      if (critHap && critUnmatchedHap && critDip && critMissVariant) {
        m_exceptList.add(except);
      }
    }
  }

  /**
   * The variants that constitute this gene
   */
  public SortedSet<Variant> getVariants() {
    return m_variants;
  }

  /**
   * Gets the match data from the gene call itself
   */
  public MatchData getMatchData() {
    return m_matchData;
  }

  /**
   * Gets the Set of Haplotypes the the matcher could not evaluate
   */
  public Set<String> getUncalledHaplotypes() {
    return m_uncalledHaplotypes;
  }

  public String toString() {
    return m_gene + " Report";
  }

  /**
   * True if this gene has entries in the <code>diplotypes</code> property called "uncalled", false if it has a call
   */
  public boolean isCalled() {
    return m_diplotypes != null && m_diplotypes.size() > 0;
  }

  @Override
  public int compareTo(@Nonnull GeneReport o) {
    int rez = Objects.compare(getGene(), o.getGene(), String.CASE_INSENSITIVE_ORDER);
    if (rez != 0) {
      return rez;
    }
    return 0;
  }

  /**
   * Gets a list of {@link DrugLink} objects that are in the same guidelines as this gene
   * @return a list of {@link DrugLink} objects
   */
  public List<DrugLink> getRelatedDrugs() {
    return m_relatedDrugs;
  }

  /**
   * Adds the drugs in the given <code>guideline</code> to this report as {@link DrugLink} objects
   * @param guideline a GuidelineReport with relatedDrugs
   */
  public void addRelatedDrugs(GuidelineReport guideline) {

    guideline.getRelatedDrugs().stream()
        .map(d -> new DrugLink(d, guideline.getId(), guideline.isRxChange()))
        .forEach(m -> m_relatedDrugs.add(m));
  }

  /**
   * True if the genotyping data in the report comes from astrolabe, false if from somewhere else
   */
  public boolean isAstrolabeCall() {
    return m_astrolabeCall;
  }

  /**
   * Gets a Set of the diplotype keys (Strings) that should be used when looking up matching annotation groups.
   * This method ensures all per-gene adjustments and fixes have been applied before it's time to match.
   * @return a Set of diplotype Strings in the form "GENE:
   */
  public Set<String> getDiplotypeLookupKeys() {
    if (!isCalled()) {

      // only HLA-B uses Other/Other, all else is Unknown/Unknown
      String defaultUncalled = getGene().equals("HLA-B") ? "Other/Other" : "Unknown/Unknown";

      return ImmutableSet.of(m_gene + ":" + defaultUncalled);
    }
    else {
      return getDiplotypes().stream().map(Diplotype::printLookupKey).collect(Collectors.toSet());
    }
  }

  public Collection<String> printDisplayCalls() {
    if (!isCalled()) {
      return ImmutableList.of(UNCALLED);
    }

    return getDiplotypes().stream().map(Diplotype::printBare).collect(Collectors.toList());
  }

  public Collection<String> printDisplayFunctions() {
    if (!isCalled()) {
      return ImmutableList.of(NA);
    }
    return getDiplotypes().stream().map(Diplotype::printFunctionPhrase).distinct().collect(Collectors.toSet());
  }

  public Collection<String> printDisplayPhenotypes() {
    if (!isCalled()) {
      return ImmutableList.of(NA);
    }
    return getDiplotypes().stream().map(Diplotype::getPhenotype).distinct().collect(Collectors.toSet());
  }
}

