package org.pharmgkb.pharmcat.reporter.model.result;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.List;
import java.util.Objects;
import java.util.Optional;
import java.util.Set;
import java.util.SortedSet;
import java.util.TreeSet;
import java.util.stream.Collectors;
import javax.annotation.Nonnull;
import com.google.common.collect.ImmutableList;
import com.google.common.collect.ImmutableSet;
import org.apache.commons.lang3.StringUtils;
import org.pharmgkb.common.comparator.HaplotypeNameComparator;
import org.pharmgkb.pharmcat.haplotype.NamedAlleleMatcher;
import org.pharmgkb.pharmcat.haplotype.model.GeneCall;
import org.pharmgkb.pharmcat.reporter.DiplotypeFactory;
import org.pharmgkb.pharmcat.reporter.VariantReportFactory;
import org.pharmgkb.pharmcat.reporter.model.AstrolabeCall;
import org.pharmgkb.pharmcat.reporter.model.DrugLink;
import org.pharmgkb.pharmcat.reporter.model.MessageAnnotation;
import org.pharmgkb.pharmcat.reporter.model.MessageVariant;
import org.pharmgkb.pharmcat.reporter.model.VariantReport;
import org.pharmgkb.pharmcat.util.Ugt1a1AlleleMatcher;


/**
 * This class is used to help collect Gene-related data for later reporting
 */
public class GeneReport implements Comparable<GeneReport> {
  private static final List<String> sf_reducibleGeneCalls = ImmutableList.of("UGT1A1");
  private static final Set<String> sf_overrideDiplotypes = ImmutableSet.of("SLCO1B1");
  private static final String UNCALLED = "not called";
  public static  final String NA = "N/A";

  private String m_gene;
  private SortedSet<String> m_uncalledHaplotypes;
  private List<MessageAnnotation> m_messages = new ArrayList<>();
  private List<DrugLink> m_relatedDrugs = new ArrayList<>();
  private boolean m_astrolabeCall = false;
  private List<Diplotype> m_matcherDiplotypes = new ArrayList<>();
  private List<Diplotype> m_reporterDiplotypes = new ArrayList<>();
  private List<VariantReport> m_variantReports = new ArrayList<>();
  private boolean m_phased = false;
  private List<String> m_highlightedVariants = new ArrayList<>();

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
  public void setCallData(@Nonnull GeneCall call) throws IOException {
    m_gene = call.getGene();
    m_uncalledHaplotypes = new TreeSet<>(HaplotypeNameComparator.getComparator());
    m_uncalledHaplotypes.addAll(call.getUncallableHaplotypes());
    m_phased = call.isPhased();

    VariantReportFactory variantReportFactory = new VariantReportFactory(m_gene);
    call.getVariants().stream()
        .map(variantReportFactory::make).forEach(a -> m_variantReports.add(a));
    call.getMatchData().getMissingPositions().stream()
        .map(variantReportFactory::make).forEach(a -> m_variantReports.add(a));
    call.getVariantsOfInterest().stream()
        .map(variantReportFactory::make).forEach(a -> m_variantReports.add(a));
  }

  /**
   * Sets data in this report based on data found in a {@link AstrolabeCall}
   * @param call a {@link AstrolabeCall} to pull data from
   */
  public void setAstrolabeData(@Nonnull AstrolabeCall call) {
    m_astrolabeCall = true;
    m_gene = call.getGene();
  }

  public void setDiplotypes(DiplotypeFactory diplotypeFactory, GeneCall geneCall) {
    m_matcherDiplotypes.addAll(diplotypeFactory.makeDiplotypes(geneCall));

    // for UGT1A1 we need to calculate diplotypes slightly differently
    if (getGene().equals("UGT1A1") && (!isPhased() || geneCall.getDiplotypes().size() > 1)) {
      Set<String> diplotypes = Ugt1a1AlleleMatcher.makeLookupCalls(this);
      m_reporterDiplotypes.addAll(diplotypeFactory.makeDiplotypes(diplotypes));
    } else {
      // some genes disregard the actual call and calculate the call by the reporter, that logic goes here
      if (sf_overrideDiplotypes.contains(getGene())) {
        m_reporterDiplotypes.addAll(diplotypeFactory.makeOverrideDiplotypes());
      }
      else {
        m_reporterDiplotypes.addAll(diplotypeFactory.makeDiplotypes(geneCall));
      }
    }
  }

  public void setDiplotypes(DiplotypeFactory diplotypeFactory, AstrolabeCall astrolabeCall) {
    diplotypeFactory.makeDiplotypes(astrolabeCall).forEach(d -> {
      m_matcherDiplotypes.add(d);
      m_reporterDiplotypes.add(d);
    });
  }

  /**
   * The gene symbol for this gene
   */
  public String getGene(){
    return m_gene;
  }

  /**
   * The list of messages that apply to this gene
   */
  public List<MessageAnnotation> getMessages(){
    return m_messages;
  }

  public void addMessages(Collection<MessageAnnotation> messages) {
    if (messages == null) {
      return;
    }

    // separate the general messages from specific genotype call messages
    messages.forEach(ma -> {
      if (ma.getExceptionType().equals(MessageAnnotation.TYPE_GENOTYPE)) {
        String rsid = ma.getMatches().getVariant();

        Optional<String> call = getVariantReports().stream()
            .filter(v -> v.getDbSnpId() != null && v.getDbSnpId().matches(rsid) && !v.isMissing())
            .map(VariantReport::getCall)
            .reduce((a,b) -> {throw new RuntimeException();});
        String genotype;
        if (!call.isPresent() || StringUtils.isBlank(call.get())) {
          genotype = "missing";
        }
        else {
          genotype = rsid + call.get().replaceAll("[\\|/]", "/"+rsid);
        }
        m_highlightedVariants.add(genotype);

      }
      else {
        m_messages.add(ma);
      }
    });

  }

  public List<VariantReport> getVariantReports() {
    return m_variantReports;
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
    return m_matcherDiplotypes != null && m_matcherDiplotypes.size() > 0;
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
        .map(d -> new DrugLink(d, guideline.getId(), guideline.isRxChange(), guideline.isRxPossible()))
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
    Set<String> results;
    if (!isCalled()) {
      if (sf_overrideDiplotypes.contains(getGene()) && !m_reporterDiplotypes.isEmpty()) {
        results = m_reporterDiplotypes.stream().map(Diplotype::printLookupKey).collect(Collectors.toSet());
      }
      else {
        // only HLA-B uses Other/Other, all else is Unknown/Unknown
        String defaultUncalled = getGene().equals("HLA-B") ? "Other/Other" : "Unknown/Unknown";
        results = ImmutableSet.of(m_gene + ":" + defaultUncalled);
      }
    }
    else {
      if (getGene().equals("UGT1A1") && (!isPhased() || m_reporterDiplotypes.size() > 1)) {
        results = Ugt1a1AlleleMatcher.makeLookupCalls(this).stream()
            .map(c -> "UGT1A1:"+c).collect(Collectors.toSet());
      }
      else {
        results = m_reporterDiplotypes.stream().map(Diplotype::printLookupKey).collect(Collectors.toSet());
      }
    }
    return results;
  }

  public Collection<String> printDisplayCalls() {
    if (!isCalled()) {
      if (sf_overrideDiplotypes.contains(getGene()) && !m_reporterDiplotypes.isEmpty()) {
        return m_reporterDiplotypes.stream().sorted().map(Diplotype::printDisplay).collect(Collectors.toList());
      }
      return ImmutableList.of(UNCALLED);
    }
    else if (isCallReducible()) {
      return ImmutableSet.of(
          m_matcherDiplotypes.stream()
              .map(Diplotype::printBare)
              .reduce(Diplotype.phasedReducer)
              .orElse(UNCALLED));
    }

    return m_matcherDiplotypes.stream().sorted().map(Diplotype::printDisplay).collect(Collectors.toList());
  }

  public Collection<String> printDisplayFunctions() {
    if (!isCalled()) {
      return ImmutableList.of(NA);
    }
    if (isCallReducible()) {
      return m_reporterDiplotypes.stream().sorted().map(Diplotype::printFunctionPhrase).collect(Collectors.toList());
    }
    return m_matcherDiplotypes.stream().sorted().map(Diplotype::printFunctionPhrase).collect(Collectors.toList());
  }

  public Collection<String> printDisplayPhenotypes() {
    if (!isCalled()) {
      return ImmutableList.of(NA);
    }
    return m_reporterDiplotypes.stream().map(Diplotype::getPhenotype).distinct().collect(Collectors.toSet());
  }

  /**
   * Does this gene contain an "incidental" allele that should be reported on
   */
  public boolean isIncidental() {
    return m_matcherDiplotypes.stream().anyMatch(Diplotype::isIncidental);
  }

  /**
   * Wether this gene has been marked as phased by the matcher
   */
  public boolean isPhased() {
    return m_phased;
  }

  /**
   * Can multiple phased diplotype calls be reduced down into the "+" notation. (e.g. *1/*2+*3)
   */
  private boolean isCallReducible() {
    return sf_reducibleGeneCalls.contains(getGene()) && isPhased();
  }

  public List<Diplotype> getMatcherDiplotypes() {
    return m_matcherDiplotypes;
  }

  public List<Diplotype> getReporterDiplotypes() {
    return m_reporterDiplotypes;
  }

  /**
   * Get variants that should be shown separately in the report
   */
  public List<String> getHighlightedVariants() {
    return m_highlightedVariants;
  }

  public void applyMessage(MessageVariant message) {
    if (!message.getGene().equals(getGene())) {
      return;
    }
    if (message.getRsid() != null) {
      m_variantReports.stream()
          .filter(r -> r.getDbSnpId() != null && r.getDbSnpId().equals(message.getRsid()))
          .forEach(r -> {
            r.addMessage(message);
          });
    }
  }
}
