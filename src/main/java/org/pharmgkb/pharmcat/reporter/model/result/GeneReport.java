package org.pharmgkb.pharmcat.reporter.model.result;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.List;
import java.util.Map;
import java.util.Objects;
import java.util.Optional;
import java.util.Set;
import java.util.SortedSet;
import java.util.TreeSet;
import java.util.function.Predicate;
import java.util.stream.Collectors;
import com.google.common.base.Preconditions;
import com.google.common.collect.ImmutableList;
import com.google.common.collect.ImmutableSet;
import com.google.gson.annotations.Expose;
import com.google.gson.annotations.SerializedName;
import org.apache.commons.lang3.StringUtils;
import org.checkerframework.checker.nullness.qual.Nullable;
import org.pharmgkb.common.comparator.HaplotypeNameComparator;
import org.pharmgkb.pharmcat.Env;
import org.pharmgkb.pharmcat.definition.DefinitionReader;
import org.pharmgkb.pharmcat.definition.model.DefinitionFile;
import org.pharmgkb.pharmcat.definition.model.NamedAllele;
import org.pharmgkb.pharmcat.haplotype.NamedAlleleMatcher;
import org.pharmgkb.pharmcat.haplotype.model.BaseMatch;
import org.pharmgkb.pharmcat.haplotype.model.GeneCall;
import org.pharmgkb.pharmcat.phenotype.Phenotyper;
import org.pharmgkb.pharmcat.reporter.DiplotypeFactory;
import org.pharmgkb.pharmcat.reporter.DpydCaller;
import org.pharmgkb.pharmcat.reporter.VariantReportFactory;
import org.pharmgkb.pharmcat.reporter.caller.Slco1b1CustomCaller;
import org.pharmgkb.pharmcat.reporter.model.DataSource;
import org.pharmgkb.pharmcat.reporter.model.DrugLink;
import org.pharmgkb.pharmcat.reporter.model.MessageAnnotation;
import org.pharmgkb.pharmcat.reporter.model.OutsideCall;
import org.pharmgkb.pharmcat.reporter.model.VariantReport;
import org.pharmgkb.pharmcat.util.DataManager;

import static org.pharmgkb.pharmcat.reporter.DpydCaller.hasTrueDiplotype;
import static org.pharmgkb.pharmcat.reporter.DpydCaller.isDpyd;
import static org.pharmgkb.pharmcat.reporter.TextConstants.*;


/**
 * This class is used to help collect Gene-related data for later reporting
 */
public class GeneReport implements Comparable<GeneReport> {
  // never display these genes in the gene call list
  private static final Set<String> IGNORED_GENES       = ImmutableSet.of("IFNL4");
  /**
   * These genes have custom callers that can potentially infer a diplotype eventhough the {@link NamedAlleleMatcher}
   * cannot make a call.
   */
  public static final Set<String> OVERRIDE_DIPLOTYPES = ImmutableSet.of("SLCO1B1");
  private static final Set<String> SINGLE_PLOIDY       = ImmutableSet.of("G6PD", "MT-RNR1");
  private static final Set<String> CHROMO_X            = ImmutableSet.of("G6PD");
  private static final Set<String> ALLELE_PRESENCE     = ImmutableSet.of("HLA-A", "HLA-B");
  public static  final String YES = "Yes";
  public static  final String NO  = "No";

  private static final DefinitionReader sf_definitionReader = new DefinitionReader();


  @Expose
  @SerializedName("phenotypeVersion")
  private String m_phenotypeVersion;
  @Expose
  @SerializedName("phenotypeSource")
  private DataSource m_phenotypeSource;

  @Expose
  @SerializedName("geneSymbol")
  private final String f_gene;
  @Expose
  @SerializedName("chr")
  private String m_chr;
  @Expose
  @SerializedName("phased")
  private boolean m_phased = false;
  @Expose
  @SerializedName("effectivelyPhased")
  private boolean m_effectivelyPhased = false;
  @Expose
  @SerializedName("callSource")
  private final CallSource m_callSource;
  @Expose
  @SerializedName("uncalledHaplotypes")
  private final SortedSet<String> m_uncalledHaplotypes = new TreeSet<>(HaplotypeNameComparator.getComparator());
  @Expose
  @SerializedName("messages")
  private final List<MessageAnnotation> m_messages = new ArrayList<>();
  @Expose
  @SerializedName("relatedDrugs")
  private SortedSet<DrugLink> m_relatedDrugs = new TreeSet<>();

  @Expose
  @SerializedName("matcherDiplotypes")
  private final List<Diplotype> m_matcherDiplotypes = new ArrayList<>();
  @Expose
  @SerializedName("componentDiplotypes")
  private final SortedSet<Diplotype> m_componentDiplotypes = new TreeSet<>();

  @Expose
  @SerializedName("matcherAlleles")
  private final List<NamedAllele> f_matcherAlleles = new ArrayList<>();
  @Expose
  @SerializedName("diplotypes")
  private final List<Diplotype> m_reporterDiplotypes = new ArrayList<>();
  @Expose
  @SerializedName("variants")
  private final List<VariantReport> m_variantReports = new ArrayList<>();
  @Expose
  @SerializedName("variantsOfInterest")
  private final List<VariantReport> m_variantOfInterestReports = new ArrayList<>();
  @Expose
  @SerializedName("highlightedVariants")
  private final List<String> m_highlightedVariants = new ArrayList<>();


  /**
   * Constructor for genes that get their data from {@link GeneCall}.
   */
  public GeneReport(GeneCall call, Env env, DataSource phenotypeSource) {
    Preconditions.checkNotNull(phenotypeSource);

    f_gene = call.getGene();
    m_phenotypeSource = phenotypeSource;
    m_phenotypeVersion = env.getPhenotypeVersion(f_gene, phenotypeSource);
    m_callSource = CallSource.MATCHER;

    m_chr = call.getChromosome();
    m_uncalledHaplotypes.clear();
    m_uncalledHaplotypes.addAll(call.getUncallableHaplotypes());
    m_phased = call.isPhased();
    m_effectivelyPhased = call.isEffectivelyPhased();
    try {
      VariantReportFactory variantReportFactory = new VariantReportFactory(f_gene, m_chr);
      call.getVariants().stream()
          .map(variantReportFactory::make).forEach(m_variantReports::add);
      call.getMatchData().getMissingPositions().stream()
          .map(variantReportFactory::make).forEach(m_variantReports::add);
      call.getVariantsOfInterest().stream()
          .map(variantReportFactory::make).forEach(m_variantOfInterestReports::add);
    } catch (IOException ex) {
      throw new RuntimeException("Could not start a gene report", ex);
    }

    // set the flag in reports for the variants with mismatched alleles
    call.getMatchData().getMismatchedPositions().stream()
        .flatMap(a -> m_variantReports.stream().filter(v -> v.getPosition() == a.getPosition()))
        .forEach(r -> r.setMismatch(true));

    call.getHaplotypes().stream()
        .map(BaseMatch::getHaplotype)
        .filter(Predicate.not(Objects::isNull))
        .forEach(f_matcherAlleles::add);

    DiplotypeFactory diplotypeFactory = new DiplotypeFactory(f_gene, env);
    //noinspection rawtypes
    Collection matches = null;
    if (isDpyd(f_gene)) {
      if (hasTrueDiplotype(call)) {
        m_componentDiplotypes.addAll(diplotypeFactory.makeComponentDiplotypes(call, m_phenotypeSource));
      } else {
        matches = call.getHaplotypes();
      }
    }
    if (matches == null) {
      matches = call.getDiplotypes();
    }
    m_matcherDiplotypes.addAll(diplotypeFactory.makeDiplotypes(matches, m_phenotypeSource));

    if (isDpyd(f_gene)) {
      m_reporterDiplotypes.addAll(DpydCaller.inferDiplotypes(matches, hasTrueDiplotype(call), env, phenotypeSource));
    } else if (Slco1b1CustomCaller.shouldBeUsedOn(this)) {
      Slco1b1CustomCaller
          .makeLookupCalls(this, diplotypeFactory, m_phenotypeSource)
          .ifPresent(m_reporterDiplotypes::add);
    } else {
      m_reporterDiplotypes.addAll(diplotypeFactory.makeDiplotypes(matches, m_phenotypeSource));
    }

    applyMatcherMessages(call);
  }


  /**
   * Constructor for genes that get their data from an {@link OutsideCall} that comes from the {@link Phenotyper}.
   */
  public GeneReport(OutsideCall call, Env env, DataSource phenotypeSource) {
    Preconditions.checkNotNull(phenotypeSource);

    f_gene = call.getGene();
    m_phenotypeSource = phenotypeSource;
    m_phenotypeVersion = env.getPhenotypeVersion(f_gene, phenotypeSource);
    m_callSource = CallSource.OUTSIDE;

    addOutsideCall(call, env);
    applyOutsideCallMessages();
  }

  /**
   * Enable adding multiple {@link OutsideCall}s for same gene as multiple diplotypes.
   */
  public void addOutsideCall(OutsideCall call, Env env) {
    Preconditions.checkState(m_callSource == CallSource.OUTSIDE);
    DiplotypeFactory diplotypeFactory = new DiplotypeFactory(f_gene, env);

    if (isDpyd(f_gene)) {
      m_reporterDiplotypes.addAll(DpydCaller.inferDiplotypes(call.getDiplotypes(), true, env,
          m_phenotypeSource));
    } else {
      m_reporterDiplotypes.addAll(diplotypeFactory.makeDiplotypes(call, m_phenotypeSource));
    }
  }


  /**
   * Constructor for unspecified {@link GeneReport}.
   */
  private GeneReport(String geneSymbol, Env env, DataSource phenotypeSource) {
    f_gene = geneSymbol;
    m_phenotypeSource = phenotypeSource;
    m_phenotypeVersion = env.getPhenotypeVersion(f_gene, phenotypeSource);
    m_callSource = CallSource.NONE;
    Diplotype unknownDiplotype = DiplotypeFactory.makeUnknownDiplotype(geneSymbol, env, phenotypeSource);
    m_matcherDiplotypes.add(unknownDiplotype);
    m_reporterDiplotypes.add(unknownDiplotype);
  }

  public static GeneReport unspecifiedGeneReport(String gene, Env env, DataSource source) {
    return new GeneReport(gene, env, source);
  }


  /**
   * Constructor for tests.
   */
  protected GeneReport(String geneSymbol, DataSource phenotypeSource, @Nullable String phenotypeVersion) {
    Preconditions.checkNotNull(phenotypeSource);

    f_gene = geneSymbol;
    m_phenotypeSource = phenotypeSource;
    m_phenotypeVersion = phenotypeVersion;
    m_callSource = CallSource.NONE;
  }

  /**
   * Should only be used by tests!
   */
  protected void addReporterDiplotype(Diplotype diplotype) {
    m_reporterDiplotypes.add(diplotype);
  }


  /**
   * The version this gene report is based on.
   */
  public @Nullable String getPhenotypeVersion() {
    return m_phenotypeVersion;
  }

  public DataSource getPhenotypeSource() {
    return m_phenotypeSource;
  }


  /**
   * The gene symbol for this gene.
   */
  public String getGene(){
    return f_gene;
  }

  /**
   * The gene symbol for display purposes, like in a report. This accounts for the IFNL3/4 problem
   */
  public String getGeneDisplay() {
    return switch (f_gene) {
      case "IFNL3", "IFNL4" -> "IFNL3/4";
      default -> f_gene;
    };
  }

  /**
   * The chromosome this gene appears on
   */
  public String getChr() {
    return m_chr;
  }

  /**
   * The list of messages that apply to this gene
   */
  public List<MessageAnnotation> getMessages(){
    return m_messages;
  }

  /**
   * Add a message annotation. Separates the general messages from specific genotype call messages
   */
  public void addMessage(MessageAnnotation ma) {
    if (ma == null) return;
    if (ma.getExceptionType().equals(MessageAnnotation.TYPE_GENOTYPE)) {
      String rsid = ma.getMatches().getVariant();

      Optional<String> call = getVariantReports().stream()
          .filter(v -> v.getDbSnpId() != null && v.getDbSnpId().matches(rsid) && !v.isMissing())
          .map(VariantReport::getCall)
          .reduce((a,b) -> {throw new RuntimeException();});
      String genotype;
      if (call.isEmpty() || StringUtils.isBlank(call.get())) {
        genotype = "missing";
      }
      else {
        genotype = rsid + call.get().replaceAll("[|/]", "/"+rsid);
      }
      m_highlightedVariants.add(genotype);

    }
    else {
      m_messages.add(ma);
    }
  }

  public void addMessages(Collection<MessageAnnotation> messages) {
    if (messages == null) {
      return;
    }
    messages.forEach(this::addMessage);
  }

  /**
   * This will add the applicable variant warnings to the {@link VariantReport} objects in this {@link GeneReport}
   * @param variantWarnings a Map of all variant warnings by "chr:position" strings
   */
  public void addVariantWarningMessages(Map<String, Collection<String>> variantWarnings) {
    if (variantWarnings != null && variantWarnings.size() > 0) {
      for (VariantReport variantReport : m_variantReports) {
        Collection<String> warnings = variantWarnings.get(variantReport.toChrPosition());
        if (warnings != null && warnings.size() > 0) {
          variantReport.setWarnings(warnings);
        }
      }
    }
  }

  public List<VariantReport> getVariantReports() {
    return m_variantReports;
  }

  /**
   * Finds the first {@link VariantReport} in this gene report that has the specified RSID. Since RSIDs should be unique
   * to variants this should never silently hide multiple matches.
   * @param rsid a dbSNP RSID that identifies a variant
   * @return an {@link Optional} {@link VariantReport}
   */
  public Optional<VariantReport> findVariantReport(String rsid) {
    return m_variantReports.stream().filter(v -> v.getDbSnpId() != null && v.getDbSnpId().contains(rsid)).findFirst();
  }

  public List<VariantReport> getVariantOfInterestReports() {
    return m_variantOfInterestReports;
  }

  /**
   * Gets the Set of Haplotypes the matcher could not evaluate
   */
  public Set<String> getUncalledHaplotypes() {
    return m_uncalledHaplotypes;
  }

  public String toString() {
    return f_gene + " Report";
  }

  /**
   * True if the {@link NamedAlleleMatcher} has returned at least one call for this gene, false otherwise
   */
  public boolean isCalled() {
    return m_matcherDiplotypes.size() > 0 && m_matcherDiplotypes.stream().noneMatch(Diplotype::isUnknown);
  }

  /**
   * True if there is a diplotype for this gene that the reporter can use, false otherwise. The reporter may be able to 
   * use diplotype calls not made by the matcher (e.g. SLCO1B1)
   */
  public boolean isReportable() {
    return m_reporterDiplotypes.size() > 0 && m_reporterDiplotypes.stream().noneMatch(Diplotype::isUnknown);
  }

  /**
   * True if this gene does not use allele function to assign phenotype but instead relies on the presence or absense of
   * alleles for its phenotypes (e.g. HLA's)
   * @return true if this gene assigns phenotype based on allele presence
   */
  public boolean isAllelePresenceType() {
    return isAllelePresenceType(f_gene);
  }

  /**
   * True if the gene does not use allele function to assign phenotype but instead relies on the presence or absense of
   * alleles for its phenotypes (e.g. HLA's)
   * @param gene the gene symbol
   * @return true if this gene assigns phenotype based on allele presence
   */
  public static boolean isAllelePresenceType(String gene) {
    return ALLELE_PRESENCE.contains(gene);
  }

  @Override
  public int compareTo(GeneReport o) {
    if (this == o) {
      return 0;
    }
    int rez = Objects.compare(getGene(), o.getGene(), String.CASE_INSENSITIVE_ORDER);
    if (rez != 0) {
      return rez;
    }
    rez = m_phenotypeSource.compareTo(o.getPhenotypeSource());
    if (rez != 0) {
      return rez;
    }
    rez = CallSource.compare(m_callSource, o.getCallSource());

    return rez;
  }

  /**
   * Gets a list of {@link DrugLink} objects that are in the same guidelines as this gene
   * @return a list of {@link DrugLink} objects
   */
  public SortedSet<DrugLink> getRelatedDrugs() {
    return m_relatedDrugs;
  }

  private void addRelatedDrug(DrugLink drug) {
    if (m_relatedDrugs == null) {
      m_relatedDrugs = new TreeSet<>();
    }
    m_relatedDrugs.add(drug);
  }

  /**
   * Adds the drugs in the given <code>drugReport</code> to this report as {@link DrugLink} objects
   * @param drugReport a {@link DrugLink} with relatedDrugs
   */
  public void addRelatedDrugs(DrugReport drugReport) {

    drugReport.getRelatedDrugs().stream()
        .map(d -> new DrugLink(d, drugReport.getId()))
        .forEach(this::addRelatedDrug);
  }

  /**
   * True if the genotyping data in the report comes from outside PharmCAT, false if match is made by PharmCAT
   */
  public boolean isOutsideCall() {
    return m_callSource == CallSource.OUTSIDE;
  }

  /**
   * Gets an enum value {@link CallSource} that describes where the diplotype call for this gene came from.
   * @return an enum value of the source of the diplotype call
   */
  public CallSource getCallSource() {
    return m_callSource;
  }


  /**
   * Gets the String representations of the genotypes in this report for display.
   * This should <em>NOT</em> be used for matching!
   *
   * @return a Collection of diplotype Strings (e.g. *2/*3, *4 (heterozygote))
   */
  public List<String> printDisplayCalls() {
    if (m_callSource == CallSource.NONE) {
      return LIST_UNCALLED_NO_DATA;

    } else if (m_callSource == CallSource.MATCHER) {
      if (!isCalled()) {
        if (m_variantReports.size() == 0 ||
            m_variantReports.stream().allMatch(VariantReport::isMissing)) {
          return LIST_UNCALLED_NO_DATA;
        }
        if (OVERRIDE_DIPLOTYPES.contains(getGeneDisplay()) && isReportable()) {
          return ImmutableList.of(m_reporterDiplotypes.get(0).printDisplay());
        }
        return LIST_UNCALLED;
      }
      return m_matcherDiplotypes.stream()
          .map(Diplotype::printDisplay)
          .toList();
    } else {
      return m_reporterDiplotypes.stream()
          .map(Diplotype::printDisplay)
          .toList();
    }
  }

  public List<String> printDisplayInferredCalls() {
    if (m_callSource == CallSource.NONE) {
      return LIST_UNCALLED_NO_DATA;
    } else if (m_callSource == CallSource.MATCHER) {
      if (!isCalled()) {
        if (m_variantReports.size() == 0 ||
            m_variantReports.stream().allMatch(VariantReport::isMissing)) {
          return LIST_UNCALLED_NO_DATA;
        }
      }
    }
    return m_reporterDiplotypes.stream()
        .map(Diplotype::printDisplay)
        .toList();
  }


  /**
   * Used in the final report template in the Genotype Summary table in the "Phenotype" column
   * @return a Collection of phenotype Strings (e.g. Poor Metabolizer)
   */
  public Collection<String> printDisplayPhenotypes() {
    if (!isReportable()) {
      return ImmutableList.of(NA);
    }
    return m_reporterDiplotypes.stream().sorted().map(Diplotype::printPhenotypes).collect(Collectors.toList());
  }

  /**
   * Whether this gene has been marked as phased by the matcher.
   */
  public boolean isPhased() {
    return m_phased;
  }

  /**
   * Whether this gene has been marked as effectively phased by the matcher.
   */
  public boolean isEffectivelyPhased() {
    return m_effectivelyPhased;
  }

  /**
   * Is this gene ignored, meaning it shouldn't be included in reports
   * @return true if this gene should be left out of output reports
   */
  public boolean isIgnored() {
    return isIgnored(getGene());
  }

  /**
   * Is the specified gene symbol ignored, e.g. it shouldn't be included in reports
   * @param gene a gene symbol
   * @return true if this gene should be left out of output reports
   */
  public static boolean isIgnored(String gene) {
    return StringUtils.isNotBlank(gene) && IGNORED_GENES.contains(gene);
  }

  /**
   * Is the specified gene single ploidy, i.e. it is on chrY or chrM which occur on a single chromosome, not a pair.
   * @param gene a gene symbol
   * @return true if the gene occurs on a single chromosome, not a pair
   */
  public static boolean isSinglePloidy(String gene) {
    return StringUtils.isNotBlank(gene) && SINGLE_PLOIDY.contains(gene);
  }

  public static boolean isXChromo(String gene) {
    return StringUtils.isNotBlank(gene) && CHROMO_X.contains(gene);
  }

  /**
   * This will test whether this has a called haplotype at the reporter level, so this will not be just matched
   * haplotypes but also outside call haplotypes and inferred haplotypes
   * @param haplotype a haplotype name to test for
   * @return true if the haplotype has been called at least once (het or hom) for this gene
   */
  public boolean hasHaplotype(String haplotype) {
    return m_reporterDiplotypes.stream()
        .anyMatch((d) -> d.hasAllele(haplotype));
  }

  /**
   * Gets the list of {@link Diplotype} objects that the {@link NamedAlleleMatcher} found from data in the input VCF
   * file. This collection may not contian data if no input VCF data was used to make diplotype calls.
   * @return the list of {@link Diplotype} objects from the {@link NamedAlleleMatcher}.
   */
  public List<Diplotype> getMatcherDiplotypes() {
    return m_matcherDiplotypes;
  }

  /**
   * Gets the list of component haplotypes as {@link Diplotype}s if using lowest function mode and
   * {@link #getMatcherDiplotypes()} has an actual diplotype.
   */
  public SortedSet<Diplotype> getComponentDiplotypes() {
    return m_componentDiplotypes;
  }

  /**
   * Gets the list of {@link NamedAllele} objects that the {@link NamedAlleleMatcher} found from data in the input VCF
   * file. This collection may not contian data if no input VCF data was used to make diplotype calls.
   * @return the list of {@link NamedAllele} objects from the {@link NamedAlleleMatcher}.
   */
  public List<NamedAllele> getMatcherAlleles() {
    return f_matcherAlleles;
  }

  /**
   * Gets the list of {@link Diplotype} objects that the {@link Phenotyper} compiled from all sources, including the
   * {@link NamedAlleleMatcher}.
   * <p>
   * This list will include matcher diplotypes and diplotypes that were determined in the {@link Phenotyper}.
   * <p>
   * This is the list of diplotypes that should be used for final reporting to the user.
   * @return the list of {@link Diplotype} objects for any final reporting
   */
  public List<Diplotype> getReporterDiplotypes() {
    return m_reporterDiplotypes;
  }

  /**
   * Get variants that should be shown separately in the report
   */
  public List<String> getHighlightedVariants() {
    return m_highlightedVariants;
  }

  /**
   * Is there any missing variant in this gene?
   * @return true if there is any missing variant
   */
  public boolean isMissingVariants() {
    if (getCallSource() == CallSource.OUTSIDE) {
      return false;
    }
    return m_variantReports.stream().anyMatch(VariantReport::isMissing);
  }

  private void applyMatcherMessages(GeneCall geneCall) {
    boolean comboOrPartialCall = geneCall.getHaplotypes().stream()
        .anyMatch((h) -> h.getHaplotype() != null && (h.getHaplotype().isCombination() || h.getHaplotype().isPartial()));
    if (comboOrPartialCall && !isDpyd(geneCall.getGene())) {
      addMessage(MessageAnnotation.loadMessage(MessageAnnotation.MSG_COMBO_NAMING));

      if (!isPhased()) {
        addMessage(MessageAnnotation.loadMessage(MessageAnnotation.MSG_COMBO_UNPHASED));
      }
    }

    if (geneCall.getGene().equals("CYP2D6")) {
      addMessage(MessageAnnotation.loadMessage(MessageAnnotation.MSG_CYP2D6_GENERAL));
      addMessage(MessageAnnotation.loadMessage(MessageAnnotation.MSG_CYP2D6_MODE));
    }

    if (!geneCall.getGene().equals("CFTR")) {
      // add reference allele message
      m_matcherDiplotypes.stream()
          .map((d) -> {
            if (showReferenceMessage(d.getAllele1())) {
              return d.getAllele1();
            }
            if (d.getAllele2() != null && showReferenceMessage(d.getAllele2())) {
              return d.getAllele2();
            }
            return null;
          })
          .filter(Objects::nonNull)
          .findFirst()
          .ifPresent((h) -> {

            StringBuilder builder = new StringBuilder()
                .append("The ")
                .append(h.getGene())
                .append(" ")
                .append(h.getName())
                .append(" allele assignment is characterized by the absence of variants at the positions that are " +
                    "included in the underlying allele definitions");
            if (this.isMissingVariants()) {
              builder.append(" either because the position is reference or missing");
            }
            builder.append(".");

            addMessage(new MessageAnnotation(MessageAnnotation.TYPE_NOTE, "reference-allele", builder.toString()));
          });
    }
  }

  private void applyOutsideCallMessages() {
    addMessage(new MessageAnnotation(MessageAnnotation.TYPE_NOTE, "outside-call",
        "The call for " + getGeneDisplay() + " comes from an outside data source which does not supply " +
            "position-level detail.  For specific disclaimers and limitations, see the original genotyping source."));
  }

  /**
   * Checks if {@link Haplotype} is reference for display purposes (don't call reference if gene only has 1 varint).
   */
  private boolean showReferenceMessage(Haplotype hap) {
    if (!hap.isReference()) {
      return false;
    }
    try {
      if (sf_definitionReader.getGenes().size() == 0) {
        sf_definitionReader.read(DataManager.DEFAULT_DEFINITION_DIR);
      }
      Optional<DefinitionFile> definitionFile = sf_definitionReader.lookupDefinitionFile(hap.getGene());
      return definitionFile.isPresent() && definitionFile.get().getVariants().length > 1;

    } catch (IOException ex) {
      throw new RuntimeException("Error checking reference", ex);
    }
  }
}
