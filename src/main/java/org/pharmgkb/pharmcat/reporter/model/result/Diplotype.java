package org.pharmgkb.pharmcat.reporter.model.result;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Objects;
import java.util.Optional;
import java.util.Set;
import java.util.TreeSet;
import java.util.stream.Collectors;
import java.util.stream.Stream;
import com.google.gson.annotations.Expose;
import com.google.gson.annotations.SerializedName;
import org.apache.commons.lang3.ObjectUtils;
import org.apache.commons.lang3.StringUtils;
import org.checkerframework.checker.nullness.qual.Nullable;
import org.pharmgkb.pharmcat.Env;
import org.pharmgkb.pharmcat.phenotype.model.DiplotypeRecord;
import org.pharmgkb.pharmcat.phenotype.model.GenePhenotype;
import org.pharmgkb.pharmcat.phenotype.model.OutsideCall;
import org.pharmgkb.pharmcat.reporter.DiplotypeFactory;
import org.pharmgkb.pharmcat.reporter.TextConstants;
import org.pharmgkb.pharmcat.reporter.model.DataSource;
import org.pharmgkb.pharmcat.reporter.model.MessageAnnotation;
import org.pharmgkb.pharmcat.reporter.model.VariantReport;
import org.pharmgkb.pharmcat.util.ActivityUtils;
import org.pharmgkb.pharmcat.util.HaplotypeNameComparator;

import static org.pharmgkb.pharmcat.reporter.caller.DpydCaller.isDpyd;


/**
 * Model class to represent a diplotype and all derived information
 *
 * @author Ryan Whaley
 */
public class Diplotype implements Comparable<Diplotype> {
  public static final String DELIMITER = "/";

  private static final String sf_toStringPattern = "%s:%s";
  private static final String sf_phenoScoreFormat = "%s (%s)";
  private static final String sf_termDelimiter = "; ";
  private static final String sf_homTemplate = "Two %s alleles";
  private static final String sf_hetTemplate = "One %s allele and one %s allele";
  private static final String sf_hetSuffix = " (heterozygous)";
  private static final String sf_homSuffix = " (homozygous)";

  @Expose
  @SerializedName("allele1")
  private Haplotype m_allele1;
  @Expose
  @SerializedName("allele2")
  private Haplotype m_allele2;
  @Expose
  @SerializedName("gene")
  private String m_gene;
  @Expose
  @SerializedName("phenotypes")
  private List<String> m_phenotypes = new ArrayList<>();
  @Expose
  @SerializedName("outsidePhenotypes")
  private boolean m_outsidePhenotypes;
  @Expose
  @SerializedName("activityScore")
  private String m_activityScore;
  @Expose
  @SerializedName("outsideActivityScore")
  private boolean m_outsideActivityScore;
  @Expose
  @SerializedName("variant")
  private VariantReport m_variant;
  @Expose
  @SerializedName("lookupKey")
  private List<String> m_lookupKeys = new ArrayList<>();
  @Expose
  @SerializedName("label")
  private String m_label;
  @Expose
  @SerializedName("observed")
  private Observation m_observed = Observation.INPUT;
  @Expose
  @SerializedName("combination")
  private boolean m_combination = false;


  /**
   * Private constructor for GSON.
   */
  private Diplotype() {
  }


  /**
   * Public constructor.
   */
  public Diplotype(String gene, Haplotype h1, @Nullable Haplotype h2, Env env, DataSource source) {
    m_gene = gene;
    m_allele1 = h1;
    m_allele2 = h2;
    m_label = printBare();
    annotateDiplotype(env, source);
  }

  public Diplotype(String gene, String hap1, @Nullable String hap2, Env env, DataSource source) {
    m_gene = gene;
    m_allele1 = env.makeHaplotype(gene, hap1, source);
    m_allele2 = hap2 == null ? null : env.makeHaplotype(gene, hap2, source);
    m_label = printBare();
    annotateDiplotype(env, source);
  }

  /**
   * Constructor for tests.
   */
  public Diplotype(String gene, String phenotype) {
    m_gene = gene;
    m_allele1 = null;
    m_allele2 = null;
    m_label = phenotype;
    m_phenotypes.add(phenotype);
    m_lookupKeys.add(phenotype);
  }

  /**
   * Constructor for a diplotype coming from an outside call that does not specify individual alleles but DOES give
   * phenotype or activity score.
   *
   * @param outsideCall an {@link OutsideCall} object
   */
  public Diplotype(GeneReport geneReport, OutsideCall outsideCall, Env env, DataSource source) {
    m_gene = outsideCall.getGene();
    if (outsideCall.getDiplotype() != null) {
      String[] alleles = DiplotypeFactory.splitDiplotype(m_gene, outsideCall.getDiplotype());
      m_allele1 = env.makeHaplotype(m_gene, alleles[0], source);
      m_allele2 = alleles.length == 2 ? env.makeHaplotype(m_gene, alleles[1], source) : null;
      m_label = printBare();

      GenePhenotype gp = env.getPhenotype(m_gene, source);

      // must set activity score before phenotype because phenotype calculation may depend on activity score
      if (outsideCall.getActivityScore() != null) {
        m_outsideActivityScore = true;
        m_activityScore = outsideCall.getActivityScore();
        // check for activity score mismatch
        String expected = lookupExpectedActivityScore(outsideCall.getPhenotype(), gp);
        if (!m_activityScore.equals(expected)) {
          if (expected == null) {
            expected = "no assigned activity score";
          }
          geneReport.addMessage(new MessageAnnotation(MessageAnnotation.TYPE_NOTE, "warn.mismatch.score",
              "Activity score from outside call (" + outsideCall.getPhenotype() + ") does not match expected " +
                  "activity score for " + m_gene + " " + m_label + " (" + expected + ")"));
        }
      }
      if (outsideCall.getPhenotype() != null) {
        m_outsidePhenotypes = true;
        m_phenotypes = List.of(outsideCall.getPhenotype());
        // check for phenotype assignment mismatch
        String expected = findExpectedPhenotype(outsideCall.getActivityScore(), gp)
            .orElse("no assigned phenotype");
        if (!outsideCall.getPhenotype().equals(expected)) {
          geneReport.addMessage(new MessageAnnotation(MessageAnnotation.TYPE_NOTE, "warn.mismatch.phenotype",
              "Phenotype from outside call  (" + outsideCall.getPhenotype() + ") does not match expected " +
                  "phenotype for " + m_gene + " " + m_label + " (" + expected + ")"));
        }
      }

      annotateDiplotype(env, source);
    }
    else {
      // no diplotypes - phenotype or activity score only

      if (outsideCall.getPhenotype() != null) {
        m_phenotypes.add(outsideCall.getPhenotype());
        m_outsidePhenotypes = true;
      }
      m_activityScore = outsideCall.getActivityScore();
      if (m_activityScore != null) {
        m_outsideActivityScore = true;
      }
      if (outsideCall.getPhenotype() != null) {
        m_label = outsideCall.getPhenotype();
      }
      else {
        m_label = outsideCall.getActivityScore();
      }
      annotateDiplotype(env, source);
    }
  }





  /**
   * Gets the gene this diplotype is for.
   * @return a HGNC gene symbol
   */
  public String getGene() {
    return m_gene;
  }

  public String getLabel() {
    return m_label;
  }

  /**
   * Gets the first {@link Haplotype} listed in this diplotype
   */
  public Haplotype getAllele1() {
    return m_allele1;
  }

  /**
   * Gets the second {@link Haplotype} listed in this diplotype
   */
  public @Nullable Haplotype getAllele2()  {
    return m_allele2;
  }

  /**
   * Is this Diplotype for a single-ploidy gene like MT-RNR1?
   * We acknowledge this isn't a "real" diplotype if this is true but for the purposes of this system we will call it a
   * diplotype.
   *
   * @return true if this Diplotype is single ploidy
   */
  private boolean isSinglePloidy() {
    return m_allele2 == null;
  }

  /**
   * Does this diplotype have an allele with the given name
   * @param alleleName an allele name, e.g. "*10"
   * @return true if this diplotype contains an allele with the given name
   */
  public boolean hasAllele(String alleleName) {
    return (m_allele1 != null && m_allele1.getName().equals(alleleName))
        || (m_allele2 != null && m_allele2.getName().equals(alleleName));
  }

  public boolean isUnknownPhenotype() {
    return m_phenotypes.size() == 0 || m_phenotypes.contains(TextConstants.NO_RESULT);
  }

  private boolean isUnknownAllele1() {
    return m_allele1 == null || m_allele1.isUnknown();
  }

  private boolean isUnknownAllele2() {
    return m_allele2 == null || m_allele2.isUnknown();
  }

  public boolean isUnknownAlleles() {
    return isUnknownAllele1() && isUnknownAllele2();
  }

  public boolean isUnknown() {
    return isUnknownPhenotype() && isUnknownAlleles();
  }

  /**
   * Gets a Sting representation of this haplotype with no gene prefix (e.g. *1/*10)
   */
  public String printBare() {
    return printOverride().orElseGet(() -> {
      if (m_allele1 == null && m_allele2 == null) {
        if (!TextConstants.isUnspecified(m_activityScore)) {
          if (isUnknownPhenotype()) {
            return m_activityScore;
          } else {
            return String.format(sf_phenoScoreFormat, String.join(DELIMITER, m_phenotypes), getActivityScore());
          }
        }
        else if (m_phenotypes.size() > 0 && !m_phenotypes.contains(TextConstants.NO_RESULT)) {
          return String.join(DELIMITER, m_phenotypes);
        }
        else return TextConstants.NA;
      }
      else {
        if (m_allele1 != null && m_allele2 != null) {
          String[] alleles = new String[]{ m_allele1.getName(), m_allele2.getName() };
          Arrays.sort(alleles, HaplotypeNameComparator.getComparator());
          return String.join(DELIMITER, alleles);
        } else if (m_allele1 != null) {
          return m_allele1.getName();
        } else {
          return TextConstants.NA;
        }
      }
    });
  }

  /**
   * Gets a String representation of this Diplotype that can be used to display in output. This should <em>NOT</em> be
   * used for matching to {@link AnnotationReport}s.
   *
   * @return a String display for this diplotype, without Gene symbol
   */
  public String printDisplay() {
    if (getVariant() != null) {
      return getVariant().printDisplay();
    }
    else {
      return printBare();
    }
  }

  /**
   * Gets a String phrase describing the individual haplotype functions, e.g. "Two low function alleles"
   * <p>
   * Will print a default N/A String if no functions exist
   */
  public String printFunctionPhrase() {

    String f1 = m_allele1 != null && m_allele1.getFunction() != null ?
        m_allele1.getFunction().toLowerCase() : null;
    String f2 = m_allele2 != null && m_allele2.getFunction() != null ?
        m_allele2.getFunction().toLowerCase() : null;

    if (!isSinglePloidy() && StringUtils.isNotBlank(f1) && StringUtils.isNotBlank(f2)) {

      if (f1.equals(f2)) {
        return String.format(sf_homTemplate, f1);
      }
      else {
        String[] functions = new String[]{f1, f2};
        Arrays.sort(functions);
        return String.format(sf_hetTemplate, functions[0], functions[1]);
      }

    } else if (isSinglePloidy() && StringUtils.isNotBlank(f1)) {
      return f1;
    }
    return TextConstants.NA;
  }

  /**
   * Gets a String representation of this diplotype with the gene as prefix, e.g. GENEX:*1/*2
   */
  public String toString() {
    return String.format(sf_toStringPattern, m_gene, printBare());
  }


  public List<String> getPhenotypes() {
    return m_phenotypes;
  }

  public String printPhenotypes() {
    if (isDpyd(m_gene)) {
      return TextConstants.SEE_DRUG;
    }
    if (m_phenotypes.size() == 0) {
      return TextConstants.NA;
    } else {
      return String.join(sf_termDelimiter, m_phenotypes);
    }
  }

  /**
   * Gets whether the phenotypes are from an outside call.
   */
  public boolean isOutsidePhenotypes() {
    return m_outsidePhenotypes;
  }


  /**
   * Print the overriding diplotype string if it exists.  Only applicable for CFTR.
   *
   * @return Optional diplotype string to override whatever the actual string would be
   */
  private Optional<String> printOverride() {
    if (m_gene.equals("CFTR")) {
      boolean refAllele1 = m_allele1.isReference();
      boolean refAllele2 = m_allele2 != null && m_allele2.isReference();
      if (refAllele1 && refAllele2 || m_allele2 == null) {
        // homozygous reference
        return Optional.of("No CPIC variants found");
      }
      if (refAllele1 || refAllele2) {
        // heterozygous
        String allele = refAllele1 ? m_allele2.getName() : m_allele1.getName();
        return Optional.of(allele + sf_hetSuffix);
      }
    }
    return Optional.empty();
  }

  /**
   * Makes a {@link Stream} of zygosity descriptors for this diplotype, e.g. *60 (heterozygous). This is a stream since
   * a single Diplotype can be described in 0, 1, or 2 Strings, depending on the particular allele.
   * @return a Stream of 0 or more zygosity Strings
   */
  public Stream<String> streamAllelesByZygosity() {
    if (m_allele1.equals(m_allele2)) {
      if (m_allele1.isReference()) {
        return Stream.empty();
      }
      return Stream.of(m_allele1.getName() + sf_homSuffix);
    }
    else {
      Set<String> alleles = new TreeSet<>(HaplotypeNameComparator.getComparator());
      if (!m_allele1.isReference()) {
        alleles.add(m_allele1.getName() + sf_hetSuffix);
      }
      if (m_allele2 != null && !m_allele2.isReference()) {
        alleles.add(m_allele2.getName() + sf_hetSuffix);
      }
      return alleles.stream();
    }
  }

  /**
   * Gets a variant used to make this diplotype call
   * @return a Variant used in this call
   */
  public VariantReport getVariant() {
    return m_variant;
  }

  public void setVariant(VariantReport variant) {
    m_variant = variant;
  }

  @Override
  public int compareTo(Diplotype o) {
    int rez = ObjectUtils.compare(getGene(), o.getGene());
    if (rez != 0) {
      return rez;
    }

    rez = ObjectUtils.compare(m_allele1, o.getAllele1());
    if (rez != 0) {
      return rez;
    }

    return ObjectUtils.compare(m_allele2, o.getAllele2());
  }

  public List<String> getLookupKeys() {
    return m_lookupKeys;
  }

  public String printLookupKeys() {
    return String.join(sf_termDelimiter, m_lookupKeys);
  }

  private static Map<String, Integer> computeLookupMap(Haplotype hap1, Haplotype hap2, String phenotype) {
    Map<String, Integer> lookupMap = new HashMap<>();

    if (hap1 == null && hap2 == null) {
      if (phenotype != null) {
        lookupMap.put(phenotype, 1);
        return lookupMap;
      }
    } else {
      if (hap1 != null) {
        if (hap2 != null) {
          if (hap1.equals(hap2)) {
            lookupMap.put(hap1.getName(), 2);
          } else {
            lookupMap.put(hap1.getName(), 1);
            lookupMap.put(hap2.getName(), 1);
          }
        } else {
          lookupMap.put(hap1.getName(), 1);
        }
      }
    }
    return lookupMap;
  }

  /**
   * Computes lookup map for this {@link Diplotype}.
   * Meant to be used internally.  This is only public for tests.
   */
  public Map<String, Integer> computeLookupMap() {
    String phenotype = null;
    if (m_phenotypes.size() == 1) {
      phenotype = m_phenotypes.get(0);
    }
    return computeLookupMap(m_allele1, m_allele2, phenotype);
  }

  /**
   * True if this diplotype does not use allele function to assign phenotype but instead relies on the presence or absense of
   * alleles for its phenotypes (e.g. HLA's)
   * @return true if this diplotype assigns phenotype based on allele presence
   */
  public boolean isAllelePresenceType() {
    return GeneReport.isAllelePresenceType(getGene());
  }

  /**
   * True if this diplotype only has the phenotype available, not the individual allele(s)
   * @return true if this diplotype only has phenotype and not individual allele(s)
   */
  public boolean isPhenotypeOnly() {
    return !isUnknownPhenotype() && isUnknownAlleles();
  }


  public String getActivityScore() {
    return m_activityScore;
  }

  public void setActivityScore(String activityScore) {
    m_activityScore = activityScore;
  }


  /**
   * Gets whether the activity score is from an outside call.
   */
  public boolean isOutsideActivityScore() {
    return m_outsideActivityScore;
  }

  public boolean hasActivityScore() {
    return !TextConstants.isUnspecified(m_activityScore);
  }

  public Observation getObserved() {
    return m_observed;
  }

  public void setObserved(Observation observed) {
    m_observed = observed;
  }

  public boolean isCombination() {
    return m_combination;
  }

  public void setCombination(boolean combination) {
    m_combination = combination;
  }


  /**
   * Fill in the phenotype and activity score depending on what information is already in the diplotype.
   */
  private void annotateDiplotype(Env env, DataSource source) {
    if (m_gene.startsWith("HLA")) {
      if (!isUnknownAlleles() && isUnknownPhenotype()) {
        m_phenotypes = DiplotypeFactory.makeHlaPhenotype(this);
        m_lookupKeys = m_phenotypes;
      }
      else if (isUnknownAlleles() && !isUnknownPhenotype()) {
        m_lookupKeys = m_phenotypes;
      }
      return;
    }

    GenePhenotype gp = env.getPhenotype(m_gene, source);
    if (gp == null) {
      return;
    }

    Map<String, Integer> lookupMap = computeLookupMap();
    if (gp.isMatchedByActivityScore()) {
      if (m_activityScore == null) {
        // diplotype -> activity score
        m_activityScore = ActivityUtils.normalize(lookupKeys("activity score", gp, lookupMap));
      }
      if (m_phenotypes.size() == 0) {
        String phenotype;
        if (!TextConstants.NA.equals(m_activityScore)) {
          phenotype = lookupPhenotypeForActivityScore(gp, m_activityScore);
          if (phenotype.equals(TextConstants.NA)) {
            phenotype = lookupPhenotypes(gp, lookupMap);
          }
        } else {
          phenotype = lookupPhenotypes(gp, lookupMap);
        }
        m_phenotypes.add(phenotype);
        // recompute lookup map because may have updated phenotypes
        lookupMap = computeLookupMap();
      }

    } else {
      // TODO(whaleyr): confirm that we don't care about activity score in this branch
      if (m_phenotypes.size() == 0) {
        // phenotype based on diplotype
        m_phenotypes = List.of(lookupPhenotypes(gp, lookupMap));
      }
    }

    if (isUnknown()) {
      m_lookupKeys = List.of(TextConstants.NO_RESULT);
    } else {
      m_lookupKeys = List.of(lookupKeys("key", gp, lookupMap));
      gp.assignActivity(m_allele1);
      gp.assignActivity(m_allele2);
    }
  }

  private String lookupKeys(String keyType, GenePhenotype gp, Map<String, Integer> lookupMap) {
    Set<String> keys = gp.getDiplotypes().stream()
        .filter(d -> d.getDiplotypeKey().equals(lookupMap))
        .map(DiplotypeRecord::getLookupKey)
        .filter(Objects::nonNull)
        .collect(Collectors.toSet());
    if (keys.size() > 1) {
      throw new IllegalStateException("More than one " + keyType + " matched for " + this + ": " +
          String.join("; ", keys));
    } else if (keys.size() == 0) {
      return TextConstants.NA;
    } else {
      return keys.iterator().next();
    }
  }


  private String lookupPhenotypes(GenePhenotype gp, Map<String, Integer> lookupMap) {
    if (isUnknownAlleles()) {
      return TextConstants.NO_RESULT;
    }
    Set<String> keys = gp.getDiplotypes().stream()
        .filter(d -> d.getDiplotypeKey().equals(lookupMap))
        .map(DiplotypeRecord::getGeneResult)
        .filter(Objects::nonNull)
        .collect(Collectors.toSet());
    if (keys.size() > 1) {
      throw new IllegalStateException("More than one phenotype matched for " + this + ": " +
          String.join("; ", keys));
    } else if (keys.size() == 0) {
      return TextConstants.NA;
    } else {
      return keys.iterator().next();
    }
  }

  private String lookupPhenotypeForActivityScore(GenePhenotype gp, String activityScore) {
    return gp.getDiplotypes().stream()
        .filter(d -> d.getLookupKey().equals(activityScore))
        .map(DiplotypeRecord::getGeneResult)
        .filter(Objects::nonNull)
        .findFirst()
        .orElse(TextConstants.NA);
  }

  /**
   * Look up activity score based on either haplotypes or phenotype.
   * Used to check outside call activity score.
   */
  private @Nullable String lookupExpectedActivityScore(@Nullable String phenotype, GenePhenotype gp) {
    if (gp != null && gp.isMatchedByActivityScore()) {
      Map<String, Integer> lookupMap = computeLookupMap(m_allele1, m_allele2, phenotype);
      return ActivityUtils.normalize(lookupKeys("activity score", gp, lookupMap));
    }
    return null;
  }


  /**
   * Look up phenotype based on either haplotypes or activity score.
   * Used to check outside call phenotype.
   */
  private Optional<String> findExpectedPhenotype(String activityScore, GenePhenotype gp) {
    if (gp != null) {
      if (gp.isMatchedByActivityScore() && activityScore != null) {
        return Optional.of(lookupPhenotypeForActivityScore(gp, activityScore));
      }
      Map<String, Integer> lookupMap = computeLookupMap(m_allele1, m_allele2, null);
      return Optional.of(lookupPhenotypes(gp, lookupMap));
    }
    return Optional.empty();
  }
}
