package org.pharmgkb.pharmcat.reporter.model.result;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Objects;
import java.util.Optional;
import java.util.Set;
import java.util.SortedSet;
import java.util.TreeSet;
import java.util.stream.Collectors;
import com.google.gson.annotations.Expose;
import com.google.gson.annotations.SerializedName;
import org.apache.commons.lang3.ObjectUtils;
import org.checkerframework.checker.nullness.qual.Nullable;
import org.pharmgkb.pharmcat.Env;
import org.pharmgkb.pharmcat.phenotype.model.DiplotypeRecord;
import org.pharmgkb.pharmcat.phenotype.model.GenePhenotype;
import org.pharmgkb.pharmcat.phenotype.model.OutsideCall;
import org.pharmgkb.pharmcat.reporter.DiplotypeFactory;
import org.pharmgkb.pharmcat.reporter.TextConstants;
import org.pharmgkb.pharmcat.reporter.model.DataSource;
import org.pharmgkb.pharmcat.reporter.model.VariantReport;
import org.pharmgkb.pharmcat.util.ActivityUtils;
import org.pharmgkb.pharmcat.util.HaplotypeNameComparator;


/**
 * Model class to represent a diplotype and all derived information
 *
 * @author Ryan Whaley
 */
public class Diplotype implements Comparable<Diplotype> {
  private static final String HETEROZYGOUS_SUFFIX = " (heterozygous)";
  private static final List<String> USE_CPIC_STYLE_DIPLOTYPE_NAMES = List.of(
      "CACNA1S",
      "CFTR",
      "RYR1"
  );
  private static final String sf_phenoScoreFormat = "%s (%s)";

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
  @SerializedName("outsidePhenotype")
  private boolean m_outsidePhenotype;
  @Expose
  @SerializedName("outsidePhenotypeMismatch")
  private String m_outsidePhenotypeMismatch;
  @Expose
  @SerializedName("activityScore")
  private String m_activityScore;
  @Expose
  @SerializedName("outsideActivityScore")
  private boolean m_outsideActivityScore;
  @Expose
  @SerializedName("outsideActivityScoreMismatch")
  private String m_outsideActivityScoreMismatch;
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
  @SerializedName("inferred")
  private boolean m_inferred = false;
  private SortedSet<Diplotype> m_inferredSourceDiplotypes;
  @Expose
  @SerializedName("combination")
  private boolean m_combination = false;
  @Expose
  @SerializedName("phenotypeDataSource")
  private DataSource m_phenotypeDataSource;


  /**
   * Private constructor for GSON.
   */
  @SuppressWarnings("unused")
  private Diplotype() {
  }

  /**
   * Public constructor.
   */
  public Diplotype(String gene, Haplotype h1, @Nullable Haplotype h2, Env env, DataSource source) {
    m_gene = gene;
    m_allele1 = h1;
    m_allele2 = h2;
    m_phenotypeDataSource = source;
    annotateDiplotype(env.getPhenotype(m_gene, source));
    m_label = buildLabel();
  }

  public Diplotype(String gene, String hap1, @Nullable String hap2, Env env, DataSource source) {
    m_gene = gene;
    m_allele1 = env.makeHaplotype(gene, hap1, source);
    m_allele2 = hap2 == null ? null : env.makeHaplotype(gene, hap2, source);
    m_phenotypeDataSource = source;
    annotateDiplotype(env.getPhenotype(m_gene, source));
    m_label = buildLabel();
  }

  /**
   * Constructor for tests.
   */
  public Diplotype(String gene, String phenotype) {
    m_gene = gene;
    m_allele1 = null;
    m_allele2 = null;
    m_phenotypes.add(phenotype);
    m_lookupKeys.add(phenotype);
    m_label = buildLabel();
  }

  /**
   * Constructor for a diplotype coming from an outside call that does not specify individual alleles but DOES give
   * phenotype or activity score.
   *
   * @param outsideCall an {@link OutsideCall} object
   */
  public Diplotype(OutsideCall outsideCall, Env env, DataSource source) {
    m_gene = outsideCall.getGene();
    m_phenotypeDataSource = source;
    GenePhenotype gp = env.getPhenotype(m_gene, source);

    if (outsideCall.getDiplotype() != null) {
      String[] alleles = DiplotypeFactory.splitDiplotype(m_gene, outsideCall.getDiplotype());
      m_allele1 = env.makeHaplotype(m_gene, alleles[0], source);
      m_allele2 = alleles.length == 2 ? env.makeHaplotype(m_gene, alleles[1], source) : null;
    }

    // must set activity score before phenotype because phenotype calculation may depend on activity score
    if (outsideCall.getActivityScore() != null) {
      m_outsideActivityScore = true;
      m_activityScore = outsideCall.getActivityScore();
      if (gp != null && gp.isMatchedByActivityScore()) {
        // check for activity score mismatch
        checkActivityScore(outsideCall.getPhenotype(), gp);
      }
    }
    if (outsideCall.getPhenotype() != null) {
      m_outsidePhenotype = true;
      m_phenotypes = List.of(outsideCall.getPhenotype());
      if (gp != null) {
        // check for phenotype assignment mismatch
        checkPhenotype(gp);
      }
    }
    annotateDiplotype(gp);
    m_label = buildLabel();
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
  public @Nullable Haplotype getAllele1() {
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

  private boolean isUnknownAllele1() {
    return m_allele1 == null || m_allele1.isUnknown();
  }

  private boolean isUnknownAllele2() {
    return m_allele2 == null || m_allele2.isUnknown();
  }

  public boolean isUnknownAlleles() {
    return isUnknownAllele1() && isUnknownAllele2();
  }


  public List<String> getPhenotypes() {
    return m_phenotypes;
  }

  /**
   * Gets whether the phenotypes are from an outside call.
   */
  public boolean isOutsidePhenotype() {
    return m_outsidePhenotype;
  }

  /**
   * Gets expected phenotype if outside phenotype (based on diplotype) does not match.
   */
  public String getOutsidePhenotypeMismatch() {
    return m_outsidePhenotypeMismatch;
  }

  public boolean isUnknownPhenotype() {
    return m_phenotypes.size() == 0 || m_phenotypes.contains(TextConstants.NO_RESULT);
  }


  public boolean isUnknown() {
    return isUnknownPhenotype() && isUnknownAlleles();
  }


  /**
   * Gets a variant used to make this diplotype call.
   *
   * @return a Variant used in this call
   */
  public VariantReport getVariant() {
    return m_variant;
  }

  public void setVariant(VariantReport variant) {
    m_variant = variant;
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

  public boolean hasActivityScore() {
    return !TextConstants.isUnspecified(m_activityScore);
  }

  /**
   * Gets whether the activity score is from an outside call.
   */
  public boolean isOutsideActivityScore() {
    return m_outsideActivityScore;
  }

  /**
   * Gets expected activity score if outside activity score (based on diplotype) does not match.
   */
  public String getOutsideActivityScoreMismatch() {
    return m_outsideActivityScoreMismatch;
  }


  public boolean isInferred() {
    return m_inferred;
  }

  public void setInferred(boolean inferred) {
    m_inferred = inferred;
  }


  public SortedSet<Diplotype> getInferredSourceDiplotypes() {
    return m_inferredSourceDiplotypes;
  }

  public void setInferredSourceDiplotypes(SortedSet<Diplotype> diplotypes) {
    m_inferredSourceDiplotypes = diplotypes;
  }

  public void setInferredSourceDiplotype(Diplotype diplotype) {
    m_inferredSourceDiplotypes = new TreeSet<>();
    m_inferredSourceDiplotypes.add(diplotype);
  }



  public boolean isCombination() {
    return m_combination;
  }

  public void setCombination(boolean combination) {
    m_combination = combination;
  }


  public List<String> getLookupKeys() {
    return m_lookupKeys;
  }


  /**
   * True if this diplotype does not use allele function to assign phenotype but instead relies on the presence or
   * absense of alleles for its phenotypes (e.g. HLA's).
   *
   * @return true if this diplotype assigns phenotype based on allele presence
   */
  public boolean isAllelePresenceType() {
    return GeneReport.isAllelePresenceType(getGene());
  }


  /**
   * Builds a String representation of this haplotype with no gene prefix (e.g. *1/*10) that can be used for lookups.
   */
  private String buildLabel() {

    // m_allele1 can be null if coming from outside calls
    if (USE_CPIC_STYLE_DIPLOTYPE_NAMES.contains(m_gene) && m_allele1 != null) {
      boolean isAllele1Ref = m_allele1.isReference();
      boolean isAllele2Ref = m_allele2 != null && m_allele2.isReference();
      if (isAllele1Ref && isAllele2Ref) {
         if (m_allele2 == null) {
           return TextConstants.REFERENCE;
         } else {
           // homozygous reference
           return TextConstants.HOMOZYGOUS_REFERENCE;
         }
      }
      if (isAllele1Ref || isAllele2Ref) {
        // heterozygous
        String allele = isAllele1Ref ? m_allele2.getName() : m_allele1.getName();
        return allele + HETEROZYGOUS_SUFFIX;
      }
      // fall through, use normal diplotype name
    }

    if (m_allele1 == null && m_allele2 == null) {
      if (!TextConstants.isUnspecified(m_activityScore)) {
        if (isUnknownPhenotype()) {
          return m_activityScore;
        } else {
          return String.format(sf_phenoScoreFormat,
              String.join(TextConstants.GENOTYPE_DELIMITER, m_phenotypes), getActivityScore());
        }
      }
      else if (!m_phenotypes.isEmpty() && !m_phenotypes.contains(TextConstants.NO_RESULT)) {
        return String.join(TextConstants.GENOTYPE_DELIMITER, m_phenotypes);
      }
      else return TextConstants.NA;
    }
    else {
      if (m_allele1 != null && m_allele2 != null) {
        String[] alleles = new String[]{ m_allele1.getName(), m_allele2.getName() };
        Arrays.sort(alleles, HaplotypeNameComparator.getComparator());
        return String.join(TextConstants.GENOTYPE_DELIMITER, alleles);
      } else if (m_allele1 != null) {
        return m_allele1.getName();
      } else {
        return TextConstants.NA;
      }
    }
  }


  /**
   * Gets a String representation of this diplotype with the gene as prefix, e.g. GENEX:*1/*2
   */
  public String toString() {
    return m_gene + ":" + buildLabel();
  }


  @Override
  public int compareTo(Diplotype o) {
    int rez = ObjectUtils.compare(getGene(), o.getGene());
    if (rez != 0) {
      return rez;
    }

    rez = compareAllele(m_allele1, o.getAllele1());
    if (rez != 0) {
      return rez;
    }

    rez = compareAllele(m_allele2, o.getAllele2());
    if (rez != 0) {
      return rez;
    }

    return ObjectUtils.compare(m_label, o.getLabel());
  }

  private int compareAllele(@Nullable Haplotype a, @Nullable Haplotype b) {

    if (a == b) {
      return 0;
    } else if (a == null) {
      return -1;
    } else if (b == null) {
      return 1;
    }
    return HaplotypeNameComparator.getComparator().compare(a.getName(), b.getName());
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
   * <p>
   * Using default scope to allow for use in tests.
   */
  Map<String, Integer> computeLookupMap() {
    String phenotype = null;
    if (m_phenotypes.size() == 1) {
      phenotype = m_phenotypes.get(0);
    }
    return computeLookupMap(m_allele1, m_allele2, phenotype);
  }

  private Map<String,Integer> computeDiplotypeKey() {
    return computeLookupMap(m_allele1, m_allele2, null);
  }


  /**
   * Fill in the phenotype, activity score, and lookup key depending on what information is already in the diplotype.
   */
  private void annotateDiplotype(@Nullable GenePhenotype gp) {
    if (m_gene.startsWith("HLA")) {
      if (!isUnknownAlleles() && isUnknownPhenotype()) {
        m_phenotypes = DiplotypeFactory.makeHlaPhenotype(this);
        m_lookupKeys = m_phenotypes;
      }
      else if (isUnknownAlleles()) {
        if (isUnknownPhenotype()) {
          m_lookupKeys.add(TextConstants.NO_RESULT);
        } else {
          m_lookupKeys = m_phenotypes;
        }
      }
      return;
    }

    if (gp == null) {
      return;
    }

    Map<String, Integer> lookupMap = computeDiplotypeKey();
    Optional<DiplotypeRecord> diplotype = gp.findDiplotype(lookupMap);

    if (
        gp.isMatchedByActivityScore()
            // TODO(whaleyr): temporary fix until allele activity score / function can be fixed for DPYD
            && !(getGene().equals("DPYD") && m_phenotypeDataSource == DataSource.DPWG)
    ) {
      if (m_activityScore == null) {
        if (m_phenotypes.isEmpty()) {
          // diplotype -> activity score
          diplotype.ifPresentOrElse(
              (d) -> m_activityScore = d.getActivityScore(),
              () -> m_activityScore = ActivityUtils.normalize(lookupActivityByDiplotype(gp, lookupMap)));
        } else {
          // leave activity score as NA but generate 1 lookupKey per activity score that matches phenotype below
          // reporter will handle displaying correct activity score
          m_activityScore = TextConstants.NA;
        }
      }
      if (m_phenotypes.isEmpty()) {
        if (m_outsideActivityScore) {
          // activity score -> phenotypes
          m_phenotypes.addAll(lookupPhenotypeByActivityScore(gp, m_activityScore));
        } else {
          // diplotype -> phenotype
          diplotype.ifPresentOrElse(
              (d) -> m_phenotypes.add(d.getPhenotype()),
              () -> m_phenotypes.add(lookupPhenotypesByDiplotype(gp, lookupMap)));
        }
      }
      if (!TextConstants.NA.equals(m_activityScore)) {
        m_lookupKeys.add(m_activityScore);
      } else {
        if (isUnknown()) {
          m_lookupKeys.add(TextConstants.NO_RESULT);
        } else {
          // if had diplotype, activity score would have been filled in
          m_lookupKeys.addAll(lookupActivityScoresByPhenotype(gp, m_phenotypes));
        }
      }
    }
    // non-activity score genes
    else {
      if (m_phenotypes.isEmpty()) {
        // phenotype based on diplotype
        m_phenotypes = List.of(lookupPhenotypesByDiplotype(gp, lookupMap));
      }
      if (isUnknown()) {
        m_lookupKeys.add(TextConstants.NO_RESULT);
      } else {
        m_lookupKeys = m_phenotypes;
      }
    }

    gp.assignActivity(m_allele1);
    gp.assignActivity(m_allele2);
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
    } else if (keys.isEmpty()) {
      return TextConstants.NA;
    } else {
      return keys.iterator().next();
    }
  }


  private String lookupPhenotypesByDiplotype(GenePhenotype gp, Map<String, Integer> lookupMap) {
    if (isUnknownAlleles()) {
      return TextConstants.NO_RESULT;
    }
    SortedSet<String> keys = gp.getDiplotypes().stream()
        .filter(d -> d.getDiplotypeKey().equals(lookupMap))
        .map(DiplotypeRecord::getGeneResult)
        .filter(Objects::nonNull)
        .collect(Collectors.toCollection(TreeSet::new));
    if (keys.size() > 1) {
      throw new IllegalStateException("More than one phenotype matched for " + this + ": " +
          String.join("; ", keys));
    } else if (keys.isEmpty()) {
      return TextConstants.NA;
    } else {
      return keys.iterator().next();
    }
  }

  private String lookupActivityByDiplotype(GenePhenotype gp, Map<String, Integer> lookupMap) {
    if (isUnknownAlleles()) {
      return TextConstants.NO_RESULT;
    }
    SortedSet<String> keys = gp.getDiplotypes().stream()
        .filter(d -> d.getDiplotypeKey().equals(lookupMap))
        .map(DiplotypeRecord::getActivityScore)
        .filter(Objects::nonNull)
        .collect(Collectors.toCollection(TreeSet::new));
    if (keys.size() > 1) {
      throw new IllegalStateException("More than one activity score matched for " + this + ": " +
          String.join("; ", keys));
    } else if (keys.isEmpty()) {
      return TextConstants.NA;
    } else {
      return keys.iterator().next();
    }
  }

  private SortedSet<String> lookupPhenotypeByActivityScore(GenePhenotype gp, String activityScore) {
    SortedSet<String> rez = gp.getDiplotypes().stream()
        .filter(d -> d.getLookupKey().equals(activityScore) && d.getGeneResult() != null)
        .map(DiplotypeRecord::getGeneResult)
        .collect(Collectors.toCollection(TreeSet::new));
    if (rez.isEmpty()) {
      rez.add(TextConstants.NA);
    }
    return rez;
  }

  private static SortedSet<String> lookupActivityScoresByPhenotype(GenePhenotype gp, Collection<String> phenotypes) {
    SortedSet<String> rez = gp.getDiplotypes().stream()
        .filter(d -> phenotypes.contains(d.getGeneResult()))
        .map(DiplotypeRecord::getActivityScore)
        .filter(Objects::nonNull)
        .collect(Collectors.toCollection(TreeSet::new));
    if (rez.isEmpty()) {
      rez.add(TextConstants.NA);
    }
    return rez;
  }

  /**
   * Look up activity score based on either haplotypes or phenotype.
   * Used to check outside call activity score.
   */
  private void checkActivityScore(@Nullable String phenotype, GenePhenotype gp) {
    if (phenotype != null) {
      SortedSet<String> expected = lookupActivityScoresByPhenotype(gp, List.of(phenotype));
      if (!expected.contains(m_activityScore)) {
        m_outsideActivityScoreMismatch = String.join(", ", expected);
        return;
      }
    }
    if (!isUnknownAlleles()) {
      Map<String, Integer> lookupMap = computeLookupMap(m_allele1, m_allele2, phenotype);
      String expected = ActivityUtils.normalize(lookupKeys("activity score", gp, lookupMap));
      if (!m_activityScore.equals(expected)) {
        m_outsideActivityScoreMismatch = expected;
      }
    }
  }


  /**
   * Look up phenotype based on either haplotypes or activity score.
   * Used to check outside call phenotype.
   */
  private void checkPhenotype(GenePhenotype gp) {

    if (gp.isMatchedByActivityScore() && hasActivityScore()) {
      SortedSet<String> expected = lookupPhenotypeByActivityScore(gp, m_activityScore);
      if (!expected.contains(m_phenotypes.get(0))) {
        m_outsidePhenotypeMismatch = String.join(", ", expected);
        return;
      }
    }
    if (!isUnknownAlleles()) {
      Map<String, Integer> lookupMap = computeLookupMap(m_allele1, m_allele2, null);
      String expected = lookupPhenotypesByDiplotype(gp, lookupMap);
      if (!expected.equals(m_phenotypes.get(0))) {
        m_outsidePhenotypeMismatch = expected;
      }
    }
  }
}
