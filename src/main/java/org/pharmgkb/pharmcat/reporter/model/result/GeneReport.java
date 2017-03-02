package org.pharmgkb.pharmcat.reporter.model.result;

import java.util.ArrayList;
import java.util.Collection;
import java.util.List;
import java.util.Objects;
import java.util.Set;
import java.util.SortedSet;
import java.util.TreeSet;
import javax.annotation.Nonnull;
import org.pharmgkb.pharmcat.UnexpectedStateException;
import org.pharmgkb.pharmcat.definition.PhenotypeMap;
import org.pharmgkb.pharmcat.definition.model.GenePhenotype;
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
  private static final String NA = "N/A";

  private String m_gene;
  private Set<String> m_diplotypes = new TreeSet<>();
  private Set<String> m_haplotypes = new TreeSet<>();
  private Set<String> m_uncalledHaplotypes;
  private SortedSet<Variant> m_variants = new TreeSet<>();
  private List<PharmcatException> m_exceptList = new ArrayList<>();
  private MatchData m_matchData;
  private List<DrugLink> m_relatedDrugs = new ArrayList<>();
  private SortedSet<String> m_functions = new TreeSet<>();
  private Set<String> m_phenotypes = new TreeSet<>();
  private boolean m_astrolabeCall = false;

  /**
   * public constructor
   */
  public GeneReport(@Nonnull String geneSymbol) {
    m_gene = geneSymbol;
    addDip(UNCALLED);
    addPhenotype(NA);
    addFunction(NA);
  }

  /**
   * public constructor
   * @param call a {@link GeneCall} that has been made by the {@link NamedAlleleMatcher}
   */
  public void setCallData(@Nonnull GeneCall call, PhenotypeMap phenotypeMap) {
    m_gene = call.getGene();
    if (!call.getDiplotypes().isEmpty()) {
      call.getDiplotypes().forEach(d -> addDip(d.getName()));
    }
    m_variants.addAll(call.getVariants());
    m_matchData = call.getMatchData();

    GenePhenotype genePhenotype = phenotypeMap.lookup(m_gene);

    // for SLCO1B1 just do a simple call and don't do a regular match
    if (m_gene.equals("SLCO1B1")) {
      callSlco1b1(genePhenotype);
      return;
    }

    // do a regular copy of diplotype, function, and phenotype
    m_uncalledHaplotypes = call.getUncallableHaplotypes();
    if (call.getDiplotypes() != null && call.getDiplotypes().size() > 0) {
      call.getDiplotypes().forEach(d -> {
        addFunction(d.getFunction());
        if (genePhenotype != null) {
          addPhenotype(genePhenotype.makePhenotype(d.getName()));
        } else {
          addPhenotype("N/A");
        }
      });
    }

    // UGT1A1 get a second chance at calling based on one position
    else {
      if (m_gene.equals("UGT1A1")) {
        callUgt1a1(genePhenotype);
      }
    }
  }

  public void setAstrolabeData(@Nonnull AstrolabeCall call, @Nonnull PhenotypeMap phenotypeMap) {
    m_astrolabeCall = true;
    m_gene = call.getGene();
    if (call.getDiplotypes() != null) {
      call.getDiplotypes().forEach(d -> {
        addDip(d);
        addFunction(phenotypeMap.lookup(m_gene).makeFunction(d));
        addPhenotype(phenotypeMap.lookup(m_gene).makePhenotype(d));
      });
    }
  }

  /**
   * Calls the allele and phenotype data specifically for SLCO1B1 which is a special case using only one variant.
   * @param genePhenotype {@link GenePhenotype} information for SLCO1B1
   */
  private void callSlco1b1(GenePhenotype genePhenotype) {
    Variant variant = m_variants.stream().filter(v -> v.getRsid() != null && v.getRsid().equals("rs4149056")).reduce((a,b) -> {
      throw new UnexpectedStateException("more than one variant found");
    }).orElse(null);

    if (variant != null) {
      String dip = null;
      switch (variant.getVcfCall()) {
        case "T|T":
          dip = "*1A/*1A";
          break;
        case "T|C":
        case "C|T":
          dip = "*1A/*5";
          break;
        case "C|C":
          dip = "*5/*5";
          break;
      }
      addDip(dip);
      addFunction(genePhenotype.makeFunction(dip));
      addPhenotype(genePhenotype.makePhenotype(dip));
    } else {
      addDip(UNCALLED);
    }
  }

  /**
   * Calls the UGT1A1 gene in the case where the matcher calls no alleles but rs887829 is available.
   * @param genePhenotype {@link GenePhenotype} information for UGT1A1
   */
  private void callUgt1a1(GenePhenotype genePhenotype) {
    Variant variant = m_variants.stream().filter(v -> v.getRsid() != null && v.getRsid().equals("rs887829")).reduce((a,b) -> {
      throw new UnexpectedStateException("more than one variant found");
    }).orElse(null);

    if (variant != null) {
      String dip = null;
      switch (variant.getVcfCall()) {
        case "T|T":
          dip = "*80/*80";
          break;
        case "T|C":
        case "C|T":
          dip = "*1/*80";
          break;
        case "C|C":
          dip = "*1/*1";
          break;
      }
      addDip(dip);
      addFunction(genePhenotype.makeFunction(dip));
      addPhenotype(genePhenotype.makePhenotype(dip));
    } else {
      addDip(UNCALLED);
    }
  }

  /**
   * Diplotypes that have been called for this gene
   */
  public Set<String> getDips(){
    return m_diplotypes;
  }

  private void addDip(String dip) {
    if (m_diplotypes.size() == 1 && m_diplotypes.contains(UNCALLED)) {
      m_diplotypes.clear();
    }

    String cleanDip = m_gene.equals("CYP2C19") ? dip.replaceAll("\\*4[AB]", "*4") : dip;

    m_diplotypes.add(cleanDip);

    if (cleanDip.contains("/")) {
      String[] haps = cleanDip.split("/");
      m_haplotypes.add(haps[0]);
      m_haplotypes.add(haps[1]);
    }
  }

  private void addFunction(String function) {
    if (m_functions.size() == 1 && m_functions.contains(NA)) {
      m_functions.clear();
    }
    m_functions.add(function);
  }

  private void addPhenotype(String phenotype) {
    if (m_phenotypes.size() == 1 && m_phenotypes.contains(NA)) {
      m_phenotypes.clear();
    }
    m_phenotypes.add(phenotype);
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
  }

  private void addException(PharmcatException except) {
    MatchLogic match = except.getMatches();

    if (match.getGene().equals(m_gene)) {
      boolean critHap = match.getHapsCalled().isEmpty() || match.getHapsCalled().stream().allMatch(h -> m_haplotypes.contains(h));
      boolean critUnmatchedHap = match.getHapsMissing().isEmpty() || match.getHapsMissing().stream().allMatch(h -> m_uncalledHaplotypes.contains(h));
      boolean critDip = match.getDips().isEmpty() || match.getDips().stream().allMatch(d -> m_diplotypes.contains(d));
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

  public MatchData getMatchData() {
    return m_matchData;
  }

  public Set<String> getUncalledHaplotypes() {
    return m_uncalledHaplotypes;
  }

  /**
   * Gets the functions in the form of "Two no function alleles"
   */
  public SortedSet<String> getFunctions() {
    return m_functions;
  }

  /**
   * Gets the gene phenotype in the form of "Intermediate Metabolizer"
   */
  public Set<String> getPhenotypes() {
    return m_phenotypes;
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

  public boolean isAstrolabeCall() {
    return m_astrolabeCall;
  }
}

