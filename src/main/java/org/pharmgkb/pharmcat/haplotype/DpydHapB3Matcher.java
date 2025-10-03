package org.pharmgkb.pharmcat.haplotype;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Objects;
import java.util.Set;
import java.util.SortedMap;
import java.util.SortedSet;
import java.util.TreeSet;
import com.google.common.base.Preconditions;
import org.jspecify.annotations.Nullable;
import org.pharmgkb.pharmcat.Env;
import org.pharmgkb.pharmcat.definition.DefinitionReader;
import org.pharmgkb.pharmcat.definition.model.NamedAllele;
import org.pharmgkb.pharmcat.definition.model.VariantLocus;
import org.pharmgkb.pharmcat.haplotype.model.BaseMatch;
import org.pharmgkb.pharmcat.haplotype.model.CombinationMatch;
import org.pharmgkb.pharmcat.haplotype.model.DiplotypeMatch;
import org.pharmgkb.pharmcat.haplotype.model.HaplotypeMatch;
import org.pharmgkb.pharmcat.reporter.MessageHelper;
import org.pharmgkb.pharmcat.reporter.TextConstants;
import org.pharmgkb.pharmcat.reporter.model.MessageAnnotation;


/**
 * This class handles DPYD's special algorithm for dealing with the HapB3 named allele.
 *
 * @author Mark Woon
 */
public class DpydHapB3Matcher {
  private static final String sf_gene = "DPYD";
  public static final String HAPB3_ALLELE = "c.1129-5923C>G, c.1236G>A (HapB3)";
  public static final String HAPB3_INTRONIC_ALLELE = "c.1129-5923C>G";
  public static final String HAPB3_EXONIC_RSID = "rs56038477";
  public static final String HAPB3_INTRONIC_RSID = "rs75017182";
  public static final String HAPB3_EXONIC_PARTIAL = "g.97573863C>T";
  public static final String HAPB3_INTRONIC_PARTIAL = "g.97579893G>C";
  public static final long HAPB3_EXONIC_POSITION = 97573863L;
  public static final long HAPB3_INTRONIC_POSITION = 97579893;

  private final Env m_env;
  private final MatchData m_origData;
  private final SortedMap<String, SampleAllele> m_alleleMap;
  // calculated in constructor
  private final VariantLocus m_hapB3ExonLocus;
  private final VariantLocus m_hapB3IntronLocus;
  private final boolean m_hasHapB3Variants;
  private final boolean m_hasNonHapB3Variants;

  private @Nullable List<String> m_hapB3IntronCall;
  private @Nullable List<String> m_hapB3Call;
  private int m_numHapB3Called;
  private final Set<MessageAnnotation> m_warnings = new HashSet<>();


  public DpydHapB3Matcher(Env env, SortedMap<String, SampleAllele> alleleMap, MatchData origData) {
    m_env = env;
    m_origData = origData;
    m_alleleMap = alleleMap;

    DefinitionReader definitionReader = env.getDefinitionReader();
    m_hapB3ExonLocus = definitionReader.getLocationsOfInterest().values().stream()
        .filter(vl -> HAPB3_EXONIC_RSID.equals(vl.getRsid()))
        .findAny()
        .orElseThrow(() -> new IllegalStateException("DPYD definition is missing exonic HapB3 variant " +
            HAPB3_EXONIC_RSID));
    m_hapB3IntronLocus = definitionReader.getLocationsOfInterest().values().stream()
        .filter(vl -> HAPB3_INTRONIC_RSID.equals(vl.getRsid()))
        .findAny()
        .orElseThrow(() -> new IllegalStateException("DPYD definition is missing intronic HapB3 variant " +
            HAPB3_INTRONIC_RSID));

    SampleAllele hapB3ExonSample = alleleMap.get(m_hapB3ExonLocus.getVcfChrPosition());
    SampleAllele hapB3IntronSample = alleleMap.get(m_hapB3IntronLocus.getVcfChrPosition());
    m_hasNonHapB3Variants = alleleMap.values().stream()
        .filter(sa -> sa != hapB3ExonSample && sa != hapB3IntronSample)
        .map(sa -> sa.getGt().replaceAll("[|/0]", "").length())
        .anyMatch(l -> l > 0);
    m_hasHapB3Variants = alleleMap.values().stream()
        .filter(sa -> sa == hapB3ExonSample || sa == hapB3IntronSample)
        .map(sa -> sa.getGt().replaceAll("[|/0]", "").length())
        .anyMatch(l -> l > 0);
  }


  /**
   * Calls HapB3 {@link HaplotypeMatch}es for unphased data.
   * While this code is not unphased-data specific, it should only be applied to unphased data since it is better to go
   * through the {@link DiplotypeMatch} code path if possible.
   */
  public void callHapB3HaplotypeMatches() {

    boolean isEffectivelyPhased = m_origData.isEffectivelyPhased();
    if (!isEffectivelyPhased && m_origData.isUsingPhaseSets()) {
      Integer exonPs = m_origData.getPhaseSet(m_hapB3ExonLocus.getPosition());
      Integer intronPs = m_origData.getPhaseSet(m_hapB3IntronLocus.getPosition());
      if (exonPs != null && exonPs.equals(intronPs)) {
        isEffectivelyPhased = true;
      }
    }

    NamedAllele hapB3Allele = m_env.getDefinitionReader().getHaplotypes(sf_gene).stream()
        .filter(na -> na.getName().equals(HAPB3_ALLELE))
        .findAny()
        .orElseThrow(() -> new IllegalStateException("DPYD definition is missing HapB3 allele (" + HAPB3_ALLELE + ")"));

    // call HapB3
    SampleAllele hapB3ExonSample = m_alleleMap.get(m_hapB3ExonLocus.getVcfChrPosition());
    SampleAllele hapB3IntronSample = m_alleleMap.get(m_hapB3IntronLocus.getVcfChrPosition());
    if (hapB3ExonSample != null || hapB3IntronSample != null) {
      if (hapB3IntronSample != null) {
        List<String> intronicHapB3 = callHapB3(hapB3Allele, m_hapB3IntronLocus, hapB3IntronSample);
        long numIntronicStrands = intronicHapB3.stream().filter(a -> !a.equals(".")).count();

        if (hapB3ExonSample == null) {
          m_hapB3IntronCall = intronicHapB3;
        } else {
          // we have both intron and exon variants
          List<String> exonicHapB3 = callHapB3(hapB3Allele, m_hapB3ExonLocus, hapB3ExonSample);
          long numExonicStrands = exonicHapB3.stream().filter(a -> !a.equals(".")).count();

          if (numIntronicStrands == 0) {
            m_hapB3Call = exonicHapB3;
            if (numExonicStrands != 0) {
              m_warnings.add(m_env.getMessage(MessageHelper.MSG_DPYD_HAPB3_INTRONIC_MISMATCH_EXONIC));
            }
          } else if (numIntronicStrands == 1) {
            if (numExonicStrands == 0) {
              m_hapB3IntronCall = intronicHapB3;
            } else {
              // ignore if not phased
              if (isEffectivelyPhased) {
                handlePhasedCall(intronicHapB3, exonicHapB3, 0);
                handlePhasedCall(intronicHapB3, exonicHapB3, 1);
              }
            }
          } else if (numIntronicStrands == 2) {
            if (numExonicStrands == 0) {
              m_hapB3IntronCall = intronicHapB3;
            } else {
              if (isEffectivelyPhased) {
                handlePhasedCall(intronicHapB3, exonicHapB3, 0);
                handlePhasedCall(intronicHapB3, exonicHapB3, 1);
              } else {
                //noinspection StatementWithEmptyBody
                if (numExonicStrands == 1) {
                  // ignore this case
                } else if (numExonicStrands == 2) {
                  m_hapB3IntronCall = new ArrayList<>();
                  m_hapB3Call = new ArrayList<>();
                  long numIntronic = intronicHapB3.stream().filter(c -> c.equals("1")).count();
                  long numExonic = exonicHapB3.stream().filter(c -> c.equals("1")).count();

                  if (numIntronic == numExonic) {
                    for (int x = 0; x < numIntronic; x += 1) {
                      m_hapB3Call.add("1");
                    }
                  } else if (numIntronic > numExonic) {
                    // 2 intron, 1 exon
                    // 2 intron, 0 exon
                    // 1 intron, 0 exon
                    m_hapB3IntronCall.add("1");
                    if (numExonic == 1) {
                      m_hapB3Call.add("1");
                    } else if (numIntronic == 2) {
                      m_hapB3IntronCall.add("1");
                    }
                  } else {
                    // intronic call trumps exonic call
                    if (numIntronic == 1) {
                      // 2 exons, 1 intron - call ref/hapB3
                      m_hapB3Call.add("1");
                      m_warnings.add(m_env.getMessage(MessageHelper.MSG_DPYD_HAPB3_INTRONIC_MISMATCH_EXONIC));
                    } else {
                      // 1 exon, 0 intron - call reference
                      // 2 exon, 0 intron - call reference
                      m_warnings.add(m_env.getMessage(MessageHelper.MSG_DPYD_HAPB3_INTRONIC_MISMATCH_EXONIC));
                    }
                  }
                }
              }
            }
          }
        }
      } else {
        m_hapB3Call = callHapB3(hapB3Allele, m_hapB3ExonLocus, hapB3ExonSample);
        long numAlleles = m_hapB3Call.stream().filter(a -> !a.equals(".")).count();
        if (numAlleles == 2) {
          m_warnings.add(m_env.getMessage(MessageHelper.MSG_DPYD_HAPB3_EXONIC_ONLY));
        }
      }
    }

    if (m_hapB3IntronCall != null) {
      m_numHapB3Called += (int)m_hapB3IntronCall.stream().filter(c -> c.equals("1")).count();
    }
    if (m_hapB3Call != null) {
      m_numHapB3Called += (int)m_hapB3Call.stream().filter(c -> c.equals("1")).count();
    }
  }

  private void handlePhasedCall(List<String> intronicCalls, List<String> exonicCalls, int idx) {

    String intronicCall = ".";
    if (intronicCalls.size() > idx) {
      intronicCall = intronicCalls.get(idx);
    }
    String exonicCall = ".";
    if (exonicCalls.size() > idx) {
      exonicCall = exonicCalls.get(idx);
    }

    if (m_hapB3IntronCall == null) {
      m_hapB3IntronCall = new ArrayList<>();
    }
    if (m_hapB3Call == null) {
      m_hapB3Call = new ArrayList<>();
    }

    if (intronicCall.equals(".") && exonicCall.equals(".")) {
      m_hapB3IntronCall.add(".");
      m_hapB3Call.add(".");
    }

    // intronic call trumps exonic call
    if (intronicCall.equals("0")) {
      m_hapB3IntronCall.add("0");
      m_hapB3Call.add("0");
      if (!exonicCall.equals("0")) {
        m_warnings.add(m_env.getMessage(MessageHelper.MSG_DPYD_HAPB3_INTRONIC_MISMATCH_EXONIC));
      }
    } else if (intronicCall.equals("1")) {
      if (exonicCall.equals("0")) {
        m_hapB3IntronCall.add("1");
        m_hapB3Call.add("0");
      } else if (exonicCall.equals(".")) {
        m_hapB3IntronCall.add("1");
        m_hapB3Call.add(".");
      } else {
        m_hapB3IntronCall.add("0");
        m_hapB3Call.add("1");
      }
    } else {
      m_hapB3IntronCall.add(".");
      m_hapB3Call.add(".");
    }
  }

  /**
   * Calls HapB3 using a single {@link VariantLocus}.
   */
  private List<String> callHapB3(NamedAllele hapB3Allele, VariantLocus locus, SampleAllele sampleAllele) {
    List<String> rez = new ArrayList<>();
    String exonAllele = Objects.requireNonNull(hapB3Allele.getAllele(locus));
    int count = 0;
    if (sampleAllele.getAllele1() == null) {
      if (sampleAllele.isPhased()) {
        // we only care if this is phased
        rez.add(".");
      }
    } else {
      rez.add(exonAllele.equals(sampleAllele.getAllele1()) ? "1" : "0");
      count += 1;
    }
    if (sampleAllele.getAllele2() != null) {
      rez.add(exonAllele.equals(sampleAllele.getAllele2()) ? "1" : "0");
      count += 1;
    }
    if (count < 2) {
      // this warning can get overwritten!
      // but we don't really care
      m_warnings.add(new MessageAnnotation(MessageAnnotation.TYPE_NOTE, "warn.alleleCount",
          "Only found " + count + " allele for " + locus.getRsid()));
    }
    return rez;
  }

  /**
   * Gets whether HapB3 variants are present in the sample.
   */
  public boolean hasHapB3Variants() {
    return m_hasHapB3Variants;
  }

  /**
   * Gets whether non-HapB3 variants are present in the sample.
   */
  public boolean hasNonHapB3Variants() {
    return m_hasNonHapB3Variants;
  }

  /**
   * Gets whether HapB3 (or HapB3 intron) named allele is called.
   */
  public boolean isHapB3Present() {
    return m_numHapB3Called > 0;
  }


  public int getNumHapB3Called() {
    return m_numHapB3Called;
  }

  /**
   * Generate {@link HaplotypeMatch}es for unphased data.
   */
  public List<HaplotypeMatch> buildHapB3HaplotypeMatches() {
    List<HaplotypeMatch> matches = new ArrayList<>();
    if (m_hapB3Call != null) {
      matches.addAll(m_hapB3Call.stream()
            .filter(c -> c.equals("1"))
            .map(c -> new HaplotypeMatch(findHapB3Allele(m_origData, HAPB3_ALLELE)))
            .toList());
    }
    if (m_hapB3IntronCall != null) {
      matches.addAll(m_hapB3IntronCall.stream()
          .filter(c -> c.equals("1"))
          .map(c -> new HaplotypeMatch(findHapB3Allele(m_origData, HAPB3_INTRONIC_ALLELE)))
          .toList());
    }
    return matches;
  }


  /**
   * Merges the HapB3 results back into the diplotype matches.
   * This requires the full MatchData (i.e., including HapB3 positions).
   */
  public SortedSet<DiplotypeMatch> mergePhasedHapB3Call(MatchData matchData,
      SortedSet<DiplotypeMatch> diplotypeMatches) {

    if (!m_hasHapB3Variants) {
      return diplotypeMatches;
    }

    SortedSet<DiplotypeMatch> finalMatches = new TreeSet<>();
    for (DiplotypeMatch dm : diplotypeMatches) {
      BaseMatch hap1 = dm.getHaplotype1();
      hap1 = callPhasedHapB3(hap1.getSequences().first(), matchData, hap1);
      BaseMatch hap2 = dm.getHaplotype2();
      if (hap2 == null) {
        // this should never happen
        throw new IllegalStateException("Haplotype 2 is null");
      }
      hap2 = callPhasedHapB3(hap2.getSequences().first(), matchData, hap2);
      finalMatches.add(new DiplotypeMatch(hap1, hap2, matchData));
    }
    return finalMatches;
  }

  /**
   * Adds the HapB3 call to phased sequence data that is Ref/Ref everywhere else.
   * This requires the full MatchData (i.e., including HapB3 positions).
   */
  public SortedSet<DiplotypeMatch> addPhasedHapB3CallToRef(MatchData matchData) {
    Preconditions.checkState(m_hasHapB3Variants && !m_hasNonHapB3Variants);

    List<BaseMatch> haps = new ArrayList<>();
    for (String seq : matchData.getPermutations()) {
      haps.add(callPhasedHapB3(seq, matchData, null));
    }

    SortedSet<DiplotypeMatch> finalMatches = new TreeSet<>();
    if (haps.size() == 2) {
      finalMatches.add(new DiplotypeMatch(haps.get(0), haps.get(1), matchData));
    } else {
      finalMatches.add(new DiplotypeMatch(haps.get(0), haps.get(0), matchData));
    }
    return finalMatches;
  }

  /**
   * Call HapB3 for a phased sequence.
   */
  private BaseMatch callPhasedHapB3(String seq, MatchData matchData, @Nullable BaseMatch baseMatch) {
    String allele = matchData.getAllele(seq, m_hapB3IntronLocus.getPosition());
    boolean hasIntronLocus = false;
    boolean hasIntron = false;
    if (allele != null && !allele.equals(".")) {
      hasIntronLocus = true;
      hasIntron = !m_hapB3IntronLocus.getRef().equals(allele);
    }
    allele = matchData.getAllele(seq, m_hapB3ExonLocus.getPosition());
    boolean hasExon = allele != null && !allele.equals(".") && !m_hapB3ExonLocus.getRef().equals(allele);

    if (hasIntron) {
      if (hasExon) {
        return buildMatch(matchData, findHapB3Allele(matchData, HAPB3_ALLELE), baseMatch);
      } else {
        return buildMatch(matchData, findHapB3Allele(matchData, HAPB3_INTRONIC_ALLELE), baseMatch);
      }
    } else {
      if (hasIntronLocus) {
        if (hasExon) {
          m_warnings.add(m_env.getMessage(MessageHelper.MSG_DPYD_HAPB3_INTRONIC_MISMATCH_EXONIC));
        }
      } else {
        if (hasExon) {
          m_warnings.add(m_env.getMessage(MessageHelper.MSG_DPYD_HAPB3_EXONIC_ONLY));
          return buildMatch(matchData, findHapB3Allele(matchData, HAPB3_ALLELE), baseMatch);
        }
      }
      return buildMatch(matchData, findHapB3Allele(matchData, TextConstants.REFERENCE), baseMatch);
    }
  }

  /**
   * Builds a match with the additional {@code hapB3Allele}.
   */
  private BaseMatch buildMatch(MatchData matchData, NamedAllele hapB3Allele, @Nullable BaseMatch baseMatch) {

    if (baseMatch == null) {
      return new HaplotypeMatch(hapB3Allele);
    }
    if (hapB3Allele.getName().equals(TextConstants.REFERENCE)) {
      return baseMatch;
    }

    if (baseMatch instanceof CombinationMatch cm) {
      SortedSet<NamedAllele> components = new TreeSet<>();
      components.add(hapB3Allele);
      for (NamedAllele c : cm.getComponentHaplotypes()) {
        NamedAllele component = matchData.getHaplotypes().stream()
            .filter(h -> h.getName().equals(c.getName()))
            .findAny()
            .orElseThrow(() -> new IllegalStateException("Cannot find DPYD allele '" + c.getName() + "'"));
        components.add(component);
      }
      return new CombinationMatch(matchData.getPositions(), cm.getSequences().first(), components, null);
    }

    HaplotypeMatch hm = (HaplotypeMatch)baseMatch;
    if (hm.getName().equals(TextConstants.REFERENCE)) {
      return new HaplotypeMatch(hapB3Allele);
    }
    SortedSet<NamedAllele> components = new TreeSet<>();
    components.add(hapB3Allele);
    NamedAllele hap = matchData.getHaplotypes().stream()
        .filter(h -> h.getName().equals(hm.getName()))
        .findAny()
        .orElseThrow(() -> new IllegalStateException("Cannot find DPYD allele '" + hm.getName() + "'"));
    components.add(hap);
    return new CombinationMatch(matchData.getPositions(), hm.getSequences().first(), components, null);
  }


  /**
   * Calling combinations without HapB3 {@link NamedAllele}s means HapB3 alleles will be called as partials.
   * This replaces those partial calls with the correct HapB3 {@code NamedAllele}s.
   */
  public SortedSet<DiplotypeMatch> fixPartials(MatchData matchData,
      SortedSet<DiplotypeMatch> diplotypeMatches) {

    SortedSet<DiplotypeMatch> updatedMatches = new TreeSet<>();
    // has matches, add HapB3 call
    for (DiplotypeMatch dm : diplotypeMatches) {
      boolean modified = false;
      BaseMatch hap1 = dm.getHaplotype1();
      if (hasHapB3Partial(hap1)) {
        hap1 = removePartials(matchData, (CombinationMatch)hap1);
        modified = true;
      }
      BaseMatch hap2 = Objects.requireNonNull(dm.getHaplotype2());
      if (hasHapB3Partial(hap2)) {
        hap2 = removePartials(matchData, (CombinationMatch)hap2);
        modified = true;
      }
      if (modified) {
        updatedMatches.add(new DiplotypeMatch(hap1, hap2, matchData));
      } else {
        updatedMatches.add(dm);
      }
    }

    return updatedMatches;
  }


  private boolean hasHapB3Partial(BaseMatch bm) {
    return bm instanceof CombinationMatch cm && cm.hasPartials() &&
        (cm.getName().contains(HAPB3_INTRONIC_PARTIAL) || cm.getName().contains(HAPB3_EXONIC_PARTIAL));
  }

  private BaseMatch removePartials(MatchData matchData, CombinationMatch cm) {
    Map<Long, String> partials = cm.getPartials();
    partials.remove(HAPB3_EXONIC_POSITION);
    partials.remove(HAPB3_INTRONIC_POSITION);
    if (partials.isEmpty()) {
      if (cm.getComponentHaplotypes().size() == 1 &&
          cm.getComponentHaplotypes().first().getName().equals(TextConstants.REFERENCE)) {
        HaplotypeMatch refMatch = new HaplotypeMatch(cm.getComponentHaplotypes().first());
        refMatch.addSequence(cm.getSequences().first());
        return refMatch;
      }
    }
    return new CombinationMatch(matchData.getPositions(), cm.getSequences()
        .first(), cm.getComponentHaplotypes(), partials);
  }


  private NamedAllele findHapB3Allele(MatchData matchData, String name) {
    return matchData.getHaplotypes().stream()
        .filter(na -> na.getName().equals(name))
        .findAny()
        .orElseThrow(() -> new IllegalStateException("DPYD definition is missing HapB3 allele (" + name + ")"));
  }


  public List<MessageAnnotation> getWarnings() {
    return new ArrayList<>(m_warnings);
  }
}
