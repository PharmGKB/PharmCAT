package org.pharmgkb.pharmcat.haplotype;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Objects;
import java.util.SortedMap;
import java.util.SortedSet;
import java.util.TreeSet;
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

  private final MatchData m_origData;
  private final SortedMap<String, SampleAllele> m_alleleMap;
  private final boolean m_isMissingHapB3;
  private List<String> m_hapB3IntronCall;
  private List<String> m_hapB3Call;
  private final boolean m_isHapB3Present;
  private int m_numHapB3Called;
  private MessageAnnotation m_warning;


  public DpydHapB3Matcher(Env env, SortedMap<String, SampleAllele> alleleMap, MatchData origData) {
    m_alleleMap = alleleMap;
    m_origData = origData;

    DefinitionReader definitionReader = env.getDefinitionReader();
    NamedAllele hapB3Allele = definitionReader.getHaplotypes(sf_gene).stream()
        .filter(na -> na.getName().equals(HAPB3_ALLELE))
        .findAny()
        .orElseThrow(() -> new IllegalStateException("DPYD definition is missing HapB3 allele (" + HAPB3_ALLELE + ")"));
    VariantLocus hapB3ExonLocus = definitionReader.getLocationsOfInterest().values().stream()
        .filter(vl -> HAPB3_EXONIC_RSID.equals(vl.getRsid()))
        .findAny()
        .orElseThrow(() -> new IllegalStateException("DPYD definition is missing exonic HapB3 variant " +
            HAPB3_EXONIC_RSID));
    VariantLocus hapB3IntronLocus = definitionReader.getLocationsOfInterest().values().stream()
        .filter(vl -> HAPB3_INTRONIC_RSID.equals(vl.getRsid()))
        .findAny()
        .orElseThrow(() -> new IllegalStateException("DPYD definition is missing intronic HapB3 variant " +
            HAPB3_INTRONIC_RSID));

    SampleAllele hapB3ExonSample = alleleMap.get(hapB3ExonLocus.getVcfChrPosition());
    SampleAllele hapB3IntronSample = alleleMap.get(hapB3IntronLocus.getVcfChrPosition());

    boolean isEffectivelyPhased = m_origData.isEffectivelyPhased();
    if (!isEffectivelyPhased && m_origData.isUsingPhaseSets()) {
      Integer exonPs = m_origData.getPhaseSet(hapB3ExonLocus.getPosition());
      Integer intronPs = m_origData.getPhaseSet(hapB3IntronLocus.getPosition());
      if (exonPs != null && exonPs.equals(intronPs)) {
        isEffectivelyPhased = true;
      }
    }

    // call HapB3
    if (hapB3ExonSample == null && hapB3IntronSample == null) {
      m_isMissingHapB3 = true;
    } else {
      m_isMissingHapB3 = false;
      if (hapB3IntronSample != null) {
        List<String> intronicHapB3 = callHapB3(hapB3Allele, hapB3IntronLocus, hapB3IntronSample);
        long numIntronicStrands = intronicHapB3.stream().filter(a -> !a.equals(".")).count();

        if (hapB3ExonSample == null) {
          m_hapB3IntronCall = intronicHapB3;
        } else {
          // we have both intron and exon variants
          List<String> exonicHapB3 = callHapB3(hapB3Allele, hapB3ExonLocus, hapB3ExonSample);
          long numExonicStrands = exonicHapB3.stream().filter(a -> !a.equals(".")).count();

          if (numIntronicStrands == 0) {
            m_hapB3Call = exonicHapB3;
            if (numExonicStrands != 0) {
              m_warning = env.getMessage(MessageHelper.MSG_DPYD_HAPB3_INTRONIC_MISMATCH_EXONIC);
            }
          } else if (numIntronicStrands == 1) {
            if (numExonicStrands == 0) {
              m_hapB3IntronCall = intronicHapB3;
            } else {
              // ignore if not phased
              if (isEffectivelyPhased) {
                handlePhasedCall(env, intronicHapB3, exonicHapB3, 0);
                handlePhasedCall(env, intronicHapB3, exonicHapB3, 1);
              }
            }
          } else if (numIntronicStrands == 2) {
            if (numExonicStrands == 0) {
              m_hapB3IntronCall = intronicHapB3;
            } else {
              if (isEffectivelyPhased) {
                handlePhasedCall(env, intronicHapB3, exonicHapB3, 0);
                handlePhasedCall(env, intronicHapB3, exonicHapB3, 1);
              } else {
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
                      m_warning = env.getMessage(MessageHelper.MSG_DPYD_HAPB3_INTRONIC_MISMATCH_EXONIC);
                    } else {
                      // 1 exon, 0 intron - call reference
                      // 2 exon, 0 intron - call reference
                      m_warning = env.getMessage(MessageHelper.MSG_DPYD_HAPB3_INTRONIC_MISMATCH_EXONIC);
                    }
                  }
                }
              }
            }
          }
        }
      } else {
        m_hapB3Call = callHapB3(hapB3Allele, hapB3ExonLocus, hapB3ExonSample);
        long numAlleles = m_hapB3Call.stream().filter(a -> !a.equals(".")).count();
        if (numAlleles == 2) {
          m_warning = env.getMessage(MessageHelper.MSG_DPYD_HAPB3_EXONIC_ONLY);
        }
      }
    }

    if (m_hapB3IntronCall != null) {
      m_numHapB3Called += (int)m_hapB3IntronCall.stream().filter(c -> c.equals("1")).count();
    }
    if (m_hapB3Call != null) {
      m_numHapB3Called += (int)m_hapB3Call.stream().filter(c -> c.equals("1")).count();
    }
    m_isHapB3Present = m_numHapB3Called > 0;
  }


  private void handlePhasedCall(Env env, List<String> intronicCalls, List<String> exonicCalls, int idx) {

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
        m_warning = env.getMessage(MessageHelper.MSG_DPYD_HAPB3_INTRONIC_MISMATCH_EXONIC);
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
      m_warning = new MessageAnnotation(MessageAnnotation.TYPE_NOTE, "warn.alleleCount",
          "Only found " + count + " allele for " + locus.getRsid());
    }
    return rez;
  }


  public static boolean isHapB3Rsid(String rsid) {
    return HAPB3_EXONIC_RSID.equals(rsid) || HAPB3_INTRONIC_RSID.equals(rsid);
  }


  /**
   * Gets whether HapB3 locations are present in sample.
   */
  public boolean isMissingHapB3Positions() {
    return m_isMissingHapB3;
  }

  /**
   * Gets whether HapB3 (or HapB3 intron) named allele is called.
   */
  public boolean isHapB3Present() {
    return m_isHapB3Present;
  }


  public int getNumHapB3Called() {
    return m_numHapB3Called;
  }

  /**
   * Generate {@code HaplotypeMatch}es for unphased data.
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


  public SortedSet<DiplotypeMatch> mergePhasedDiplotypeMatch(MatchData matchData, SortedSet<DiplotypeMatch> diplotypeMatches) {

    if (!m_isHapB3Present) {
      return diplotypeMatches;
    }

    SortedSet<DiplotypeMatch> finalMatches = new TreeSet<>();
    for (DiplotypeMatch dm : diplotypeMatches) {
      finalMatches.add(buildDiplotype(matchData, dm));
    }
    return finalMatches;
  }

  private DiplotypeMatch buildDiplotype(MatchData matchData, DiplotypeMatch dm) {
    if (dm.getHaplotype2() == null) {
      // should never happen
      throw new IllegalStateException("Single stranded DPYD diplotype!");
    }

    BaseMatch h1 = dm.getHaplotype1();
    BaseMatch h2 = dm.getHaplotype2();

    boolean isHapB3homozygous = m_numHapB3Called == 2 && (
        (m_hapB3Call == null || m_hapB3Call.stream().filter(c -> c.equals("1")).count() == 2) ||
            (m_hapB3IntronCall == null || m_hapB3IntronCall.stream().filter(c -> c.equals("1")).count() == 2)
    );

    if (dm.getHaplotype1().getName().equals(dm.getHaplotype2().getName()) || isHapB3homozygous) {
      // homozygous
      if (h1 == h2 && h2 instanceof CombinationMatch && m_numHapB3Called == 1) {
        // NamedAlleleMatcherTest.testDpydPhasedDouble tests this code path
        h2 = new CombinationMatch((CombinationMatch)h2);
      }
      return buildDiplotype(matchData, h1, h2);
    }

    int h1s1 = checkStrand(matchData, h1, true);
    int h1s2 = checkStrand(matchData, h1, false);
    int h2s1 = checkStrand(matchData, h2, true);
    int h2s2 = checkStrand(matchData, h2, false);

    if (h1s1 >= h1s2) {
      // hap 1 is strand 1 (but can be strand 2 if h1s1 == h1s2)
      if (h2s1 > h2s2) {
        if (h1s1 == h1s2) {
          return buildDiplotype(matchData, h2, h1);
        }
        throw new IllegalStateException("STRAND MISMATCH 1");
      }
      return buildDiplotype(matchData, h1, h2);
    }
    // hap 1 is strand 2
    if (h2s1 >= h2s2) {
      return buildDiplotype(matchData, h2, h1);
    }
    throw new IllegalStateException("STRAND MISMATCH 2");
  }

  /**
   * Build {@link DiplotypeMatch} for phased data.
   */
  private DiplotypeMatch buildDiplotype(MatchData matchData, BaseMatch h1, BaseMatch h2) {
    return new DiplotypeMatch(
        updateHapB3Haplotype(matchData, h1, 0),
        updateHapB3Haplotype(matchData, h2, 1),
        matchData);
  }


  /**
   * Calculates a score for how well a strand matches a {@link BaseMatch}.
   */
  private int checkStrand(MatchData matchData, BaseMatch bm, boolean strand1) {
    int total = 0;
    int match = 0;
    int noData = 0;
    for (VariantLocus vl : matchData.getPositions()) {
      if (isHapB3Rsid(vl.getRsid())) {
        continue;
      }
      total += 1;
      SampleAllele sa = m_alleleMap.get(vl.getVcfChrPosition());
      String sampleAllele = strand1 ? sa.getComputedAllele1() : sa.getComputedAllele2();
      if (sampleAllele == null || sampleAllele.equals(".")) {
        //System.out.println("Sample has no allele @ " + vl);
        noData += 1;
        continue;
      }
      String matchAllele = bm.getHaplotype().getAllele(vl);
      if (matchAllele == null) {
        //System.out.println(bm + " has no allele @ " + vl);
        noData += 1;
        continue;
      }
      if (!matchAllele.equals(sampleAllele)) {
        //System.out.println("Strand " + (strand1 ? "1 " : "2 ") + bm + " mismatch @ " + vl + " (" + vl.getRef() + ">" +
        //    vl.getAlts().get(0) + ") - sample is " + sampleAllele + " but looking for " + matchAllele);
        continue;
      }
      match += 1;
    }
    return match;
  }

  private BaseMatch updateHapB3Haplotype(MatchData matchData, BaseMatch bm, int alleleIndex) {
    String alleleName;
    if (m_hapB3Call != null && m_hapB3Call.size() > alleleIndex && m_hapB3Call.get(alleleIndex).equals("1")) {
      alleleName = HAPB3_ALLELE;
    } else if (m_hapB3IntronCall != null && m_hapB3IntronCall.size() > alleleIndex &&
        m_hapB3IntronCall.get(alleleIndex).equals("1")) {
      alleleName = HAPB3_INTRONIC_ALLELE;
    } else {
      return bm;
    }

    if (bm.getSequences().size() > 1) {
      throw new IllegalStateException("Phased DPYD match should only have 1 sequence");
    }

    NamedAllele hapB3 = findHapB3Allele(matchData, alleleName);
    if (bm instanceof CombinationMatch) {
      SortedSet<NamedAllele> components = new TreeSet<>();
      components.add(hapB3);
      for (NamedAllele c : ((CombinationMatch)bm).getComponentHaplotypes()) {
        NamedAllele component = matchData.getHaplotypes().stream()
            .filter(h -> h.getName().equals(c.getName()))
            .findAny()
            .orElseThrow(() -> new IllegalStateException("Cannot find DPYD allele '" + c.getName() + "'"));
        components.add(component);
      }
      return new CombinationMatch(matchData.getPositions(), bm.getSequences().first(), components, null);
    }
    HaplotypeMatch hm = (HaplotypeMatch)bm;
    if (hm.getName().equals(TextConstants.REFERENCE)) {
      return new HaplotypeMatch(hapB3);
    }
    SortedSet<NamedAllele> components = new TreeSet<>();
    components.add(hapB3);
    NamedAllele hap = matchData.getHaplotypes().stream()
        .filter(h -> h.getName().equals(hm.getName()))
        .findAny()
        .orElseThrow(() -> new IllegalStateException("Cannot find DPYD allele '" + hm.getName() + "'"));
    components.add(hap);
    // warning: the sequence will not match haplotype because the sequence won't have HapB3 positions
    return new CombinationMatch(matchData.getPositions(), hm.getSequences().first(), components, null);
  }


  private NamedAllele findHapB3Allele(MatchData matchData, String name) {
    return matchData.getHaplotypes().stream()
        .filter(na -> na.getName().equals(name))
        .findAny()
        .orElseThrow(() -> new IllegalStateException("DPYD definition is missing HapB3 allele (" + name + ")"));
  }


  public List<MessageAnnotation> getWarnings() {
    if (m_warning == null) {
      return Collections.emptyList();
    }
    return List.of(m_warning);
  }
}
