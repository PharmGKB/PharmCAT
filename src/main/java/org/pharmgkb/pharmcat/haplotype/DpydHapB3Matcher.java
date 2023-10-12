package org.pharmgkb.pharmcat.haplotype;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Objects;
import java.util.SortedMap;
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
  public static final String HAPB3_EXONIC_RSID = "rs56038477";
  public static final String HAPB3_INTRONIC_RSID = "rs75017182";

  private final SortedMap<String, SampleAllele> m_alleleMap;
  private final boolean m_isMissingHapB3;
  private List<String> m_hapB3Call;
  private MessageAnnotation m_warning;


  public DpydHapB3Matcher(Env env, SortedMap<String, SampleAllele> alleleMap,
      boolean isEffectivelyPhased) {
    m_alleleMap = alleleMap;
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

    // call HapB3
    if (hapB3ExonSample == null && hapB3IntronSample == null) {
      m_isMissingHapB3 = true;
    } else {
      m_isMissingHapB3 = false;
      if (hapB3IntronSample != null) {
        m_hapB3Call = callHapB3(hapB3Allele, hapB3IntronLocus, hapB3IntronSample);

        if (hapB3ExonSample != null) {
          if (m_hapB3Call.size() == 2) {
            List<String> exonicHapB3 = callHapB3(hapB3Allele, hapB3ExonLocus, hapB3ExonSample);
            if (!m_hapB3Call.contains("1")) {
              if (exonicHapB3.contains("1")) {
                m_warning = env.getMessage(MessageHelper.MSG_DPYD_HAPB3_INTRONIC_MISMATCH_EXONIC);
              }
            } else {
              if (exonicHapB3.size() == 2) {
                if (isEffectivelyPhased) {
                  if (!(m_hapB3Call.get(0).equals(exonicHapB3.get(0)) && m_hapB3Call.get(1).equals(exonicHapB3.get(1)))) {
                    m_warning = env.getMessage(MessageHelper.MSG_DPYD_HAPB3_INTRONIC_MISMATCH_EXONIC);
                  }
                } else {
                  if (m_hapB3Call.stream().filter(x -> x.equals("1")).count() !=
                      exonicHapB3.stream().filter(x -> x.equals("1")).count()) {
                    m_warning = env.getMessage(MessageHelper.MSG_DPYD_HAPB3_INTRONIC_MISMATCH_EXONIC);
                  }
                }
              } else {
                m_warning = env.getMessage(MessageHelper.MSG_DPYD_HAPB3_INTRONIC_MISMATCH_EXONIC);
              }
            }
          }
        }
      } else {
        m_hapB3Call = callHapB3(hapB3Allele, hapB3ExonLocus, hapB3ExonSample);
        m_warning = env.getMessage(MessageHelper.MSG_DPYD_HAPB3_EXONIC_ONLY);
      }
    }
  }


  /**
   * Calls HapB3 using a single {@link VariantLocus}.
   */
  private List<String> callHapB3(NamedAllele hapB3Allele, VariantLocus locus, SampleAllele sampleAllele) {
    List<String> rez = new ArrayList<>();
    String exonAllele = Objects.requireNonNull(hapB3Allele.getAllele(locus));
    rez.add(exonAllele.equals(sampleAllele.getAllele1()) ? "1" : "0");
    if (sampleAllele.getAllele2() != null) {
      rez.add(exonAllele.equals(sampleAllele.getAllele2()) ? "1" : "0");
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
   * Gets whether HapB3 named allele is called.
   */
  public boolean isHapB3Present() {
    return m_hapB3Call != null && m_hapB3Call.size() == 2 && m_hapB3Call.contains("1");
  }


  public int getNumHapB3Called() {
    if (m_hapB3Call == null) {
      return 0;
    }
    return (int)m_hapB3Call.stream()
        .filter(c -> c.equals("1"))
        .count();
  }

  public List<HaplotypeMatch> buildHapB3HaplotypeMatches(MatchData matchData) {
    return m_hapB3Call.stream()
        .filter(c -> c.equals("1"))
        .map(c -> new HaplotypeMatch(findHapB3(matchData)))
        .toList();
  }


  public List<DiplotypeMatch> mergePhasedDiplotypeMatch(MatchData matchData, List<DiplotypeMatch> diplotypeMatches) {

    if (!isHapB3Present()) {
      return diplotypeMatches;
    }
    if (diplotypeMatches.size() > 1) {
      throw new IllegalStateException("Should only have a single diplotype match!");
    }
    DiplotypeMatch dm = diplotypeMatches.get(0);
    if (dm.getHaplotype2() == null) {
      // should never happen
      throw new IllegalStateException("Single stranded DPYD diplotype!");
    }

    BaseMatch h1 = dm.getHaplotype1();
    BaseMatch h2 = dm.getHaplotype2();

    if (dm.getHaplotype1().getName().equals(dm.getHaplotype2().getName()) || getNumHapB3Called() == 2) {
      Collections.sort(m_hapB3Call);
      if (h1 == h2 && h2 instanceof CombinationMatch && getNumHapB3Called() == 1) {
        // NamedAlleleMatcherTest.testDpydPhasedDouble tests this code path
        h2 = new CombinationMatch((CombinationMatch)h2);
      }
      return List.of(new DiplotypeMatch(
          updateHapB3Haplotype(matchData, h1, m_hapB3Call.get(0)),
          updateHapB3Haplotype(matchData, h2, m_hapB3Call.get(1)),
          matchData
      ));
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

  private List<DiplotypeMatch> buildDiplotype(MatchData matchData, BaseMatch h1, BaseMatch h2) {
    return List.of(new DiplotypeMatch(
        updateHapB3Haplotype(matchData, h1, m_hapB3Call.get(0)),
        updateHapB3Haplotype(matchData, h2, m_hapB3Call.get(1)),
        matchData));
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
      if (sampleAllele == null) {
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

  private BaseMatch updateHapB3Haplotype(MatchData matchData, BaseMatch bm, String isHapB3) {
    if (isHapB3.equals("0")) {
      return bm;
    }

    if (bm.getSequences().size() > 1) {
      throw new IllegalStateException("Phased DPYD match should only have 1 sequence");
    }

    NamedAllele hapB3 = findHapB3(matchData);
    if (bm instanceof CombinationMatch) {
      CombinationMatch cm = null;
      for (NamedAllele c : ((CombinationMatch)bm).getComponentHaplotypes()) {
        NamedAllele component = matchData.getHaplotypes().stream()
            .filter(h -> h.getName().equals(c.getName()))
            .findAny()
            .orElseThrow(() -> new IllegalStateException("Cannot find DPYD allele '" + c.getName() + "'"));
        if (cm == null) {
          cm = new CombinationMatch(matchData.getPositions(), component, bm.getSequences().first());
        } else {
          cm.merge(component);
        }
      }
      Objects.requireNonNull(cm).merge(hapB3);
      return cm;
    }
    HaplotypeMatch hm = (HaplotypeMatch)bm;
    if (hm.getName().equals(TextConstants.REFERENCE)) {
      return new HaplotypeMatch(hapB3);
    }
    NamedAllele hap = matchData.getHaplotypes().stream()
        .filter(h -> h.getName().equals(hm.getName()))
        .findAny()
        .orElseThrow(() -> new IllegalStateException("Cannot find DPYD allele '" + hm.getName() + "'"));
    // warning: sequence will not match haplotype because sequence won't have HapB3 positions
    CombinationMatch cm = new CombinationMatch(matchData.getPositions(), hap, hm.getSequences().first());
    cm.merge(hapB3);
    return cm;
  }


  private NamedAllele findHapB3(MatchData matchData) {
    return matchData.getHaplotypes().stream()
        .filter(na -> na.getName().equals(HAPB3_ALLELE))
        .findAny()
        .orElseThrow(() -> new IllegalStateException("DPYD definition is missing HapB3 allele (" + HAPB3_ALLELE + ")"));
  }


  public List<MessageAnnotation> getWarnings() {
    if (m_warning == null) {
      return Collections.emptyList();
    }
    return List.of(m_warning);
  }
}
