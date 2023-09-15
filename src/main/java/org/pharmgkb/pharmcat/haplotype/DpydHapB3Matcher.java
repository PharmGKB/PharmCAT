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

    if (dm.getHaplotype1().getName().equals(dm.getHaplotype2().getName())) {
      Collections.sort(m_hapB3Call);
      BaseMatch h1 = dm.getHaplotype1();
      BaseMatch h2 = dm.getHaplotype2();
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

    if (checkStrand(matchData, dm.getHaplotype1(), true)) {
      if (!checkStrand(matchData, dm.getHaplotype2(), false)) {
        throw new IllegalStateException("STRAND MISMATCH 1");
      }
      return List.of(new DiplotypeMatch(
          updateHapB3Haplotype(matchData, dm.getHaplotype1(), m_hapB3Call.get(0)),
          updateHapB3Haplotype(matchData, dm.getHaplotype2(), m_hapB3Call.get(1)),
          matchData
      ));
    } else {
      if (!checkStrand(matchData, dm.getHaplotype1(), false)) {
        throw new IllegalStateException("STRAND MISMATCH 2");
      }
      if (!checkStrand(matchData, dm.getHaplotype2(), true)) {
        throw new IllegalStateException("STRAND MISMATCH 2");
      }
      return List.of(new DiplotypeMatch(
          updateHapB3Haplotype(matchData, dm.getHaplotype2(), m_hapB3Call.get(0)),
          updateHapB3Haplotype(matchData, dm.getHaplotype1(), m_hapB3Call.get(1)),
          matchData
      ));
    }
  }

  private boolean checkStrand(MatchData matchData, BaseMatch bm, boolean strand1) {
    for (VariantLocus vl : matchData.getPositions()) {
      if (isHapB3Rsid(vl.getRsid())) {
        continue;
      }
      SampleAllele sa = m_alleleMap.get(vl.getVcfChrPosition());
      String expectedAllele = strand1 ? sa.getAllele1() : sa.getAllele2();
      if (expectedAllele == null) {
        //System.out.println("Sample has no allele @ " + vl);
        continue;
      }
      String actualAllele = bm.getHaplotype().getAllele(vl);
      if (actualAllele == null) {
        //System.out.println(bm + " has no allele @ " + vl);
        continue;
      }
      if (!actualAllele.equals(expectedAllele)) {
        return false;
      }
    }
    return true;
  }

  private BaseMatch updateHapB3Haplotype(MatchData matchData, BaseMatch bm, String isHapB3) {
    if (isHapB3.equals("0")) {
      return bm;
    }
    if (bm instanceof CombinationMatch) {
      CombinationMatch cm = (CombinationMatch)bm;
      cm.merge(findHapB3(matchData));
      return cm;
    }
    HaplotypeMatch hm = (HaplotypeMatch)bm;
    if (hm.getName().equals("Reference")) {
      return isHapB3.equals("1") ? new HaplotypeMatch(findHapB3(matchData)) : bm;
    }
    if (hm.getSequences().size() > 1) {
      throw new IllegalStateException("Phased DPYD haplotype match should only have 1 sequence");
    }
    NamedAllele hap = matchData.getHaplotypes().stream()
        .filter(h -> h.getName().equals(hm.getName()))
        .findAny()
        .orElseThrow(() -> new IllegalStateException("Cannot find DPYD allele '" + hm.getName() + "'"));
    // warning: sequence will not match haplotype because sequence won't have HapB3 positions
    CombinationMatch cm = new CombinationMatch(matchData.getPositions(), hap, hm.getSequences().first());
    cm.merge(findHapB3(matchData));
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
