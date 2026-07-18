package org.pharmgkb.pharmcat.haplotype.benchmark;

import java.util.ArrayList;
import java.util.List;
import org.pharmgkb.pharmcat.haplotype.DpydHapB3Matcher;


/**
 * Static catalog of benchmark scenarios that together exercise every {@code NamedAlleleMatcher} code path:
 * <ul>
 *   <li>Standard-gene reference/single-het/multi-het (phased and unphased).</li>
 *   <li>Combination mode across gene sizes (CYP2C9, CYP2D6, CFTR).</li>
 *   <li>Lowest-function gene paths (DPYD, RYR1) including phase-set constraints.</li>
 *   <li>DPYD HapB3 code path variants (HapB3-only, HapB3 + others, intron/exon mismatch).</li>
 *   <li>Missing-position and undocumented-variation branches.</li>
 * </ul>
 */
final class Scenarios {

  private Scenarios() {}


  static List<Scenario> all() {
    List<Scenario> list = new ArrayList<>();

    // --- Standard-gene path (non-lowest-function) ---
    list.add(allRef("CYP2C19"));
    list.add(allRef("CYP2D6"));
    list.add(singleHet("CYP2C19"));
    list.add(singleHet("CYP2D6"));
    list.add(multiHet("CYP2C19", 5, false));
    list.add(multiHet("CYP2D6", 10, false));
    list.add(multiHet("CYP2C19", 5, true));
    list.add(multiHet("CYP2D6", 10, true));
    list.add(missingPositions("CYP2D6", 5, 3));

    // --- Combination research mode (CombinationMatcher) ---
    list.add(combinations("CYP2C9", 4));
    list.add(combinations("CYP2D6", 6));
    list.add(combinations("CFTR", 4));

    // --- Lowest-function-gene path ---
    list.add(allRef("DPYD"));
    list.add(allRef("RYR1"));
    list.add(multiHet("DPYD", 4, true));    // Stage 1: exact
    list.add(multiHet("DPYD", 4, false));   // Stage 3: fallback to haplotypes
    list.add(multiHet("RYR1", 6, false));   // Stage 3 stress on largest gene
    list.add(phaseSets("DPYD", 6));         // Stage 2: phase-set constrained
    list.add(phaseSets("RYR1", 8));

    // --- DPYD HapB3 (5-stage path) ---
    list.add(dpydHapB3Only(true));
    list.add(dpydHapB3WithOthers(true));
    list.add(dpydHapB3WithOthers(false));
    list.add(dpydHapB3IntronExonMismatch());

    return list;
  }


  // -----------------
  // Scenario builders
  // -----------------

  private static Scenario allRef(String gene) {
    return new Scenario("all-ref/" + gene, gene, false, true, true, b -> { /* default 0/0 */ });
  }

  private static Scenario singleHet(String gene) {
    return new Scenario("single-het/" + gene, gene, false, true, true,
        b -> b.set(0, "0/1"));
  }

  private static Scenario multiHet(String gene, int hetCount, boolean phased) {
    String tag = phased ? "phased" : "unphased";
    String gt = phased ? "0|1" : "0/1";
    return new Scenario("multi-het-" + tag + "/" + gene + "/" + hetCount, gene, false, true, true, b -> {
      if (phased) {
        b.setDefaultPhased(true);
      }
      int n = Math.min(hetCount, b.size());
      for (int i = 0; i < n; i += 1) {
        b.set(i, gt);
      }
    });
  }

  private static Scenario missingPositions(String gene, int missingCount, int hetCount) {
    return new Scenario("missing/" + gene + "/miss=" + missingCount + "/het=" + hetCount,
        gene, false, true, true, b -> {
      int missN = Math.min(missingCount, b.size());
      for (int i = 0; i < missN; i += 1) {
        b.set(i, "./.");
      }
      int hetN = Math.min(hetCount, b.size() - missN);
      for (int i = 0; i < hetN; i += 1) {
        b.set(missN + i, "0/1");
      }
    });
  }

  private static Scenario combinations(String gene, int hetCount) {
    return new Scenario("combinations/" + gene + "/het=" + hetCount, gene, true, false, true, b -> {
      int n = Math.min(hetCount, b.size());
      for (int i = 0; i < n; i += 1) {
        b.set(i, "0/1");
      }
    });
  }

  private static Scenario phaseSets(String gene, int hetCount) {
    return new Scenario("phase-sets/" + gene + "/het=" + hetCount, gene, false, true, true, b -> {
      int n = Math.min(hetCount, b.size());
      // split hets across two phase sets
      int half = Math.max(1, n / 2);
      for (int i = 0; i < half; i += 1) {
        b.set(i, "0|1", 1);
      }
      for (int i = half; i < n; i += 1) {
        b.set(i, "0|1", 2);
      }
    });
  }


  // --- DPYD HapB3 scenarios ---

  private static Scenario dpydHapB3Only(boolean phased) {
    String tag = phased ? "phased" : "unphased";
    String gt = phased ? "0|1" : "0/1";
    return new Scenario("dpyd-hapb3-only/" + tag, "DPYD", false, true, true, b -> {
      if (phased) {
        b.setDefaultPhased(true);
      }
      b.setRsid(DpydHapB3Matcher.HAPB3_EXONIC_RSID, gt);
      b.setRsid(DpydHapB3Matcher.HAPB3_INTRONIC_RSID, gt);
    });
  }

  private static Scenario dpydHapB3WithOthers(boolean phased) {
    String tag = phased ? "phased" : "unphased";
    String gt = phased ? "0|1" : "0/1";
    return new Scenario("dpyd-hapb3+others/" + tag, "DPYD", false, true, true, b -> {
      if (phased) {
        b.setDefaultPhased(true);
      }
      b.setRsid(DpydHapB3Matcher.HAPB3_EXONIC_RSID, gt);
      b.setRsid(DpydHapB3Matcher.HAPB3_INTRONIC_RSID, gt);
      // add 3 other DPYD hets, skipping HapB3 positions
      int added = 0;
      for (int i = 0; i < b.size() && added < 3; i += 1) {
        String rsid = b.entry(i).rsid();
        if (DpydHapB3Matcher.HAPB3_EXONIC_RSID.equals(rsid) ||
            DpydHapB3Matcher.HAPB3_INTRONIC_RSID.equals(rsid)) {
          continue;
        }
        b.set(i, gt);
        added += 1;
      }
    });
  }

  private static Scenario dpydHapB3IntronExonMismatch() {
    // exon het, intron ref: partial HapB3 evidence
    return new Scenario("dpyd-hapb3-intron-exon-mismatch", "DPYD", false, true, true,
        b -> b.setRsid(DpydHapB3Matcher.HAPB3_EXONIC_RSID, "0/1"));
  }
}
