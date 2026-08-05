package org.pharmgkb.pharmcat.util;

import org.junit.jupiter.api.Test;
import org.pharmgkb.pharmcat.definition.model.DefinitionFile;

import static org.junit.jupiter.api.Assertions.assertNotNull;
import static org.junit.jupiter.api.Assertions.assertThrows;


/**
 * JUnit test for {@link DataManager}.
 *
 * @author Mark Woon
 */
class DataManagerTest {
  /**
   * Minimal synthetic SLCO1B1 definition file with the 3 core positions {@code fixSlco1b1()} expects
   * (N130D/rs2306283, V174A/rs4149056, R580X/rs71581941). *1 is reference; *15 and *45 differ from *1 at every
   * position {@code fixSlco1b1()} validates.
   */
  private static final String VALID_JSON = """
      {
        "gene": "SLCO1B1",
        "variants": [
          {"chromosome": "chr1", "position": 100, "rsid": "rs2306283"},
          {"chromosome": "chr1", "position": 200, "rsid": "rs4149056"},
          {"chromosome": "chr1", "position": 300, "rsid": "rs71581941"}
        ],
        "positionToLocus": {"100": 0, "200": 1, "300": 2},
        "namedAlleles": [
          {"id": "SLCO1B1.001", "name": "*1", "alleles": ["A", "C", "G"], "cpicAlleles": ["A", "C", "G"], "reference": true},
          {"id": "SLCO1B1.015", "name": "*15", "alleles": ["G", "T", "C"], "cpicAlleles": ["G", "T", "C"], "reference": false},
          {"id": "SLCO1B1.045", "name": "*45", "alleles": ["G", "T", "A"], "cpicAlleles": ["G", "T", "A"], "reference": false}
        ]
      }
      """;

  /**
   * Same as {@link #VALID_JSON}, except *45's R580X (3rd position) allele is identical to *1's - a data regression
   * that {@code fixSlco1b1()} is supposed to catch.
   */
  private static final String UNCHANGED_R580X_JSON = """
      {
        "gene": "SLCO1B1",
        "variants": [
          {"chromosome": "chr1", "position": 100, "rsid": "rs2306283"},
          {"chromosome": "chr1", "position": 200, "rsid": "rs4149056"},
          {"chromosome": "chr1", "position": 300, "rsid": "rs71581941"}
        ],
        "positionToLocus": {"100": 0, "200": 1, "300": 2},
        "namedAlleles": [
          {"id": "SLCO1B1.001", "name": "*1", "alleles": ["A", "C", "G"], "cpicAlleles": ["A", "C", "G"], "reference": true},
          {"id": "SLCO1B1.015", "name": "*15", "alleles": ["G", "T", "C"], "cpicAlleles": ["G", "T", "C"], "reference": false},
          {"id": "SLCO1B1.045", "name": "*45", "alleles": ["G", "T", "G"], "cpicAlleles": ["G", "T", "G"], "reference": false}
        ]
      }
      """;


  private DefinitionFile loadAndInitialize(String json) {
    DefinitionFile definitionFile = DataSerializer.GSON.fromJson(json, DefinitionFile.class);
    for (var namedAllele : definitionFile.getNamedAlleles()) {
      namedAllele.initialize(definitionFile.getVariants());
    }
    return definitionFile;
  }


  @Test
  void fixSlco1b1CreatesSuballelesForValidData() {
    DefinitionFile definitionFile = loadAndInitialize(VALID_JSON);

    DataManager.fixSlco1b1(definitionFile);

    assertNotNull(definitionFile.getNamedAllele("*45.001"));
    assertNotNull(definitionFile.getNamedAllele("*45.002"));
  }

  /**
   * Regression test: {@code fixSlco1b1()} must compare *45's R580X allele (the 3rd core position) against *1's
   * R580X allele, not against *1's V174A allele (the 2nd position). With the bug, this data regression (*45's
   * R580X identical to *1's) went undetected.
   */
  @Test
  void fixSlco1b1DetectsUnchangedR580X() {
    DefinitionFile definitionFile = loadAndInitialize(UNCHANGED_R580X_JSON);

    IllegalStateException ex = assertThrows(IllegalStateException.class,
        () -> DataManager.fixSlco1b1(definitionFile));
    assertNotNull(ex.getMessage());
  }
}
