package org.pharmgkb.pharmcat.util;

import java.io.IOException;
import java.nio.file.Path;
import java.util.Map;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.TestInfo;
import org.pharmgkb.common.util.PathUtils;
import org.pharmgkb.pharmcat.TestUtils;
import org.pharmgkb.pharmcat.definition.model.DefinitionExemption;

import static org.junit.jupiter.api.Assertions.*;


/**
 * JUnit test for {@link DataSerializer}.
 *
 * @author Mark Woon
 */
class DataSerializerTest {

  @Test
  void exemptionsOnly(TestInfo testInfo) throws Exception {

    Path tsvFile = PathUtils.getPathToResource("org/pharmgkb/pharmcat/definition/DataSerializerTest-exemptions.tsv");
    Path jsonFile = TestUtils.createTestFile(testInfo, ".json");
    Path refJsonFile = PathUtils.getPathToResource("org/pharmgkb/pharmcat/definition/DataSerializerTest-exemptions.json");

    try {
      DataSerializer dataSerializer = new DataSerializer();
      Map<String, DefinitionExemption> definitionExemptions = dataSerializer.deserializeExemptionsFromTsv(tsvFile, null);
      checkExemptions(definitionExemptions, false);

      // write it out
      DataSerializer.serializeToJson(definitionExemptions, jsonFile);
      // compare with expected
      TestUtils.assertEqual(refJsonFile, jsonFile);

      // read it back in
      definitionExemptions = dataSerializer.deserializeExemptionsFromJson(jsonFile);
      checkExemptions(definitionExemptions, false);

    } catch (IOException ex) {
      if (!TestUtils.isContinuousIntegration() || TestUtils.isIgnorableExternalServiceException(ex.getCause())) {
        throw ex;
      }
    }
  }

  @Test
  void exemptionsWithUnphasedPriorities(TestInfo testInfo) throws Exception {

    Path exemptionsTsvFile = PathUtils.getPathToResource("org/pharmgkb/pharmcat/definition/DataSerializerTest-exemptions.tsv");
    Path unphasedPrioritiesTsvFile = PathUtils.getPathToResource("org/pharmgkb/pharmcat/definition/DataSerializerTest-unphasedPriorities.tsv");
    Path jsonFile = TestUtils.createTestFile(testInfo, ".json");
    Path refJsonFile = PathUtils.getPathToResource("org/pharmgkb/pharmcat/definition/DataSerializerTest-exemptionsWithUnphasedPriorities.json");

    try {
      DataSerializer dataSerializer = new DataSerializer();
      Map<String, DefinitionExemption> definitionExemptions =
          dataSerializer.deserializeExemptionsFromTsv(exemptionsTsvFile, unphasedPrioritiesTsvFile);
      checkExemptions(definitionExemptions, true);

      // write it out
      DataSerializer.serializeToJson(definitionExemptions, jsonFile);
      // compare with expected
      TestUtils.assertEqual(refJsonFile, jsonFile);

      // read it back in
      definitionExemptions = dataSerializer.deserializeExemptionsFromJson(jsonFile);
      checkExemptions(definitionExemptions, true);

    } catch (IOException ex) {
      if (!TestUtils.isContinuousIntegration() || TestUtils.isIgnorableExternalServiceException(ex.getCause())) {
        throw ex;
      }
    }
  }


  private void checkExemptions(Map<String, DefinitionExemption> definitionExemptions, boolean withPriorities) {
    assertEquals(4, definitionExemptions.size());

    assertNull(definitionExemptions.get("TPMT"));

    DefinitionExemption cyp2c9 = definitionExemptions.get("CYP2C9");
    assertNotNull(cyp2c9);
    assertEquals(0, cyp2c9.getRequiredPositions().size());
    assertEquals(0, cyp2c9.getIgnoredPositions().size());
    assertEquals(1, cyp2c9.getExtraPositions().size());
    assertEquals(0, cyp2c9.getIgnoredAlleles().size());

    DefinitionExemption g6pd = definitionExemptions.get("G6PD");
    assertNotNull(g6pd);
    assertEquals(0, g6pd.getRequiredPositions().size());
    assertEquals(2, g6pd.getIgnoredPositions().size());
    assertEquals(0, g6pd.getExtraPositions().size());
    assertEquals(1, g6pd.getIgnoredAlleles().size());

    DefinitionExemption nat2 = definitionExemptions.get("NAT2");
    assertNotNull(nat2);
    assertEquals(4, nat2.getRequiredPositions().size());
    assertEquals(0, nat2.getIgnoredPositions().size());
    assertEquals(0, nat2.getExtraPositions().size());
    assertEquals(0, nat2.getIgnoredAlleles().size());
    if (withPriorities) {
      assertEquals(4, nat2.getUnphasedDiplotypePriorities().size());
    }
  }
}
