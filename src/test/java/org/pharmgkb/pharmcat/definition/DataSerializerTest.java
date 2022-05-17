package org.pharmgkb.pharmcat.definition;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.Optional;
import java.util.Set;
import com.google.common.base.Charsets;
import org.apache.commons.io.FileUtils;
import org.junit.jupiter.api.Test;
import org.pharmgkb.common.util.PathUtils;
import org.pharmgkb.pharmcat.TestUtils;
import org.pharmgkb.pharmcat.definition.model.DefinitionExemption;
import org.pharmgkb.pharmcat.util.DataSerializer;

import static org.junit.jupiter.api.Assertions.*;


/**
 * JUnit test for {@link DataSerializer}.
 *
 * @author Mark Woon
 */
class DataSerializerTest {

  @Test
  void testExemptions() throws Exception {

    Path tsvFile = PathUtils.getPathToResource("org/pharmgkb/pharmcat/definition/exemptions.tsv");
    Path jsonFile = TestUtils.createTempFile("transformExemptions", ".json");
    Path refJsonFile = PathUtils.getPathToResource("org/pharmgkb/pharmcat/definition/exemptions.json");

    try {
      DataSerializer dataSerializer = new DataSerializer();
      Set<DefinitionExemption> definitionExemptions = dataSerializer.deserializeExemptionsFromTsv(tsvFile);
      assertEquals(3, definitionExemptions.size());

      Optional<DefinitionExemption> cyp2c9 = definitionExemptions.stream()
          .filter((de) -> de.getGene().equals("CYP2C9"))
          .findFirst();
      assertTrue(cyp2c9.isPresent());
      assertEquals(0, cyp2c9.get().getIgnoredPositions().size());
      assertEquals(1, cyp2c9.get().getExtraPositions().size());
      assertEquals(0, cyp2c9.get().getIgnoredAlleles().size());
      assertNull(cyp2c9.get().isAllHits());
      assertNull(cyp2c9.get().isAssumeReference());

      Optional<DefinitionExemption> cyp2c19 = definitionExemptions.stream()
          .filter((de) -> de.getGene().equals("CYP2C19"))
          .findFirst();
      assertTrue(cyp2c19.isPresent());
      assertEquals(1, cyp2c19.get().getIgnoredPositions().size());
      assertEquals(0, cyp2c19.get().getExtraPositions().size());
      assertEquals(0, cyp2c19.get().getIgnoredAlleles().size());
      assertNull(cyp2c19.get().isAllHits());
      assertNull(cyp2c19.get().isAssumeReference());

      // write it out
      dataSerializer.serializeToJson(definitionExemptions, jsonFile);
      // compare with expected
      FileUtils.contentEqualsIgnoreEOL(refJsonFile.toFile(), jsonFile.toFile(), Charsets.UTF_8.displayName());

    } catch (IOException ex) {
      if (!TestUtils.isContinuousIntegration() || TestUtils.isIgnorableExternalServiceException(ex.getCause())) {
        throw ex;
      }
    } finally {
      Files.deleteIfExists(jsonFile);
    }
  }

  @Test
  void testValuedBooleanExemptions() throws Exception {

    Path tsvFile = PathUtils.getPathToResource("org/pharmgkb/pharmcat/definition/exemptions.tsv");
    Path jsonFile = TestUtils.createTempFile("transformExemptions", ".json");
    Path refJsonFile = PathUtils.getPathToResource("org/pharmgkb/pharmcat/definition/exemptions.json");

    try {
      DataSerializer dataSerializer = new DataSerializer();
      Set<DefinitionExemption> definitionExemptions = dataSerializer.deserializeExemptionsFromTsv(tsvFile);
      assertEquals(3, definitionExemptions.size());

      Optional<DefinitionExemption> cyp2c9 = definitionExemptions.stream()
          .filter((de) -> de.getGene().equals("CYP2C9"))
          .findFirst();
      assertTrue(cyp2c9.isPresent());
      assertEquals(0, cyp2c9.get().getIgnoredPositions().size());
      assertEquals(1, cyp2c9.get().getExtraPositions().size());
      assertEquals(0, cyp2c9.get().getIgnoredAlleles().size());
      assertNull(cyp2c9.get().isAllHits());
      assertNull(cyp2c9.get().isAssumeReference());

      Optional<DefinitionExemption> cyp2c19 = definitionExemptions.stream()
          .filter((de) -> de.getGene().equals("CYP2C19"))
          .findFirst();
      assertTrue(cyp2c19.isPresent());
      assertEquals(1, cyp2c19.get().getIgnoredPositions().size());
      assertEquals(0, cyp2c19.get().getExtraPositions().size());
      assertEquals(0, cyp2c19.get().getIgnoredAlleles().size());
      assertNull(cyp2c19.get().isAllHits());
      assertNull(cyp2c19.get().isAssumeReference());

      Optional<DefinitionExemption> cyp3a5 = definitionExemptions.stream()
          .filter((de) -> de.getGene().equals("CYP3A5"))
          .findFirst();
      assertTrue(cyp3a5.isPresent());
      assertEquals(0, cyp3a5.get().getIgnoredPositions().size());
      assertEquals(0, cyp3a5.get().getExtraPositions().size());
      assertEquals(0, cyp3a5.get().getIgnoredAlleles().size());
      assertEquals(Boolean.TRUE, cyp3a5.get().isAllHits());
      assertFalse(cyp3a5.get().isAssumeReference());

      // write it out
      dataSerializer.serializeToJson(definitionExemptions, jsonFile);
      // compare with expected
      FileUtils.contentEqualsIgnoreEOL(refJsonFile.toFile(), jsonFile.toFile(), Charsets.UTF_8.displayName());

    } catch (IOException ex) {
      if (!TestUtils.isContinuousIntegration() || TestUtils.isIgnorableExternalServiceException(ex.getCause())) {
        throw ex;
      }
    } finally {
      Files.deleteIfExists(jsonFile);
    }
  }
}
