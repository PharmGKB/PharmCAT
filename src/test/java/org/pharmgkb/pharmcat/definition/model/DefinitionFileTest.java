package org.pharmgkb.pharmcat.definition.model;

import java.util.List;
import org.junit.jupiter.api.Test;
import org.pharmgkb.pharmcat.definition.DefinitionReader;
import org.pharmgkb.pharmcat.util.DataManager;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertNotSame;
import static org.junit.jupiter.api.Assertions.assertSame;


/**
 * @author Mark Woon
 */
class DefinitionFileTest {

  /**
   * The shared {@link HaplotypeCandidateIndex} is cached per matching-flag key and must be discarded when the named
   * alleles are mutated, otherwise a mutation leaves a stale index built from the old definitions.
   */
  @Test
  void candidateIndexCachedAndInvalidatedOnMutation() throws Exception {
    DefinitionReader reader = new DefinitionReader(
        List.of(DataManager.getDefinitionFilePath("CYP3A5")),
        DataManager.DEFAULT_EXEMPTIONS_FILE);
    DefinitionFile df = reader.getDefinitionFile("CYP3A5");

    HaplotypeCandidateIndex first = df.getCandidateIndex(false, false);
    assertSame(first, df.getCandidateIndex(false, false), "same flags should return the cached index");

    NamedAllele toRemove = df.getNamedAlleles().stream()
        .filter(na -> !na.isReference())
        .findFirst()
        .orElseThrow();
    df.removeNamedAllele(toRemove);

    HaplotypeCandidateIndex afterRemove = df.getCandidateIndex(false, false);
    assertNotSame(first, afterRemove, "index should be rebuilt after a definition mutation");
    assertEquals(first.numHaplotypes() - 1, afterRemove.numHaplotypes(),
        "rebuilt index should reflect the removed allele");
  }
}
