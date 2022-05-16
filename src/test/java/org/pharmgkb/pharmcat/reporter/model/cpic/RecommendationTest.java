package org.pharmgkb.pharmcat.reporter.model.cpic;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import org.junit.jupiter.api.Test;
import org.pharmgkb.pharmcat.reporter.model.result.Diplotype;
import org.pharmgkb.pharmcat.reporter.model.result.Genotype;

import static org.junit.jupiter.api.Assertions.assertTrue;


class RecommendationTest {

  @Test
  void testMatchesLookupKey() {
    List<Map<String,String>> genotypes = new ArrayList<>();
    Map<String,String> genotype = new HashMap<>();
    genotype.put("GENEX", "TEST");
    genotypes.add(genotype);

    Recommendation recommendation = new Recommendation();
    recommendation.setGenotypes(genotypes);

    Map<String,String> testGenotype = new HashMap<>();
    testGenotype.put("GENEX", "TEST");
    assertTrue(recommendation.matchesLookupKey.test(testGenotype));
  }

  @Test
  void testMatchesGenotype() {
    List<Map<String,String>> genotypes = new ArrayList<>();
    Map<String,String> sourceGenotype = new HashMap<>();
    sourceGenotype.put("GENEX", "TEST");
    genotypes.add(sourceGenotype);

    Recommendation recommendation = new Recommendation();
    recommendation.setGenotypes(genotypes);

    Diplotype diplotype = new Diplotype("GENEX", "TEST");
    Genotype genotype = new Genotype(diplotype);
    Diplotype diplotype2 = new Diplotype("GENEY", "TOAST");
    genotype.addDiplotype(diplotype2);
    assertTrue(recommendation.matchesGenotype(genotype));
  }
}
