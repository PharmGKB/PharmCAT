package org.pharmgkb.pharmcat.reporter.model.pgkb;

import java.util.HashMap;
import java.util.List;
import java.util.Map;
import org.junit.jupiter.api.Test;
import org.pharmgkb.pharmcat.reporter.model.result.Diplotype;
import org.pharmgkb.pharmcat.reporter.model.result.Genotype;

import static org.junit.jupiter.api.Assertions.assertTrue;


class RecommendationAnnotationTest {

  @Test
  void testMatchesGenotype() {
    Map<String,Object> sourceGenotype = new HashMap<>();
    sourceGenotype.put("GENEX", "TEST");

    RecommendationAnnotation recommendation = new RecommendationAnnotation();
    recommendation.setLookupKey(sourceGenotype);

    Diplotype diplotype = new Diplotype("GENEX", "TEST");
    Diplotype diplotype2 = new Diplotype("GENEY", "TOAST");
    Genotype genotype = Genotype.forTest(List.of(diplotype, diplotype2));
    assertTrue(recommendation.matchesGenotype(genotype));
  }
}
