package org.pharmgkb.pharmcat.reporter.model.cpic;

import java.util.HashMap;
import java.util.Map;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.Assertions.assertTrue;


class RecommendationTest {

  @Test
  void testMatchLookupKey() {
    Map<String, String> lookupKey = new HashMap<>();
    lookupKey.put("CYP2D6", "1");
    lookupKey.put("CYP2C19", "Poor Metabolizer");

    Recommendation recommendation = new Recommendation();
    recommendation.setLookupKey(lookupKey);

    // typical, well-formed value with sorted gene keys
    assertTrue(recommendation.matchLookupKey("CYP2C19:Poor Metabolizer;CYP2D6:1"));
    // same value but all lower-case should still pass
    assertTrue(recommendation.matchLookupKey("cyp2c19:poor metabolizer;cyp2d6:1"));
    // lookup key with unsorted gene keys will not work
    assertFalse(recommendation.matchLookupKey("CYP2D6:1;CYP2C19:Poor Metabolizer"));
  }
}
