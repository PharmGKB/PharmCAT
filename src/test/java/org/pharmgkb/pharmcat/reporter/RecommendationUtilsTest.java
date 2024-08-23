package org.pharmgkb.pharmcat.reporter;

import java.util.HashMap;
import java.util.Map;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.Assertions.assertTrue;


class RecommendationUtilsTest {

  @Test
  void testMapContains() {
    Map<String, Object> mapA = new HashMap<>();
    Map<String, Object> mapB = new HashMap<>();

    // empty maps don't contain each other
    assertFalse(RecommendationUtils.mapContains(mapA, mapB));

    // empty map cannot be contained by another map
    mapA.put("key1", "value1");
    assertFalse(RecommendationUtils.mapContains(mapA, mapB));

    // same exact values should contain
    mapB.put("key1", "value1");
    assertTrue(RecommendationUtils.mapContains(mapA, mapB));

    // adding an extra value to the parent should contain
    mapA.put("key2", "value2");
    assertTrue(RecommendationUtils.mapContains(mapA, mapB));

    // changing the value for key that exists in both maps should break containment
    mapA.put("key1", "valueChanged");
    assertFalse(RecommendationUtils.mapContains(mapA, mapB));

    Map<String, Object> mapC = new HashMap<>();
    Map<String, Object> mapD = new HashMap<>();
    mapC.put("key1", Map.of("a1", 2));
    mapC.put("key2", Map.of("a2", 2));
    mapD.put("key1", Map.of("a1", 2));

    assertTrue(RecommendationUtils.mapContains(mapC, mapD));
  }
}
