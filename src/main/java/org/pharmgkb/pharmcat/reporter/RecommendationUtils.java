package org.pharmgkb.pharmcat.reporter;

import java.util.Map;


public class RecommendationUtils {

  public static boolean mapContains(Map<String, Object> superSet, Map<String, Object> subSet) {
    return superSet != null && superSet.size() > 0 && subSet != null && subSet.size() > 0 &&
        superSet.size() >= subSet.size() &&
        superSet.entrySet().containsAll(subSet.entrySet());
  }
}
