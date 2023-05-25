package org.pharmgkb.pharmcat;

import java.util.Arrays;
import java.util.Collection;
import java.util.Set;
import java.util.TreeSet;
import java.util.function.BinaryOperator;
import java.util.function.Function;
import java.util.stream.Collectors;
import org.checkerframework.checker.nullness.qual.Nullable;
import org.pharmgkb.pharmcat.reporter.TextConstants;
import org.pharmgkb.pharmcat.util.HaplotypeNameComparator;


/**
 * Helper functions for working with diplotypes in tests.
 *
 * @author Mark Woon
 */
public class DiplotypeUtils {

  /**
   * This Function can be used in reduce() calls
   */
  public static final BinaryOperator<String> PhasedReducer = (a,b)-> {
    if (a == null) {
      return b;
    }
    if (b == null) {
      return a;
    }

    int left = 0;
    int right = 1;

    String[] aHaps = a.split(TextConstants.GENOTYPE_DELIMITER);
    String[] bHaps = b.split(TextConstants.GENOTYPE_DELIMITER);

    Set<String> finalLeft = new TreeSet<>(HaplotypeNameComparator.getComparator());
    Set<String> finalRight = new TreeSet<>(HaplotypeNameComparator.getComparator());

    addHaps(finalLeft, aHaps[left]);
    addHaps(finalRight, aHaps[right]);

    if (finalLeft.contains(bHaps[right]) || finalRight.contains(bHaps[left])) {
      addHaps(finalLeft, bHaps[right]);
      addHaps(finalRight, bHaps[left]);
    }
    else {
      addHaps(finalLeft, bHaps[left]);
      addHaps(finalRight, bHaps[right]);
    }

    String joined0 = String.join("+", finalLeft);
    String joined1 = String.join("+", finalRight);

    return joined0+ TextConstants.GENOTYPE_DELIMITER +joined1;
  };

  /**
   * Takes diplotype names and encloses them in parentheses if they contain a "+"
   */
  private static final Function<String,String> hapNameEncloser = h -> {
    if (h.contains("+")) return "(" + h + ")"; else return h;
  };

  private static void addHaps(Set<String> hapSet, String hap) {
    hapSet.addAll(Arrays.asList(hap.split("\\+")));
  }


  /**
   * When multiple diplotype calls are phased they can be combined into one String in certain circumstances. This method
   * combines things like a/b and a/c into a/b+c.
   * @param dips a {@link Collection} of diplotype strings in the form "a/b", "a/c", etc...
   * @return a single String combining the diplotypes, e.g. "a/b+c"
   */
  public static String reducePhasedDiplotypes(@Nullable Collection<String> dips) {
    if (dips == null || dips.size() == 0) return "";

    final Set<String> leftBucket = new TreeSet<>(HaplotypeNameComparator.getComparator());
    final Set<String> rightBucket = new TreeSet<>(HaplotypeNameComparator.getComparator());

    dips.forEach(d -> {
      String[] haps = d.split(TextConstants.GENOTYPE_DELIMITER);
      if (!rightBucket.contains(haps[0])) {
        leftBucket.add(haps[0]);
        rightBucket.add(haps[1]);
      } else {
        leftBucket.add(haps[1]);
        rightBucket.add(haps[0]);
      }
    });

    Function<String,String> leftMapper = leftBucket.size() > 1 ? hapNameEncloser : Function.identity();
    String left = leftBucket.stream().map(leftMapper).collect(Collectors.joining("+"));
    Function<String,String> rightMapper = rightBucket.size() > 1 ? hapNameEncloser : Function.identity();
    String right = rightBucket.stream().map(rightMapper).collect(Collectors.joining("+"));

    return left + TextConstants.GENOTYPE_DELIMITER + right;
  }
}
