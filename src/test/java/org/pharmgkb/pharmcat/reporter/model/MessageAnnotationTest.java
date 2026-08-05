package org.pharmgkb.pharmcat.reporter.model;

import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.assertEquals;


/**
 * Test the {@link MessageAnnotation} object.
 *
 * @author Mark Woon
 */
class MessageAnnotationTest {

  /**
   * Regression test: {@code hashCode()} must be consistent with {@code equals()}. Two rows that differ only in
   * their match logic (gene/haplotypes/etc., not part of {@code equals()}) are still equal per {@code equals()}
   * (same name/version/exceptionType/message), so they must also produce equal hash codes.
   */
  @Test
  void testHashCodeConsistentWithEquals() {
    String rowA = "TestRule\tCYP2D6\t\t\t\t\t\t\t\tnote\tTest message\tv1";
    String rowB = "TestRule\tCYP2C19\t\t\t\t\t\t\t\tnote\tTest message\tv1";
    MessageAnnotation a = new MessageAnnotation(rowA);
    MessageAnnotation b = new MessageAnnotation(rowB);

    assertEquals(a, b, "should be equal - only differ in match logic, which equals() doesn't check");
    assertEquals(a.hashCode(), b.hashCode(), "equal objects must have equal hash codes");
  }
}
