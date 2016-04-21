package org.pharmgkb.pharmcat;

import java.net.URISyntaxException;
import java.net.URL;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.List;
import java.util.Set;
import java.util.TreeSet;
import java.util.stream.Collectors;
import javax.annotation.Nonnull;
import com.google.common.base.Preconditions;
import org.pharmgkb.pharmcat.haplotype.model.DiplotypeMatch;

import static org.junit.Assert.*;


/**
 * Utility methods for test cases.
 *
 * @author Mark Woon
 */
public class TestUtil {


  /**
   * Converts a file name into a {@link Path}.
   * We expect the file to be in the {@code test/resources} folder.
   *
   * @param filename a relative filename (e.g. {@code org/cpic/haplotype/CYP2C19.tsv})
   */
  public static @Nonnull Path getFile(@Nonnull String filename) throws URISyntaxException {

    Preconditions.checkNotNull(filename);
    URL url = TestUtil.class.getClassLoader().getResource(filename);
    assertNotNull(url);
    return Paths.get(url.toURI());
  }


  /**
   * Checks that the list of diplotype matches are what we expect.
   *
   * @param expectedPairs the set of expected diplotypes in "*1/*2" format
   * @param matches the list of haplotype matches
   */
  public static void assertDiplotypePairs(@Nonnull Set<String> expectedPairs,
      @Nonnull List<DiplotypeMatch> matches) {

    Preconditions.checkNotNull(expectedPairs);
    Preconditions.checkNotNull(matches);

    Set<String> pairs = matches.stream()
        .map(DiplotypeMatch::getName)
        .collect(Collectors.toCollection(TreeSet::new));
    assertEquals("Incoming matches has non-unique pairs", matches.size(), pairs.size());

    if (expectedPairs.size() != pairs.size() || !expectedPairs.equals(pairs)) {
      System.out.println("Expected: " + new TreeSet<>(expectedPairs));
      System.out.println("Got:      " + pairs);
      fail("Did not get expected matches");
    }
  }
}
