package org.pharmgkb.pharmcat;

import java.net.URISyntaxException;
import java.net.URL;
import java.nio.file.Path;
import java.nio.file.Paths;
import javax.annotation.Nonnull;
import com.google.common.base.Preconditions;

import static org.junit.Assert.assertNotNull;


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
    assertNotNull(filename + " is null", url);
    return Paths.get(url.toURI());
  }
}
