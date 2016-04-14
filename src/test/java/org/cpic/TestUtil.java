package org.cpic;

import java.net.URISyntaxException;
import java.net.URL;
import java.nio.file.Path;
import java.nio.file.Paths;

import static org.junit.Assert.assertNotNull;


/**
 * @author Mark Woon
 */
public class TestUtil {

  public static Path getFile(String file) throws URISyntaxException {

    URL url = TestUtil.class.getClassLoader().getResource(file);
    assertNotNull(url);
    return Paths.get(url.toURI());
  }
}
