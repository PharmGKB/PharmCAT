package org.pharmgkb.pharmcat.definition;

import java.io.IOException;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.assertEquals;

/**
 * JUnit test for {@link AssemblyMap}.
 *
 * @author Ryan Whaley
 */
class AssemblyMapTest {

  @Test
  void testGet() throws IOException {
    AssemblyMap assemblyMap = new AssemblyMap();

    assertEquals(AssemblyMap.GRCH37, assemblyMap.get("NC_000010.10"));
    assertEquals(AssemblyMap.GRCH38, assemblyMap.get("NC_000010.11"));

    assertEquals(AssemblyMap.GRCH37, assemblyMap.get("NC_000004.11"));
    assertEquals(AssemblyMap.GRCH38, assemblyMap.get("NC_000004.12"));
  }
}
