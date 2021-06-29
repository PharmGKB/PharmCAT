package org.pharmgkb.pharmcat.reporter;

import java.io.IOException;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.assertEquals;


class DrugCollectionTest {
  @Test
  void testLoad() throws IOException {
    DrugCollection drugCollection = new DrugCollection();
    assertEquals(66, drugCollection.size());
  }
}
