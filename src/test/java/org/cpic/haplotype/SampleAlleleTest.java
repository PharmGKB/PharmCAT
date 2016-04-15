package org.cpic.haplotype;

import com.google.common.collect.Sets;
import org.junit.Test;

import static org.junit.Assert.assertTrue;


/**
 * JUnit test for {@link SampleAllele}.
 *
 * @author Mark Woon
 */
public class SampleAlleleTest {


  @Test
  public void testCompare() {

    SampleAllele sa1 = new SampleAllele("chr1", 1, "A", "A", true);
    SampleAllele sa2 = new SampleAllele("chr1", 2, "A", "A", true);
    SampleAllele sa3 = new SampleAllele("chr2", 1, "A", "A", true);
    SampleAllele sa4 = new SampleAllele("chrX", 1, "A", "A", true);
    SampleAllele sa5 = new SampleAllele("chrY", 1, "A", "A", true);
    SampleAllele sa6 = new SampleAllele("chr11", 1, "A", "A", true);

    assertTrue(Sets.newTreeSet(Sets.newHashSet(sa1, sa2)).first() == sa1);
    assertTrue(Sets.newTreeSet(Sets.newHashSet(sa2, sa1)).first() == sa1);
    assertTrue(Sets.newTreeSet(Sets.newHashSet(sa2, sa3)).first() == sa2);
    assertTrue(Sets.newTreeSet(Sets.newHashSet(sa6, sa2)).first() == sa2);
    assertTrue(Sets.newTreeSet(Sets.newHashSet(sa1, sa4)).first() == sa1);
    assertTrue(Sets.newTreeSet(Sets.newHashSet(sa4, sa5)).first() == sa4);
    assertTrue(Sets.newTreeSet(Sets.newHashSet(sa5, sa4)).first() == sa4);
  }

}
