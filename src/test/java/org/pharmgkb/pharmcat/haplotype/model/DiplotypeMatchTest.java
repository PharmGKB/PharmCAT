package org.pharmgkb.pharmcat.haplotype.model;

import java.util.SortedSet;
import java.util.TreeMap;
import java.util.TreeSet;
import com.google.common.collect.Lists;
import org.junit.Test;
import org.pharmgkb.pharmcat.definition.model.NamedAllele;
import org.pharmgkb.pharmcat.definition.model.VariantLocus;
import org.pharmgkb.pharmcat.haplotype.MatchData;

import static org.junit.Assert.assertEquals;


/**
 * JUnit test for {@link DiplotypeMatch}.
 *
 * @author Mark Woon
 */
public class DiplotypeMatchTest {


  @Test
  public void testCompareTo() {

    VariantLocus var1 = new VariantLocus("chr1", 1, "g.1T>A");
    VariantLocus var2 = new VariantLocus("chr1", 2, "g.2T>A");
    VariantLocus var3 = new VariantLocus("chr1", 3, "g.3T>A");
    VariantLocus[] variants = new VariantLocus[] { var1, var2, var3 };


    NamedAllele hap1 = new NamedAllele("*1", "*1", new String[] { "T", "T", "T" });
    hap1.initialize(variants);
    NamedAllele hap2 = new NamedAllele("*4", "*4", new String[] { "A", "A", "A" });
    hap2.initialize(variants);
    NamedAllele hap3 = new NamedAllele("*3", "*3", new String[] { "T", "A", null });
    hap3.initialize(variants);

    HaplotypeMatch hm1 = new HaplotypeMatch(hap1);
    HaplotypeMatch hm2 = new HaplotypeMatch(hap2);
    HaplotypeMatch hm3 = new HaplotypeMatch(hap3);

    MatchData dataset = new MatchData(new TreeMap<>(), variants, null);

    DiplotypeMatch dm1 = new DiplotypeMatch(hm1, hm1, dataset);
    DiplotypeMatch dm2 = new DiplotypeMatch(hm1, hm2, dataset);
    DiplotypeMatch dm3 = new DiplotypeMatch(hm2, hm2, dataset);
    DiplotypeMatch dm4 = new DiplotypeMatch(hm3, hm2, dataset);

    SortedSet<DiplotypeMatch> matches = new TreeSet<>(Lists.newArrayList(dm1, dm2));
    assertEquals(dm1, matches.first());

    matches = new TreeSet<>(Lists.newArrayList(dm3, dm2));
    assertEquals(dm2, matches.first());

    matches = new TreeSet<>(Lists.newArrayList(dm3, dm4));
    assertEquals(dm3, matches.first());
  }
}
