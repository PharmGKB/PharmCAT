package org.pharmgkb.pharmcat.haplotype;

import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import org.junit.jupiter.api.Test;
import org.pharmgkb.pharmcat.definition.model.VariantLocus;
import org.pharmgkb.pharmcat.definition.model.VariantType;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertSame;


/**
 * JUnit test for {@link SampleAllele}.
 *
 * @author Mark Woon
 */
class SampleAlleleTest {


  @Test
  void testCompare() {

    SampleAllele sa1 = new SampleAllele("chr1", 1, "A", "A", true, Lists.newArrayList("A", "G"));
    SampleAllele sa2 = new SampleAllele("chr1", 2, "A", "A", true, Lists.newArrayList("A", "G"));
    SampleAllele sa3 = new SampleAllele("chr2", 1, "A", "A", true, Lists.newArrayList("A", "G"));
    SampleAllele sa4 = new SampleAllele("chrX", 1, "A", "A", true, Lists.newArrayList("A", "G"));
    SampleAllele sa5 = new SampleAllele("chrY", 1, "A", "A", true, Lists.newArrayList("A", "G"));
    SampleAllele sa6 = new SampleAllele("chr11", 1, "A", "A", true, Lists.newArrayList("A", "G"));

    assertSame(Sets.newTreeSet(Sets.newHashSet(sa1, sa2)).first(), sa1);
    assertSame(Sets.newTreeSet(Sets.newHashSet(sa2, sa1)).first(), sa1);
    assertSame(Sets.newTreeSet(Sets.newHashSet(sa2, sa3)).first(), sa2);
    assertSame(Sets.newTreeSet(Sets.newHashSet(sa6, sa2)).first(), sa2);
    assertSame(Sets.newTreeSet(Sets.newHashSet(sa1, sa4)).first(), sa1);
    assertSame(Sets.newTreeSet(Sets.newHashSet(sa4, sa5)).first(), sa4);
    assertSame(Sets.newTreeSet(Sets.newHashSet(sa5, sa4)).first(), sa4);
  }

  @Test
  void testForVariant() {

    VariantLocus insVariant = new VariantLocus("chr1", 1, "g.1A>AT");

    // inserts
    insVariant.setType(VariantType.INS);

    SampleAllele ins1 = new SampleAllele("chr1", 1, "TC", "TCA", true, Lists.newArrayList("TC", "TCA"));
    SampleAllele rez = ins1.forVariant(insVariant);
    assertEquals("del", rez.getAllele1());
    assertEquals("insA", rez.getAllele2());

    ins1 = new SampleAllele("chr1", 1, "TC", "TCATA", true, Lists.newArrayList("TC", "TCATA"));
    rez = ins1.forVariant(insVariant);
    assertEquals("del", rez.getAllele1());
    assertEquals("insATA", rez.getAllele2());

    // deletes
    insVariant.setType(VariantType.DEL);

    ins1 = new SampleAllele("chr1", 1, "TC", "T", true, Lists.newArrayList("TC", "T"));
    rez = ins1.forVariant(insVariant);
    assertEquals("C", rez.getAllele1());
    assertEquals("delC", rez.getAllele2());

    ins1 = new SampleAllele("chr1", 1, "TCAT", "T", true, Lists.newArrayList("TCAT", "T"));
    rez = ins1.forVariant(insVariant);
    assertEquals("CAT", rez.getAllele1());
    assertEquals("delCAT", rez.getAllele2());


    // repeats
    insVariant.setType(VariantType.REPEAT);
    insVariant.setReferenceRepeat("A(TA)6TAA");

    ins1 = new SampleAllele("chr1", 1, "ATATAA", "ATATATATATATAA", true, Lists.newArrayList("ATATAA", "ATATATATATATAA"));
    rez = ins1.forVariant(insVariant);
    assertEquals("A(TA)1TAA", rez.getAllele1());
    assertEquals("A(TA)5TAA", rez.getAllele2());

    ins1 = new SampleAllele("chr1", 1, "ATATATATATATATATAA", "ATATATATATATAA", true,
        Lists.newArrayList("ATATATATATATATATAA", "ATATATATATATAA"));
    rez = ins1.forVariant(insVariant);
    assertEquals("A(TA)7TAA", rez.getAllele1());
    assertEquals("A(TA)5TAA", rez.getAllele2());
  }
}
