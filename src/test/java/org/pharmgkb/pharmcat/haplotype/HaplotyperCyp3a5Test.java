package org.pharmgkb.pharmcat.haplotype;

import java.nio.file.Path;
import java.util.List;
import java.util.SortedMap;
import com.google.common.collect.Lists;
import org.junit.Before;
import org.junit.Test;
import org.pharmgkb.pharmcat.TestUtil;
import org.pharmgkb.pharmcat.haplotype.model.DiplotypeMatch;


/**
 * JUnit test for {@link Haplotyper#callDiplotypes(SortedMap, String)}.
 *
 * @author Lester Carter
 */
public class HaplotyperCyp3a5Test {
  private Path m_tsvFile;

  @Before
  public void before() throws Exception {
    m_tsvFile =  TestUtil.getFile("org/pharmgkb/pharmcat/haplotype/CYP3A5.tsv");
  }


  /*
  TODO Lester - CYP3A5 tsv file contains a range. Does input vcf need just first position or all?  Check with Mark
  This example shows two substantial issues that will effect all the more complex examples

   1) How are insertions and deletions handled
   2) How are characters like Y handled that need some translation
   3) '-' in tsv file for star 1

  Also repeated code above
   */

  @Test
  public void cyp3a5s3d9() throws Exception {
    // Test *3/*9.  Note Y in 99672916 position
    /*

    TODO: Fails due to '-' in tsv file: Allele	CA100125.1	*1	Normal function	G	T	A	T	C	C	-	G
    Part of 'how ranges are dealt with' issue.  If a range should the vcf contain all positions?  And should *1 actually
    Have the reference here?  Or even better blank for wild card?

     */
    Path vcfFile = TestUtil.getFile("org/pharmgkb/pharmcat/haplotype/cyp3a5/s3s9.vcf");
    List<DiplotypeMatch> matches = HaplotyperTest.testCallHaplotype(m_tsvFile, vcfFile);

    List<String> expectedMatches = Lists.newArrayList("*3/*9");
    TestUtil.assertDiplotypePairs(expectedMatches, matches);
  }

  @Test
  public void cyp3a5s1s7() throws Exception {
    // Test *1/*7.  Het in rs41303343 position, a single insertion. Will currently fail to '-' issue in tsv file

    /*

     This comes down to what we expect in the ref and alt position.  In the 1000G data it is encoded lke this:
      grep -i "rs41303343" CYP3A5.recode.hg38.vcf
      99652770 rs41303343 T TA 100 PASS AA=A|A|AA|insertion...

      Again I expect this to fail as T and TA are not what it is in the tsv file, though this is what I've placed in the
      vcf for now.

     */

    Path vcfFile = TestUtil.getFile("org/pharmgkb/pharmcat/haplotype/cyp3a5/s1s7.vcf");
    List<DiplotypeMatch> matches = HaplotyperTest.testCallHaplotype(m_tsvFile, vcfFile);

    List<String> expectedMatches = Lists.newArrayList("*1/*7");
    TestUtil.assertDiplotypePairs(expectedMatches, matches);
  }


}
