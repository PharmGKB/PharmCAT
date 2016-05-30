package org.pharmgkb.pharmcat.haplotype;

import java.nio.file.Path;
import java.util.List;
import com.google.common.collect.Lists;
import org.junit.Before;
import org.junit.Test;
import org.pharmgkb.common.util.PathUtils;
import org.pharmgkb.pharmcat.haplotype.model.HaplotyperResult;

import static org.pharmgkb.pharmcat.haplotype.HaplotyperTest.assertDiplotypePairs;
import static org.pharmgkb.pharmcat.haplotype.HaplotyperTest.testCallHaplotype;


/**
 * JUnit test for {@link Haplotyper#callDiplotypes(MatchData)}.
 *
 * @author Lester Carter
 */
public class RealDataTest {
  private Path m_definitionFile;

  @Before
  public void before() throws Exception {
    m_definitionFile =  PathUtils.getPathToResource("org/pharmgkb/pharmcat/haplotype/CYP3A5_translation.json");
  }

  @Test
  public void cyp3a5s3s3() throws Exception {
    /*
    Call is *3/*3 from dmet data, but we get *3/*9.  *3/*3 should be an option based on lab discussion
    I have checked in  a very small version of this vcf file, that only contains the positions found in the
    supplied data.
    */
    Path vcfFile = PathUtils.getPathToResource("org/pharmgkb/pharmcat/haplotype/realtestdata/wgstestcyp3a5.vcf");
    List<String> expectedMatches = Lists.newArrayList("*3/*3");

    HaplotyperResult result = testCallHaplotype(m_definitionFile, vcfFile, true, false, true);
    assertDiplotypePairs(expectedMatches, result);
  }

}
