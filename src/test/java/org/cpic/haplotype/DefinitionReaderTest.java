package org.cpic.haplotype;

import java.io.File;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.List;
import com.google.common.collect.ListMultimap;
import org.cpic.TestUtil;
import org.junit.Test;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

public class DefinitionReaderTest {

	@Test
	public void testReader() throws Exception {

		System.out.println("DefinitionReaderTest");

		DefinitionReader dr = new DefinitionReader();
		//ClassLoader classLoader = getClass().getClassLoader();
        File file = new File(DefinitionReader.class.getResource("CYP2C19.tsv").getFile());
        Path path = Paths.get(file.getAbsolutePath());
		dr.read(path);

		ListMultimap<String, Variant> m_haplotypePositions=dr.getHaplotypePositions();
		ListMultimap<String, Haplotype> m_haplotypes=dr.getHaplotypes();

		System.out.println(m_haplotypes);
		System.out.println(m_haplotypePositions);

		List<Variant> v_list = m_haplotypePositions.get("CYP2C19");
		List<Haplotype> h_list = m_haplotypes.get("CYP2C19");

		System.out.println(v_list.size());
		System.out.println(h_list.size());

		assertEquals(39,v_list.size());
		assertEquals(33,h_list.size());

		for (Variant v : v_list) {
			//System.out.println(v.getREF());
		    //System.out.println(v.getALTs());
		}

		for (Haplotype h : h_list) {
			//System.out.println(v.getREF());
		    //System.out.println(v.getALTs());
		}
	}


  @Test
  public void testReadAllDefinitions() throws Exception {

    Path file = TestUtil.getFile("org/cpic/haplotype/CYP2C19.tsv");
    DefinitionReader reader = new DefinitionReader();
    reader.read(file.getParent());

    for (String gene : reader.getHaplotypePositions().keySet()) {
      for (Variant variant : reader.getHaplotypePositions().get(gene)) {
        assertTrue(variant.getCHROM().startsWith("chr"));
        assertTrue(variant.getPOS() > 0);
      }
    }
  }
}
