package org.cpic.haplotype;

import static org.junit.Assert.assertEquals;

import java.io.File;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.List;

import org.junit.Test;

import com.google.common.collect.ListMultimap;

public class DefinitionReaderTest {

	@Test
	public void testReader() throws Exception {
		
		System.out.println("DefinitionReaderTest");
		
		DefinitionReader dr = new DefinitionReader();
		//ClassLoader classLoader = getClass().getClassLoader();
        File file = new File(DefinitionReader.class.getResource("SLCO1B1.tsv").getFile());
        Path path = Paths.get(file.getAbsolutePath());
		dr.read(path);
		
		ListMultimap<String, Variant> m_haplotypePositions=dr.getHaplotypePositions();
		ListMultimap<String, Haplotype> m_haplotypes=dr.getHaplotypes();
		
		System.out.println(m_haplotypes);
		System.out.println(m_haplotypePositions);
		
		List<Variant> v_list = m_haplotypePositions.get("SLCO1B1");
		List<Haplotype> h_list = m_haplotypes.get("SLCO1B1");

		System.out.println(v_list.size());
		System.out.println(h_list.size());
		
		assertEquals(29,v_list.size());
		assertEquals(37,h_list.size());
		
		for (Variant v : v_list) {
			//System.out.println(v.getREF());
		    //System.out.println(v.getALTs());
		}

		for (Haplotype h : h_list) {
			//System.out.println(v.getREF());
		    //System.out.println(v.getALTs());
		}
		
		
	
		
				
		
		
	}
	
	
}
