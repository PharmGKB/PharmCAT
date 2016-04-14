package org.cpic.haplotype;

import java.io.File;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.List;

import org.junit.Test;

import com.google.common.collect.ListMultimap;

public class DefinitionReaderTest {

	@Test
	public void testReader() throws Exception {
		
		System.out.println("DefinitionReaderTest");
		
		DefinitionReader dr = new DefinitionReader();
		//ClassLoader classLoader = getClass().getClassLoader();
        File file = new File(DefinitionReader.class.getResource("template.allele.translation.tsv").getFile());
        Path path = Paths.get(file.getAbsolutePath());
		dr.read(path);
		ListMultimap<String, Variant> m_haplotypePositions=dr.getHaplotypePositions();
		List<Variant> v_list = m_haplotypePositions.get("ABCD");
		System.out.println(v_list.size());
		for (Variant v : v_list) {
		    System.out.println(v.getREF());
		}
		System.out.println(m_haplotypePositions);
		
	
		
				
		
		
	}
	
	
}
