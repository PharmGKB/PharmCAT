package org.cpic.haplotype;

import com.google.common.base.Preconditions;
import com.google.common.collect.Multimap;

import java.io.BufferedReader;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;

/**
 * Read in haplotype definition files.
 *
 * @author Mark Woon
 */
public class DefinitionReader {
  private Multimap<String, Variant> m_haplotypePositions;


  public Multimap<String, Variant> getHaplotypePositions() {
    return m_haplotypePositions;
  }



  public void read(Path path) throws IOException {

    if (Files.isDirectory(path)) {
      Files.list(path)
          .filter(f -> f.toString().endsWith(".tsv"))
          .forEach(this::readFile);
    } else {
      readFile(path);
    }
  }


  private void readFile(Path file)  {

    Preconditions.checkArgument(Files.isRegularFile(file));
    System.out.println(file);
    try (BufferedReader bufferedReader = Files.newBufferedReader(file)) {
    	
    	
    } catch (Exception ex) {
    	throw new RuntimeException(ex);
    }
  }


  public static void main(String[] args) {

    try {
      Path path = Paths.get(args[0]);
      DefinitionReader r = new DefinitionReader();
      r.read(path);
    } catch (Exception ex) {
      ex.printStackTrace();
    }
  }
}
