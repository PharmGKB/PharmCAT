package org.pharmgkb.pharmcat.definition;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.Arrays;
import java.util.List;
import java.util.Optional;
import com.google.gson.Gson;
import org.pharmgkb.pharmcat.definition.model.GenePhenotype;


/**
 * This class loads and manages the data from the gene phenotypes file
 *
 * @author Ryan Whaley
 */
public class PhenotypeMap {

  private static final String FILE_NAME = "gene_phenotypes.json";
  private final List<GenePhenotype> m_genes;

  /**
   * public constructor, loads the data from a local file
   */
  public PhenotypeMap() {
    try (BufferedReader reader = new BufferedReader(new InputStreamReader(getClass().getResourceAsStream(FILE_NAME)))) {
      Gson gson = new Gson();
      m_genes = Arrays.asList(gson.fromJson(reader, GenePhenotype[].class));
    } catch (IOException e) {
      throw new RuntimeException("Error reading phenotype definitions", e);
    }
  }

  public PhenotypeMap(Path filePath) {
    try (BufferedReader reader = new BufferedReader(new InputStreamReader(Files.newInputStream(filePath)))) {
      Gson gson = new Gson();
      m_genes = Arrays.asList(gson.fromJson(reader, GenePhenotype[].class));
    } catch (IOException e) {
      throw new RuntimeException("Error reading phenotype definitions", e);
    }
  }

  protected List<GenePhenotype> getGenes() {
    return m_genes;
  }

  /**
   * Lookup and return the {@link Optional} {@link GenePhenotype} object for the given gene symbol
   * @param gene an HGNC gene symbol
   */
  public Optional<GenePhenotype> lookup(String gene) {
    return m_genes.stream().filter(p -> gene.equals(p.getGene())).findFirst();
  }
}
