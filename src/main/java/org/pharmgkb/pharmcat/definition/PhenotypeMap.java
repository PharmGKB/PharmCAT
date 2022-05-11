package org.pharmgkb.pharmcat.definition;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Objects;
import java.util.Optional;
import com.google.gson.Gson;
import org.apache.commons.lang3.StringUtils;
import org.pharmgkb.pharmcat.definition.model.GenePhenotype;


/**
 * This class loads and manages the data from the gene phenotypes file.
 *
 * @author Ryan Whaley
 */
public class PhenotypeMap {

  private static final String FILE_NAME = "gene_phenotypes.json";
  private static final String DPWG_NAME = "dpwg_phenotypes.json";
  private final List<GenePhenotype> m_genes;

  /**
   * public constructor, loads the data from a local file
   */
  public PhenotypeMap() {
    try (
        BufferedReader reader = new BufferedReader(new InputStreamReader(Objects.requireNonNull(getClass().getResourceAsStream(FILE_NAME))));
        BufferedReader dpwgReader = new BufferedReader(new InputStreamReader(Objects.requireNonNull(getClass().getResourceAsStream(DPWG_NAME))))
    ) {
      Gson gson = new Gson();
      m_genes = new ArrayList<>();
      m_genes.addAll(Arrays.asList(gson.fromJson(reader, GenePhenotype[].class)));
      m_genes.addAll(Arrays.asList(gson.fromJson(dpwgReader, GenePhenotype[].class)));
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

  protected List<GenePhenotype> getCpicGenes() {
    return m_genes.stream().filter(g -> Objects.isNull(g.getDrugId())).toList();
  }

  /**
   * Lookup and return the {@link Optional} {@link GenePhenotype} object for the given gene symbol
   * @param gene an HGNC gene symbol
   */
  public Optional<GenePhenotype> lookup(String gene) {
    return m_genes.stream().filter(p -> p.getDrugId() == null && gene.equals(p.getGene())).findFirst();
  }
}
