package org.pharmgkb.pharmcat.definition;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.Arrays;
import java.util.List;
import java.util.Optional;
import com.google.gson.Gson;
import org.apache.commons.lang3.StringUtils;
import org.pharmgkb.pharmcat.definition.model.GenePhenotype;


/**
 * This class loads and manages the data from the gene phenotypes file
 *
 * @author Ryan Whaley
 */
public class PhenotypeMap {

  private static final String FILE_NAME = "gene.phenotypes.json";
  private final List<GenePhenotype> m_genes;

  /**
   * public constructor, loads the data from a local file
   */
  public PhenotypeMap() throws Exception {
    try (BufferedReader reader = new BufferedReader(new InputStreamReader(getClass().getResourceAsStream(FILE_NAME)))) {
      Gson gson = new Gson();
      m_genes = Arrays.asList(gson.fromJson(reader, GenePhenotype[].class));
    } catch (IOException e) {
      throw new Exception("Error reading phenotype definitions", e);
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

  public Optional<String> lookupPhenotype(String genotype) {
    if (StringUtils.isBlank(genotype) || !genotype.contains(":")) return Optional.empty();
    String[] tokens = genotype.split(":");

    GenePhenotype genePhenotype = lookup(tokens[0]).orElse(null);
    if (genePhenotype != null) {
      String result = genePhenotype.getDiplotypeResults().get(genotype);
      if (result != null) {
        return Optional.of(tokens[0] + ":" + result);
      } else {
        return Optional.empty();
      }
    } else {
      return Optional.empty();
    }
  }
}
