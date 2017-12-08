package org.pharmgkb.pharmcat.definition;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.Arrays;
import java.util.List;
import java.util.Optional;
import java.util.Set;
import javax.annotation.Nullable;
import com.google.gson.Gson;
import org.pharmgkb.pharmcat.UnexpectedStateException;
import org.pharmgkb.pharmcat.definition.model.GenePhenotype;
import org.pharmgkb.pharmcat.haplotype.model.GeneCall;
import org.pharmgkb.pharmcat.haplotype.model.HaplotypeMatch;


/**
 * This class loads and manages the data from the gene phenotypes file
 *
 * @author Ryan Whaley
 */
public class PhenotypeMap {

  private List<GenePhenotype> m_genes;

  /**
   * public constructor, loads the data from a local file
   */
  public PhenotypeMap(@Nullable List<GeneCall> calls) throws Exception {
    try (BufferedReader reader = new BufferedReader(new InputStreamReader(getClass().getResourceAsStream("gene.phenotypes.json")))) {

      Gson gson = new Gson();
      m_genes = Arrays.asList(gson.fromJson(reader, GenePhenotype[].class));

      if (calls != null) {
        calls.forEach(c -> {
          String gene = c.getGene();
          Set<HaplotypeMatch> haplotypematches = c.getHaplotypes();
          Optional<GenePhenotype> lookupResult = lookup(gene);

          if (lookupResult.isPresent()) {
            GenePhenotype genePhenotype = lookupResult.get();
            for (HaplotypeMatch haplotypeMatch : haplotypematches) {
              String hap = haplotypeMatch.getName();
              String newFunction = haplotypeMatch.getFunction();
              String existingFunction = genePhenotype.lookupHaplotype(hap);

              if (existingFunction != null && !existingFunction.equalsIgnoreCase(newFunction)) {
                throw new UnexpectedStateException("Function mismatch for " + gene + " " + hap + " > " + existingFunction + " != " + newFunction);
              }
              if (existingFunction == null) {
                genePhenotype.addHaplotypeFunction(hap, newFunction);
              }
            }
          }
        });
      }

    } catch (IOException e) {
      throw new Exception("Error reading phenotype definitions", e);
    }
  }

  protected List getGenes() {
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
