package org.pharmgkb.pharmcat.definition;

import java.io.BufferedReader;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.Collection;
import java.util.HashMap;
import java.util.Map;
import java.util.Optional;
import com.google.common.base.Preconditions;
import com.google.gson.Gson;
import org.checkerframework.checker.nullness.qual.Nullable;
import org.pharmgkb.common.util.PathUtils;
import org.pharmgkb.pharmcat.definition.model.GenePhenotype;
import org.pharmgkb.pharmcat.reporter.model.DataSource;


/**
 * This class loads and manages the data from the gene phenotypes file.
 *
 * @author Ryan Whaley
 */
public class PhenotypeMap {
  public static final String CPIC_PHENOTYPES_JSON_FILE_NAME = "gene_phenotypes.json";
  private static final Path sf_genePhenotypesFile =
      PathUtils.getPathToResource("org/pharmgkb/pharmcat/definition/" + CPIC_PHENOTYPES_JSON_FILE_NAME);
  public static final String DPWG_PHENOTYPES_JSON_FILE_NAME = "dpwg_phenotypes.json";
  private static final Path sf_dpwgPhenotypesFile =
      PathUtils.getPathToResource("org/pharmgkb/pharmcat/definition/" + DPWG_PHENOTYPES_JSON_FILE_NAME);
  private final Map<String, GenePhenotype> m_cpicMap = new HashMap<>();
  private final Map<String, GenePhenotype> m_dpwgMap = new HashMap<>();


  /**
   * public constructor, loads the data from a local file
   */
  public PhenotypeMap() {
    try (BufferedReader cpicReader = Files.newBufferedReader(sf_genePhenotypesFile);
         BufferedReader dpwgReader = Files.newBufferedReader(sf_dpwgPhenotypesFile)
    ) {
      initialize(cpicReader, dpwgReader);
    } catch (IOException e) {
      throw new RuntimeException("Error reading phenotype definitions", e);
    }
  }

  public PhenotypeMap(Path dir) {
    Preconditions.checkArgument(Files.isDirectory(dir));
    try (BufferedReader cpicReader = Files.newBufferedReader(dir.resolve(CPIC_PHENOTYPES_JSON_FILE_NAME));
         BufferedReader dpwgReader = Files.newBufferedReader(dir.resolve(DPWG_PHENOTYPES_JSON_FILE_NAME))
    ) {
      initialize(cpicReader, dpwgReader);
    } catch (IOException e) {
      throw new RuntimeException("Error reading phenotype definitions", e);
    }
  }


  private void initialize(BufferedReader cpicReader, BufferedReader dpwgReader) {
    Gson gson = new Gson();
    GenePhenotype[] rez = gson.fromJson(cpicReader, GenePhenotype[].class);
    for (GenePhenotype gp : rez) {
      if (m_cpicMap.put(gp.getGene(), gp) != null) {
        throw new IllegalStateException("Multiple CPIC GenePhenotypes for " + gp.getGene());
      }
    }
    rez = gson.fromJson(dpwgReader, GenePhenotype[].class);
    for (GenePhenotype gp : rez) {
      if (m_dpwgMap.put(gp.getGene(), gp) != null) {
        throw new IllegalStateException("Multiple DPWG GenePhenotypes for " + gp.getGene());
      }
    }
  }

  public @Nullable String getVersion(String gene, DataSource source) {
    GenePhenotype gp = lookupPhenotype(gene, source);
    if (gp != null) {
      return gp.getVersion();
    }
    return null;

  }

  public Collection<GenePhenotype> getCpicGenes() {
    return m_cpicMap.values();
  }

  public Collection<GenePhenotype> getDpwgGenes() {
    return m_dpwgMap.values();
  }


  /**
   * Lookup the phenotypes for the specified {@code gene} and {@code source}.
   */
  public @Nullable GenePhenotype lookupPhenotype(String gene, DataSource source) {
    switch (source) {
      case CPIC -> {
        return m_cpicMap.get(gene);
      }
      case DPWG -> {
        return m_dpwgMap.get(gene);
      }
    }
    return null;
  }


  /**
   * Lookup and return the {@link Optional} {@link GenePhenotype} object for the given gene symbol
   * @param gene an HGNC gene symbol
   */
  public Optional<GenePhenotype> lookupCpic(String gene) {
    return Optional.ofNullable(m_cpicMap.get(gene));
  }
}
