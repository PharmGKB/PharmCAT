package org.pharmgkb.pharmcat.phenotype;

import java.io.BufferedReader;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Optional;
import java.util.stream.Stream;
import com.google.common.base.Preconditions;
import org.jspecify.annotations.Nullable;
import org.pharmgkb.common.util.PathUtils;
import org.pharmgkb.pharmcat.phenotype.model.GenePhenotype;
import org.pharmgkb.pharmcat.util.DataSerializer;


/**
 * This class loads and manages the data from the gene phenotypes file.
 *
 * @author Ryan Whaley
 */
public class PhenotypeMap {
  private static final Path sf_phenotypesDir =
      // have to resolve specific file, then backtrack to get directory
      PathUtils.getPathToResource("org/pharmgkb/pharmcat/phenotype/CYP2C19.json").getParent();
  private final Map<String, GenePhenotype> m_geneMap = new HashMap<>();


  /**
   * public constructor, loads the data from a local file
   */
  public PhenotypeMap() {
    this(sf_phenotypesDir);
  }

  public PhenotypeMap(Path dir) {
    Preconditions.checkArgument(Files.isDirectory(dir));
    try {
      initialize(dir, m_geneMap);
    } catch (IOException ex) {
      throw new RuntimeException("Error reading phenotype data", ex);
    }
  }

  private static void initialize(Path sourceDir, Map<String, GenePhenotype> sourceMap) throws IOException {
    List<Path> phenotypeFiles = new ArrayList<>();
    try (Stream<Path> stream = Files.list(sourceDir)) { {
      stream.filter(f -> f.getFileName().toString().endsWith(".json"))
          .forEach(phenotypeFiles::add);
    }}
    if (phenotypeFiles.isEmpty()) {
      throw new IOException("Cannot find phenotype files");
    }
    for (Path phenotypeFile : phenotypeFiles) {
      try (BufferedReader br = Files.newBufferedReader(phenotypeFile)) {
        GenePhenotype gp = DataSerializer.GSON.fromJson(br, GenePhenotype.class);
        if (sourceMap.put(gp.getGene(), gp) != null) {
          throw new IllegalStateException("Multiple GenePhenotypes for " + gp.getGene());
        }
      }
    }
  }


  public @Nullable String getVersion(String gene) {
    GenePhenotype gp = getPhenotype(gene);
    if (gp != null) {
      return gp.getVersion();
    }
    return null;
  }

  public Collection<GenePhenotype> getGenePhenotypes() {
    return m_geneMap.values();
  }


  /**
   * Lookup the phenotypes for the specified {@code gene} and {@code source}.
   */
  public @Nullable GenePhenotype getPhenotype(String gene) {
    return m_geneMap.get(gene);
  }


  /**
   * Lookup and return the {@link Optional} {@link GenePhenotype} object for the given gene symbol
   * @param gene an HGNC gene symbol
   */
  public Optional<GenePhenotype> getPhenotypeOptional(String gene) {
    return Optional.ofNullable(m_geneMap.get(gene));
  }
}
