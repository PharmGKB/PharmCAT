package org.pharmgkb.pharmcat.haplotype;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.SortedMap;
import java.util.SortedSet;
import java.util.TreeMap;
import com.google.common.base.Preconditions;
import org.checkerframework.checker.nullness.qual.Nullable;
import org.pharmgkb.pharmcat.definition.model.DefinitionExemption;
import org.pharmgkb.pharmcat.definition.model.DefinitionFile;
import org.pharmgkb.pharmcat.definition.model.NamedAllele;
import org.pharmgkb.pharmcat.definition.model.VariantLocus;
import org.pharmgkb.pharmcat.util.DataManager;
import org.pharmgkb.pharmcat.util.DataSerializer;


/**
 * Read in haplotype definition files.
 *
 * @author Mark Woon
 */
public class DefinitionReader {
  private final DataSerializer m_definitionSerializer = new DataSerializer();
  private final SortedMap<String, DefinitionFile> m_definitionFiles = new TreeMap<>();
  private final Map<String, DefinitionExemption> m_exemptions = new TreeMap<>();
  private String m_genomeBuild;


  /**
   * Gets the genome build used by the allele definitions.
   * This should be called <em>after</em> all allele definitions have been read.
   */
  public String getGenomeBuild() {
    Preconditions.checkState(m_definitionFiles.size() > 0);

    if (m_genomeBuild == null) {
      for (DefinitionFile definitionFile : m_definitionFiles.values()) {
        if (m_genomeBuild == null) {
          m_genomeBuild = definitionFile.getGenomeBuild();
        } else if (!m_genomeBuild.equalsIgnoreCase(definitionFile.getGenomeBuild())) {
          throw new IllegalStateException("Definition files use different genome builds (" + m_genomeBuild + " vs " +
              definitionFile.getGenomeBuild() + " for " + definitionFile.getGeneSymbol() + ")");
        }
      }
    }
    return m_genomeBuild;
  }


  public Set<String> getGenes() {
    return m_definitionFiles.keySet();
  }

  public Map<String,Integer> getGeneAlleleCount() {
    Map<String,Integer> countMap = new TreeMap<>();
    for (String gene : getGenes()) {
      DefinitionFile file = m_definitionFiles.get(gene);
      countMap.put(gene, file.getNamedAlleles().size());
    }
    return countMap;
  }


  public DefinitionFile getDefinitionFile(String gene) {
    Preconditions.checkArgument(m_definitionFiles.containsKey(gene));
    return m_definitionFiles.get(gene);
  }


  public VariantLocus[] getPositions(String gene) {
    return m_definitionFiles.get(gene).getVariants();
  }


  public SortedSet<NamedAllele> getHaplotypes(String gene) {
    Preconditions.checkArgument(m_definitionFiles.containsKey(gene));
    return m_definitionFiles.get(gene).getNamedAlleles();
  }

  public @Nullable DefinitionExemption getExemption(String gene) {
    return m_exemptions.get(gene.toLowerCase());
  }



  public void read(Path path) throws IOException {

    if (Files.isDirectory(path)) {
      List<Path> files = Files.list(path)
          .filter(f -> f.toString().endsWith("_translation.json"))
          .toList();
      for (Path file : files) {
        readFile(file);
      }
      readExemptions(path);
    } else {
      readFile(path);
    }
  }


  private void readFile(Path file) throws IOException {

    Preconditions.checkNotNull(file);
    Preconditions.checkArgument(Files.isRegularFile(file), "%s is not a file", file);
    DefinitionFile definitionFile = m_definitionSerializer.deserializeDefinitionsFromJson(file);

    String gene = definitionFile.getGeneSymbol();
    m_definitionFiles.put(gene, definitionFile);
  }


  public void readExemptions(Path path) throws IOException {

    Preconditions.checkNotNull(path);
    Path file;
    if (Files.isDirectory(path)) {
      file = path.resolve(DataManager.EXEMPTIONS_JSON_FILE_NAME);
    } else {
      file = path;
    }
    Preconditions.checkArgument(Files.isRegularFile(file), "Not a file: %s", file);

    Set<DefinitionExemption> exemptions = m_definitionSerializer.deserializeExemptionsFromJson(file);
    for (DefinitionExemption de : exemptions) {
      m_exemptions.put(de.getGene().toLowerCase(), de);
    }
  }
}
