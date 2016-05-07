package org.pharmgkb.pharmcat.haplotype;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.List;
import java.util.Set;
import java.util.SortedMap;
import java.util.TreeMap;
import java.util.stream.Collectors;
import javax.annotation.Nonnull;
import com.google.common.base.Preconditions;
import org.pharmgkb.pharmcat.definition.GeneratedDefinitionSerializer;
import org.pharmgkb.pharmcat.definition.model.DefinitionFile;
import org.pharmgkb.pharmcat.definition.model.NamedAllele;
import org.pharmgkb.pharmcat.definition.model.VariantLocus;


/**
 * Read in haplotype definition files.
 *
 * @author Mark Woon
 */
public class DefinitionReader {
  private GeneratedDefinitionSerializer m_definitionSerializer = new GeneratedDefinitionSerializer();
  private SortedMap<String, DefinitionFile> m_definitionFiles = new TreeMap<>();


  public @Nonnull Set<String> getGenes() {
    return m_definitionFiles.keySet();
  }


  public @Nonnull DefinitionFile getDefinitionFile(String gene) {
    Preconditions.checkArgument(m_definitionFiles.containsKey(gene));
    return m_definitionFiles.get(gene);
  }


  public @Nonnull VariantLocus[] getPositions(String gene) {
    return m_definitionFiles.get(gene).getVariants();
  }


  public @Nonnull List<NamedAllele> getHaplotypes(String gene) {
    Preconditions.checkArgument(m_definitionFiles.containsKey(gene));
    return m_definitionFiles.get(gene).getNamedAlleles();
  }



  public void read(Path path) throws IOException {

    if (Files.isDirectory(path)) {
      List<Path> files = Files.list(path)
          .filter(f -> f.toString().endsWith(".json"))
          .collect(Collectors.toList());
      for (Path file : files) {
        readFile(file);
      }
    } else {
      readFile(path);
    }
  }


  private void readFile(@Nonnull Path file) throws IOException {

    Preconditions.checkNotNull(file);
    Preconditions.checkArgument(Files.isRegularFile(file));
    DefinitionFile definitionFile = m_definitionSerializer.deserializeFromJson(file);

    String gene = definitionFile.getGeneSymbol();
    m_definitionFiles.put(gene, definitionFile);
  }
}
