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
import java.util.stream.Collectors;
import javax.annotation.Nonnull;
import javax.annotation.Nullable;
import com.google.common.base.Preconditions;
import com.google.common.base.Splitter;
import com.google.common.collect.Sets;
import org.pharmgkb.pharmcat.definition.DefinitionManager;
import org.pharmgkb.pharmcat.definition.GeneratedDefinitionSerializer;
import org.pharmgkb.pharmcat.definition.model.DefinitionExemption;
import org.pharmgkb.pharmcat.definition.model.DefinitionFile;
import org.pharmgkb.pharmcat.definition.model.NamedAllele;
import org.pharmgkb.pharmcat.definition.model.VariantLocus;


/**
 * Read in haplotype definition files.
 *
 * @author Mark Woon
 */
public class DefinitionReader {
  private static final Splitter sf_commaSplitter = Splitter.on(",").trimResults().omitEmptyStrings();
  private GeneratedDefinitionSerializer m_definitionSerializer = new GeneratedDefinitionSerializer();
  private SortedMap<String, DefinitionFile> m_definitionFiles = new TreeMap<>();
  private Map<String, DefinitionExemption> m_exemptions = new TreeMap<>();
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

  public @Nullable DefinitionExemption getExemption(@Nonnull String gene) {
    return m_exemptions.get(gene.toLowerCase());
  }



  public void read(Path path) throws IOException {

    if (Files.isDirectory(path)) {
      List<Path> files = Files.list(path)
          .filter(f -> f.toString().endsWith(".json"))
          .collect(Collectors.toList());
      for (Path file : files) {
        readFile(file);
      }
      readExemptions(path);
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


  public void readExemptions(@Nonnull Path path) throws IOException {

    Preconditions.checkNotNull(path);
    Path file;
    if (Files.isDirectory(path)) {
      file = path.resolve(DefinitionManager.EXEMPTIONS_FILE_NAME);
    } else {
      file = path;
    }
    Preconditions.checkArgument(Files.isRegularFile(file), "Not a file: %s", file);

     Set<DefinitionExemption> exemptions = Files.lines(file)
        .skip(1) // skip the header
        .map(line -> {
          String[] data = line.split("\t");
          String gene = data[0];
          SortedSet<String> ignoreAlleles = Sets.newTreeSet(sf_commaSplitter.splitToList(data[1]));
          boolean allHits = data.length > 2 && Boolean.parseBoolean(data[2]);
          return new DefinitionExemption(gene, ignoreAlleles, allHits);
        })
        .collect(Collectors.toSet());
    for (DefinitionExemption de : exemptions) {
      m_exemptions.put(de.getGene().toLowerCase(), de);
    }
  }
}
