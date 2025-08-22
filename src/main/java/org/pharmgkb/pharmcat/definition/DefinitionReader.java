package org.pharmgkb.pharmcat.definition;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.*;
import java.util.stream.Stream;
import com.google.common.base.Preconditions;
import com.google.common.collect.ImmutableMap;
import org.jspecify.annotations.Nullable;
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
  private final Map<String, DefinitionExemption> m_exemptions;
  private String m_genomeBuild;
  private ReferenceAlleleMap m_referenceAlleleMap;
  /** Map of {@code <chr:position>} Strings to {@link VariantLocus} */
  private ImmutableMap<String, VariantLocus> m_locationsOfInterest;
  /** Map of {@code <chr:position>} Strings to gene */
  private ImmutableMap<String, String> m_locationsByGene;


  public DefinitionReader() throws IOException {
    this(DataManager.DEFAULT_DEFINITION_DIR);
  }

  public DefinitionReader(Path dir) throws IOException {
    Preconditions.checkArgument(Files.isDirectory(dir));

    try (Stream<Path> fileStream = Files.list(dir)) {
      List<Path> files = fileStream.filter(f -> f.toString().endsWith("_translation.json"))
          .toList();
      for (Path file : files) {
        readFile(file);
      }
    }
    m_exemptions = readExemptions(dir);
    generateMetadata();
  }

  public DefinitionReader(Path definitionFile, @Nullable Path exemptionsFile) throws IOException {
    this(List.of(definitionFile), exemptionsFile);
  }

  public DefinitionReader(List<Path> definitionFiles, @Nullable Path exemptionsFile) throws IOException {
    Preconditions.checkNotNull(definitionFiles);

    for (Path file : definitionFiles) {
      if (!Files.isRegularFile(file)) {
        throw new IllegalArgumentException(file + " is not a file");
      }
      readFile(file);
    }
    if (exemptionsFile != null) {
      m_exemptions = readExemptions(exemptionsFile);
    } else {
      m_exemptions = Collections.emptyMap();
    }
    generateMetadata();
  }


  /**
   * Gets the genome build used by the allele definitions.
   * This should be called <em>after</em> all allele definitions have been read.
   */
  public String getGenomeBuild() {
    Preconditions.checkState(!m_definitionFiles.isEmpty());

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
    Preconditions.checkArgument(m_definitionFiles.containsKey(gene), "No definition file for " + gene);
    return m_definitionFiles.get(gene);
  }

  public Optional<DefinitionFile> lookupDefinitionFile(String gene) {
    return Optional.ofNullable(m_definitionFiles.get(gene));
  }


  public VariantLocus[] getPositions(String gene) {
    return m_definitionFiles.get(gene).getVariants();
  }


  public ImmutableMap<String, VariantLocus> getLocationsOfInterest() {
    return m_locationsOfInterest;
  }

  public ImmutableMap<String, String> getLocationsByGene() {
    return m_locationsByGene;
  }


  public SortedSet<NamedAllele> getHaplotypes(String gene) {
    Preconditions.checkArgument(m_definitionFiles.containsKey(gene));
    return m_definitionFiles.get(gene).getNamedAlleles();
  }

  public NamedAllele getReferenceHaplotype(String gene) {
    return m_definitionFiles.get(gene).getReferenceNamedAllele();
  }

  public @Nullable DefinitionExemption getExemption(String gene) {
    return m_exemptions.get(gene);
  }


  public ReferenceAlleleMap getReferenceAlleleMap() {
    if (m_referenceAlleleMap == null) {
      m_referenceAlleleMap = new ReferenceAlleleMap(this);
    }
    return m_referenceAlleleMap;
  }


  private void generateMetadata() {

    Set<String> data = new HashSet<>();
    ImmutableMap.Builder<String, VariantLocus> vlMapBuilder = ImmutableMap.builder();
    ImmutableMap.Builder<String, String> geneMapBuilder = ImmutableMap.builder();
    for (String gene : m_definitionFiles.keySet()) {
      Arrays.stream(m_definitionFiles.get(gene).getVariants())
          .forEach(v -> {
            String vcp = v.getVcfChrPosition();
            data.add(vcp);
            vlMapBuilder.put(vcp, v);
            geneMapBuilder.put(vcp, gene);
          });
      DefinitionExemption exemption = m_exemptions.get(gene);
      if (exemption != null) {
        exemption.getExtraPositions()
            .forEach(v -> {
              String vcp = v.getVcfChrPosition();
              if (!data.contains(vcp)) {
                data.add(vcp);
                vlMapBuilder.put(vcp, v);
                geneMapBuilder.put(vcp, gene);
              }
            });
      }
    }
    m_locationsOfInterest = vlMapBuilder.build();
    m_locationsByGene = geneMapBuilder.build();
  }


  private void readFile(Path file) throws IOException {

    Preconditions.checkNotNull(file);
    Preconditions.checkArgument(Files.isRegularFile(file), "%s is not a file", file);
    DefinitionFile definitionFile = m_definitionSerializer.deserializeDefinitionsFromJson(file);

    String gene = definitionFile.getGeneSymbol();
    m_definitionFiles.put(gene, definitionFile);
  }


  private Map<String, DefinitionExemption> readExemptions(Path path) throws IOException {

    Preconditions.checkNotNull(path);
    Path file;
    if (Files.isDirectory(path)) {
      file = path.resolve(DataManager.EXEMPTIONS_JSON_FILE_NAME);
    } else {
      file = path;
    }
    Preconditions.checkArgument(Files.isRegularFile(file), "Not a file: %s", file);

    return m_definitionSerializer.deserializeExemptionsFromJson(file);
  }


  private static DefinitionReader s_defaultReader;

  public static DefinitionReader defaultReader() throws IOException {
    if (s_defaultReader == null) {
      s_defaultReader = new DefinitionReader();
    }
    return s_defaultReader;
  }
}
