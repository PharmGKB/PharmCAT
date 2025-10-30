package org.pharmgkb.pharmcat;

import java.io.BufferedReader;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.HashMap;
import java.util.Map;
import java.util.Optional;
import java.util.Set;
import java.util.TreeMap;
import com.google.common.collect.HashMultimap;
import com.google.common.collect.Multimap;
import org.jspecify.annotations.Nullable;
import org.pharmgkb.pharmcat.definition.DefinitionReader;
import org.pharmgkb.pharmcat.definition.model.DefinitionFile;
import org.pharmgkb.pharmcat.phenotype.PhenotypeMap;
import org.pharmgkb.pharmcat.phenotype.model.GenePhenotype;
import org.pharmgkb.pharmcat.reporter.MessageHelper;
import org.pharmgkb.pharmcat.reporter.PgkbGuidelineCollection;
import org.pharmgkb.pharmcat.reporter.caller.Cyp2d6CopyNumberCaller;
import org.pharmgkb.pharmcat.reporter.model.DataSource;
import org.pharmgkb.pharmcat.reporter.model.MessageAnnotation;
import org.pharmgkb.pharmcat.reporter.model.result.Haplotype;
import org.pharmgkb.pharmcat.util.DataManager;


/**
 * Global environment for PharmCAT.
 *
 * @author Mark Woon
 */
public class Env {
  private final DefinitionReader m_definitionReader;
  private final PhenotypeMap m_phenotypeMap;
  private final PgkbGuidelineCollection m_drugs;
  private final MessageHelper m_messageHelper;
  private final Map<String, Map<String, Haplotype>> m_haplotypeCache = new HashMap<>();
  private final Multimap<String, String> m_validHaplotypes = HashMultimap.create();
  private final Map<Path, Map<String, Map<String, String>>> m_sampleDataMap = new HashMap<>();


  public Env() throws IOException, ReportableException {
    this(null, null);
  }

  public Env(@Nullable Path definitionDir, @Nullable Set<String> genes) throws IOException, ReportableException {
    if (definitionDir == null) {
      definitionDir = DataManager.DEFAULT_DEFINITION_DIR;
    }
    m_definitionReader = new DefinitionReader(definitionDir, genes);
    if (m_definitionReader.getGenes().isEmpty()) {
      throw new ReportableException("Did not find any allele definitions at " + definitionDir);
    }
    try {
      m_messageHelper = new MessageHelper();
    } catch (IOException ex) {
      throw new RuntimeException("Error loading messages", ex);
    }

    m_phenotypeMap = new PhenotypeMap();
    m_drugs = new PgkbGuidelineCollection();

    // initialize dependent classes
    Cyp2d6CopyNumberCaller.initialize(this);
  }


  public DefinitionReader getDefinitionReader() {
    return m_definitionReader;
  }

  public String getReferenceAllele(String gene) {
    return m_definitionReader.getReferenceAlleleMap().get(gene);
  }


  public @Nullable String getPhenotypeVersion(String gene) {
    return m_phenotypeMap.getVersion(gene);
  }

  public @Nullable GenePhenotype getPhenotype(String gene) {
    return m_phenotypeMap.getPhenotype(gene);
  }

  /**
   * Checks if the specified allele is used in either definition files or phenotype.
   */
  public boolean isValidNamedAllele(String gene, String allele) {

    if (m_validHaplotypes.containsEntry(gene, allele)) {
      return true;
    }

    if (gene.startsWith("HLA-")) {
      // HLA's are a special case
      m_validHaplotypes.put(gene, allele);
      return true;
    }

    Optional<DefinitionFile> opt = m_definitionReader.lookupDefinitionFile(gene);
    if (opt.isPresent() && opt.get().getNamedAllele(allele) != null) {
      m_validHaplotypes.put(gene, allele);
      return true;
    }

    String inferredAllele = allele;
    if (gene.equals("CYP2D6")) {
      inferredAllele = Cyp2d6CopyNumberCaller.inferHaplotypeName(allele);
    }
    GenePhenotype gp = m_phenotypeMap.getPhenotype(gene);
    if (gp != null) {
      if (gp.getHaplotypes().containsKey(inferredAllele) || gp.getActivityValues().containsKey(inferredAllele)) {
        m_validHaplotypes.put(gene, allele);
        return true;
      }
    }
    return false;
  }


  public PgkbGuidelineCollection getDrugs() {
    return m_drugs;
  }


  /**
   * Checks if gene is used in any guideline from the specified {@code source}.
   */
  public boolean hasGene(DataSource source, String gene) {
    return m_drugs.getGenesUsedInSource(source).contains(gene);
  }

  /**
   * Checks if gene can be called by NamedAlleleMatcher or is used in any drug recommendation from any source.
   */
  @SuppressWarnings("BooleanMethodIsAlwaysInverted")
  public boolean hasGene(String gene) {
    return m_definitionReader.getGenes().contains(gene) ||
        m_drugs.getGenesWithRecommendations().contains(gene);
  }

  /**
   * Checks if this gene is matched by activity score in any source.
   */
  public boolean isActivityScoreGene(String gene) {
    GenePhenotype gp = getPhenotype(gene);
    return gp != null && gp.isMatchedByActivityScore();
  }


  public MessageHelper getMessageHelper() {
    return m_messageHelper;
  }

  public @Nullable MessageAnnotation getMessage(String key) {
    return m_messageHelper.getMessage(key);
  }


  /**
   * Make or retrieve a cached {@link Haplotype} object that corresponds to the given allele name.
   */
  public synchronized Haplotype makeHaplotype(String gene, String name) {
    return m_haplotypeCache
        .computeIfAbsent(gene, (g) -> new HashMap<>())
        .computeIfAbsent(name, (n) -> {
          Haplotype haplotype = new Haplotype(gene, name);
          GenePhenotype gp = getPhenotype(gene);
          if (gp != null) {
            haplotype.setFunction(gp.getHaplotypeFunction(name));
            haplotype.setActivityValue(gp.getHaplotypeActivity(name));
          }
          haplotype.setReference(name.equals(getReferenceAllele(gene)));
          return haplotype;
        });
  }


  public synchronized @Nullable Map<String, String> getSampleMetadata(Path sampleMetadataFile, String sampleId,
      boolean cache) throws IOException {

    Map<String, Map<String, String>> fileMap;
    if (cache) {
      fileMap = m_sampleDataMap.computeIfAbsent(sampleMetadataFile, f -> new TreeMap<>());
      if (!fileMap.isEmpty()) {
        return fileMap.get(sampleId);
      }
    } else {
      fileMap = new TreeMap<>();
    }
    try (BufferedReader reader = Files.newBufferedReader(sampleMetadataFile)) {
      String line;
      while ((line = reader.readLine()) != null) {
        String[] row = line.split("\t");
        if (row.length >= 3) {
          String sid = row[0];
          if (cache) {
            fileMap.computeIfAbsent(sid, k -> new HashMap<>())
                .put(row[1], row[2]);
          } else {
            if (sid.equals(sampleId)) {
              fileMap.computeIfAbsent(sid, k -> new HashMap<>())
                  .put(row[1], row[2]);
            } else if (fileMap.containsKey(sid)) {
              // all values for a single sample must be consecutive
              break;
            }
          }
        }
      }
    }
    return fileMap.get(sampleId);
  }
}
