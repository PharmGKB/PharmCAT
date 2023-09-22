package org.pharmgkb.pharmcat;

import java.io.IOException;
import java.nio.file.Path;
import java.util.HashMap;
import java.util.Map;
import org.checkerframework.checker.nullness.qual.Nullable;
import org.pharmgkb.pharmcat.definition.DefinitionReader;
import org.pharmgkb.pharmcat.phenotype.PhenotypeMap;
import org.pharmgkb.pharmcat.phenotype.model.GenePhenotype;
import org.pharmgkb.pharmcat.reporter.MessageHelper;
import org.pharmgkb.pharmcat.reporter.PgkbGuidelineCollection;
import org.pharmgkb.pharmcat.reporter.caller.Cyp2d6CopyNumberCaller;
import org.pharmgkb.pharmcat.reporter.model.DataSource;
import org.pharmgkb.pharmcat.reporter.model.MessageAnnotation;
import org.pharmgkb.pharmcat.reporter.model.result.Haplotype;


/**
 * Global environment for PharmCAT.
 *
 * @author Mark Woon
 */
public class Env {
  private final DefinitionReader m_definitionReader;
  private final PhenotypeMap m_phenotypeMap;
  private final PgkbGuidelineCollection m_drugs;
  private MessageHelper m_messageHelper;
  private final Map<DataSource, Map<String, Map<String, Haplotype>>> m_haplotypeCache = new HashMap<>();


  public Env() throws IOException, ReportableException {
    this(null);
  }

  public Env(@Nullable Path definitionDir) throws IOException, ReportableException {
    if (definitionDir != null) {
      m_definitionReader = new DefinitionReader(definitionDir);
      if (m_definitionReader.getGenes().isEmpty()) {
        throw new ReportableException("Did not find any allele definitions at " + definitionDir);
      }
    } else {
      m_definitionReader = DefinitionReader.defaultReader();
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


  public @Nullable String getPhenotypeVersion(String gene, DataSource source) {
    return m_phenotypeMap.getVersion(gene, source);
  }

  public @Nullable GenePhenotype getPhenotype(String gene, DataSource source) {
    return m_phenotypeMap.getPhenotype(gene, source);
  }


  public PgkbGuidelineCollection getDrugs() {
    return m_drugs;
  }


  public boolean hasGene(DataSource source, String gene) {
    return m_drugs.getGenesUsedInSource(source).contains(gene);
  }

  /**
   * Checks if gene can be called by NamedAlleleMatcher or is used in any drug recommendation from any source.
   */
  public boolean hasGene(String gene) {
    return m_definitionReader.getGenes().contains(gene) ||
        m_drugs.getGenesWithRecommendations().contains(gene);
  }

  /**
   * Checks if this gene is matched by activity score in any source.
   */
  public boolean isActivityScoreGene(String gene) {
    GenePhenotype gp = getPhenotype(gene, DataSource.CPIC);
    if (gp != null && gp.isMatchedByActivityScore()) {
      return true;
    }
    gp = getPhenotype(gene, DataSource.DPWG);
    return gp != null && gp.isMatchedByActivityScore();
  }


  public MessageHelper getMessageHelper() {
    if (m_messageHelper == null) {
      try {
        m_messageHelper = new MessageHelper();
      } catch (IOException ex) {
        throw new RuntimeException("Error loading messages", ex);
      }
    }
    return m_messageHelper;
  }

  public MessageAnnotation getMessage(String key) {
    return getMessageHelper().getMessage(key);
  }


  /**
   * Make or retrieve a cached {@link Haplotype} object that corresponds to the given allele name.
   */
  public synchronized Haplotype makeHaplotype(String gene, String name, DataSource source) {
    return m_haplotypeCache.computeIfAbsent(source, (s) -> new HashMap<>())
        .computeIfAbsent(gene, (g) -> new HashMap<>())
        .computeIfAbsent(name, (n) -> {
          Haplotype haplotype = new Haplotype(gene, name);
          GenePhenotype gp = getPhenotype(gene, source);
          if (gp != null) {
            haplotype.setFunction(gp.getHaplotypeFunction(name));
            haplotype.setActivityValue(gp.getHaplotypeActivity(name));
          }
          haplotype.setReference(name.equals(getReferenceAllele(gene)));
          return haplotype;
        });
  }
}
