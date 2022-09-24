package org.pharmgkb.pharmcat;

import java.io.IOException;
import java.nio.file.Path;
import java.util.HashMap;
import java.util.Map;
import java.util.SortedSet;
import org.checkerframework.checker.nullness.qual.Nullable;
import org.pharmgkb.pharmcat.definition.DefinitionReader;
import org.pharmgkb.pharmcat.phenotype.PhenotypeMap;
import org.pharmgkb.pharmcat.phenotype.model.GenePhenotype;
import org.pharmgkb.pharmcat.reporter.DrugCollection;
import org.pharmgkb.pharmcat.reporter.MessageHelper;
import org.pharmgkb.pharmcat.reporter.PgkbGuidelineCollection;
import org.pharmgkb.pharmcat.reporter.model.DataSource;
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
  private final DrugCollection m_cpicDrugs;
  private final PgkbGuidelineCollection m_dpwgDrugs;
  private MessageHelper m_messageHelper;
  private final Map<DataSource, Map<String, Map<String, Haplotype>>> m_haplotypeCache = new HashMap<>();


  public Env() throws IOException, ReportableException {
    this(null);
  }

  public Env(@Nullable Path definitionDir) throws IOException, ReportableException {
    m_definitionReader = new DefinitionReader();
    if (definitionDir != null) {
      m_definitionReader.read(definitionDir);
      if (m_definitionReader.getGenes().size() == 0) {
        throw new ReportableException("Did not find any allele definitions at " + definitionDir);
      }
    } else {
      m_definitionReader.read(DataManager.DEFAULT_DEFINITION_DIR);
    }
    m_phenotypeMap = new PhenotypeMap();
    m_cpicDrugs = new DrugCollection();
    m_dpwgDrugs = new PgkbGuidelineCollection();
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


  public DrugCollection getCpicDrugs() {
    return m_cpicDrugs;
  }

  public PgkbGuidelineCollection getDpwgDrugs() {
    return m_dpwgDrugs;
  }


  public SortedSet<String> getCpicGenes() {
    return m_cpicDrugs.getAllReportableGenes();
  }

  public SortedSet<String> getDpwgGenes() {
    return m_dpwgDrugs.getGenes();
  }

  public boolean hasGene(DataSource source, String gene) {
    if (source == DataSource.CPIC) {
      return getCpicGenes().contains(gene);
    } else {
      return getDpwgGenes().contains(gene);
    }
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


  /**
   * Make or retrieve a cached {@link Haplotype} object that corresponds to the given allele name.
   */
  public Haplotype makeHaplotype(String gene, String name, DataSource source) {
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
