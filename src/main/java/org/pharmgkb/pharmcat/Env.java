package org.pharmgkb.pharmcat;

import java.io.IOException;
import java.nio.file.Path;
import java.util.HashMap;
import java.util.Map;
import org.checkerframework.checker.nullness.qual.Nullable;
import org.pharmgkb.pharmcat.definition.DefinitionReader;
import org.pharmgkb.pharmcat.definition.PhenotypeMap;
import org.pharmgkb.pharmcat.definition.model.GenePhenotype;
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
