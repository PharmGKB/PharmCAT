package org.pharmgkb.pharmcat.haplotype.benchmark;

import java.io.IOException;
import java.io.PrintWriter;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.HashMap;
import java.util.List;
import org.jspecify.annotations.Nullable;
import org.pharmgkb.pharmcat.definition.DefinitionReader;
import org.pharmgkb.pharmcat.definition.model.DefinitionFile;
import org.pharmgkb.pharmcat.util.VcfHelper;


/**
 * Writes a per-gene benchmark VCF from {@link PositionsIndex} entries.
 *
 * <p>Every gene position gets emitted; the caller sets genotypes and optional phase sets by index or rsid. Anything
 * left unset defaults to {@code 0/0} (or {@code 0|0} when {@link #setDefaultPhased(boolean)} is on).</p>
 */
public final class BenchmarkVcfBuilder {
  private final String m_gene;
  private final List<PositionsIndex.Entry> m_entries;
  private final DefinitionFile m_definitionFile;
  private final GtSpec[] m_specs;
  private boolean m_defaultPhased;


  public BenchmarkVcfBuilder(String gene, PositionsIndex positionsIndex, DefinitionReader definitionReader) {
    m_gene = gene;
    m_entries = positionsIndex.forGene(gene);
    m_definitionFile = definitionReader.getDefinitionFile(gene);
    m_specs = new GtSpec[m_entries.size()];
  }


  public int size() {
    return m_entries.size();
  }

  PositionsIndex.Entry entry(int index) {
    return m_entries.get(index);
  }

  int indexOfRsid(String rsid) {
    for (int i = 0; i < m_entries.size(); i += 1) {
      if (rsid.equals(m_entries.get(i).rsid())) {
        return i;
      }
    }
    throw new IllegalArgumentException("No position with rsid " + rsid + " for " + m_gene);
  }


  BenchmarkVcfBuilder setDefaultPhased(boolean phased) {
    m_defaultPhased = phased;
    return this;
  }

  public BenchmarkVcfBuilder set(int index, String gt) {
    m_specs[index] = new GtSpec(gt, null);
    return this;
  }

  BenchmarkVcfBuilder set(int index, String gt, int phaseSet) {
    m_specs[index] = new GtSpec(gt, phaseSet);
    return this;
  }

  BenchmarkVcfBuilder setRsid(String rsid, String gt) {
    return set(indexOfRsid(rsid), gt);
  }

  BenchmarkVcfBuilder setRsid(String rsid, String gt, int phaseSet) {
    return set(indexOfRsid(rsid), gt, phaseSet);
  }


  public Path write(Path outFile) throws IOException {
    HashMap<String, String> contigs = new HashMap<>();
    contigs.put(m_definitionFile.getChromosome(), m_definitionFile.getGenomeBuild());
    try (PrintWriter writer = new PrintWriter(Files.newBufferedWriter(outFile))) {
      VcfHelper.printVcfHeaders(writer, "PharmCAT benchmark (" + m_gene + ")", contigs);
      String defaultGt = m_defaultPhased ? "0|0" : "0/0";
      for (int i = 0; i < m_entries.size(); i += 1) {
        PositionsIndex.Entry e = m_entries.get(i);
        GtSpec spec = m_specs[i];
        String gt = spec == null ? defaultGt : spec.gt();
        Integer ps = spec == null ? null : spec.phaseSet();
        String rsid = e.hasRsid() ? e.rsid() : null;
        VcfHelper.printVcfLine(writer, e.chrom(), e.position(), rsid, e.ref(), e.alt(),
            "PX=" + m_gene, gt, ps);
      }
    }
    return outFile;
  }


  private record GtSpec(String gt, @Nullable Integer phaseSet) {}
}
