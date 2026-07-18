package org.pharmgkb.pharmcat.haplotype.benchmark;

import java.io.BufferedReader;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.SortedSet;
import java.util.TreeSet;
import java.util.stream.Collectors;


/**
 * Reads {@code pharmcat_positions.vcf} once and exposes per-gene position data for benchmark VCF construction.
 *
 * <p>Positions are indexed by gene via the {@code PX=<gene>} INFO field. Only the first ALT is retained; positions with
 * an rsid of "." fall back to a {@code chrom:pos} key that {@link org.pharmgkb.pharmcat.TestVcfBuilder} accepts.</p>
 */
public final class PositionsIndex {
  private static final String POSITIONS_VCF_FILENAME = "pharmcat_positions.vcf";
  private static PositionsIndex s_instance;

  private final Map<String, List<Entry>> m_byGene;


  public static synchronized PositionsIndex getInstance() throws IOException {
    if (s_instance == null) {
      s_instance = new PositionsIndex(locatePositionsFile());
    }
    return s_instance;
  }

  private PositionsIndex(Path positionsFile) throws IOException {
    Map<String, List<Entry>> byGene = new HashMap<>();
    try (BufferedReader reader = Files.newBufferedReader(positionsFile)) {
      String line;
      while ((line = reader.readLine()) != null) {
        if (line.isEmpty() || line.startsWith("#")) {
          continue;
        }
        String[] cols = line.split("\t");
        if (cols.length < 8) {
          continue;
        }
        String gene = extractGene(cols[7]);
        if (gene == null) {
          continue;
        }
        String chrom = cols[0];
        long position = Long.parseLong(cols[1]);
        String rsid = cols[2];
        String ref = cols[3];
        String firstAlt = cols[4].split(",", 2)[0];
        byGene.computeIfAbsent(gene, g -> new ArrayList<>())
            .add(new Entry(chrom, position, rsid, ref, firstAlt));
      }
    }
    Map<String, List<Entry>> immutable = new HashMap<>();
    byGene.forEach((g, list) -> immutable.put(g, Collections.unmodifiableList(list)));
    m_byGene = Collections.unmodifiableMap(immutable);
  }

  private static String extractGene(String info) {
    for (String field : info.split(";")) {
      if (field.startsWith("PX=")) {
        return field.substring(3);
      }
    }
    return null;
  }


  List<Entry> forGene(String gene) {
    List<Entry> entries = m_byGene.get(gene);
    if (entries == null) {
      throw new IllegalArgumentException("No positions for " + gene + " in " + POSITIONS_VCF_FILENAME);
    }
    return entries;
  }

  SortedSet<String> genes() {
    return m_byGene.keySet().stream().collect(Collectors.toCollection(TreeSet::new));
  }


  private static Path locatePositionsFile() throws IOException {
    Path start = Path.of("").toAbsolutePath();
    Path cur = start;
    while (cur != null) {
      Path candidate = cur.resolve(POSITIONS_VCF_FILENAME);
      if (Files.isRegularFile(candidate)) {
        return candidate;
      }
      cur = cur.getParent();
    }
    throw new IOException("Could not locate " + POSITIONS_VCF_FILENAME + " from " + start);
  }


  /** Single position record: chrom, position, rsid ("." if absent), REF, first ALT. */
  record Entry(String chrom, long position, String rsid, String ref, String alt) {
    boolean hasRsid() {
      return !rsid.isEmpty() && !rsid.equals(".");
    }
  }
}
