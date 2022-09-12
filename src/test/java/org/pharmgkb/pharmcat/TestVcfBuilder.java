package org.pharmgkb.pharmcat;

import java.io.IOException;
import java.io.PrintWriter;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import org.checkerframework.checker.nullness.qual.Nullable;
import org.junit.jupiter.api.TestInfo;
import org.pharmgkb.pharmcat.definition.DefinitionReader;
import org.pharmgkb.pharmcat.definition.model.DefinitionExemption;
import org.pharmgkb.pharmcat.definition.model.DefinitionFile;
import org.pharmgkb.pharmcat.definition.model.VariantLocus;
import org.pharmgkb.pharmcat.util.DataManager;
import org.pharmgkb.pharmcat.util.VcfHelper;


/**
 * Builds test VCF files.
 *
 * @author Mark Woon
 */
public class TestVcfBuilder {
  private final Map<String, Map<String, VcfEdit>> m_edits = new HashMap<>();
  private final Map<String, Map<String, VcfEdit>> m_extraPositions = new HashMap<>();
  private final TestInfo m_testInfo;
  private final String m_name;
  private final List<Path> m_definitionFiles = new ArrayList<>();
  private boolean m_isPhased;
  private boolean m_deleteOnExit = true;


  public TestVcfBuilder(TestInfo testInfo) {
    m_testInfo = testInfo;
    m_name = TestUtils.getTestName(testInfo);
  }

  public TestVcfBuilder(TestInfo testInfo, String name) {
    m_testInfo = testInfo;
    m_name = name;
  }


  public TestVcfBuilder phased() {
    m_isPhased = true;
    return this;
  }

  public TestVcfBuilder saveFile() {
    m_deleteOnExit = false;
    return this;
  }

  /**
   * Specify a specific definition file.
   * If not provided, will use default.
   */
  public TestVcfBuilder withDefinition(Path file) {
    if (!Files.isRegularFile(file)) {
      throw new IllegalArgumentException("Not a file: " + file);
    }
    m_definitionFiles.add(file);
    return this;
  }


  public TestVcfBuilder variation(String gene, String rsid, String... cpicAlleles) {
    VcfEdit edit = new VcfEdit(rsid, cpicAlleles);
    m_edits.computeIfAbsent(gene, g -> new HashMap<>())
        .put(edit.id, edit);
    return this;
  }

  public TestVcfBuilder variation(String gene, String chrom, long position, String... cpicAlleles) {
    VcfEdit edit = new VcfEdit(chrom, position, cpicAlleles);
    m_edits.computeIfAbsent(gene, g -> new HashMap<>())
        .put(edit.id, edit);
    return this;
  }


  public TestVcfBuilder missing(String gene, String... rsids) {
    for (String rsid : rsids) {
      VcfEdit edit = new VcfEdit(rsid, null);
      m_edits.computeIfAbsent(gene, g -> new HashMap<>())
          .put(edit.id, edit);
    }
    return this;
  }

  public TestVcfBuilder missing(String gene, String chrom, long position) {
    VcfEdit edit = new VcfEdit(chrom, position, null);
    m_edits.computeIfAbsent(gene, g -> new HashMap<>())
        .put(edit.id, edit);
    return this;
  }


  public TestVcfBuilder extraPosition(String gene, String rsid, @Nullable String... cpicAlleles) {
    VcfEdit edit = new VcfEdit(rsid, cpicAlleles);
    m_extraPositions.computeIfAbsent(gene, g -> new HashMap<>())
        .put(edit.id, edit);
    return this;
  }


  /**
   * Adds reference for specified gene.
   * This is only necessary to add gene positions to VCF when there is no variation.
   */
  public TestVcfBuilder reference(String gene) {
    m_edits.computeIfAbsent(gene, g -> new HashMap<>());
    return this;
  }


  public Path generate() throws IOException {
    DefinitionReader definitionReader = new DefinitionReader();
    if (m_definitionFiles.size() > 0) {
      for (Path file : m_definitionFiles) {
        definitionReader.read(file);
      }
    } else {
      definitionReader.read(DataManager.DEFAULT_DEFINITION_DIR);
    }
    definitionReader.readExemptions(DataManager.DEFAULT_DEFINITION_DIR.resolve(DataManager.EXEMPTIONS_JSON_FILE_NAME));

    String filename = m_name.replaceAll("\\*", "s")
        .replaceAll("/", "-")
        .replaceAll("[^a-zA-Z0-9_\\-]", "_") + ".vcf";

    Path dir = TestUtils.getTestOutputDir(m_testInfo, false);
    Path file = dir.resolve(filename);
    try (PrintWriter writer = new PrintWriter(Files.newBufferedWriter(file))) {
      VcfHelper.printVcfHeaders(writer, "PharmCAT test (" + m_name + ")",
          VcfHelper.getContigs(m_edits.keySet(), definitionReader));

      for (String gene : m_edits.keySet()) {
        DefinitionFile definitionFile = definitionReader.getDefinitionFile(gene);
        generateForGene(definitionFile, List.of(definitionFile.getVariants()), m_edits.get(gene), writer,
            "PX=" + gene);

        DefinitionExemption exemption = definitionReader.getExemption(gene);
        if (exemption != null && exemption.getExtraPositions() != null) {
          generateForGene(definitionFile, exemption.getExtraPositions(), m_extraPositions.get(gene), writer, "POI");
        }
      }
    }
    if (m_deleteOnExit) {
      file.toFile().deleteOnExit();
    }
    return file;
  }


  private void generateForGene(DefinitionFile definitionFile, Collection<VariantLocus> variants,
      Map<String, VcfEdit> editMap, PrintWriter writer, String info) {

    List<String> matched = new ArrayList<>();
    for (VariantLocus vl : variants) {
      VcfEdit edit = null;
      if (editMap != null) {
        if (vl.getRsid() != null) {
          edit = editMap.get(vl.getRsid());
        } else {
          edit = editMap.get(vl.getChromosome() + ":" + vl.getPosition());
        }
      }
      String sample = m_isPhased ? "0|0" : "0/0";
      if (edit != null) {
        matched.add(edit.id);
        if (edit.cpicAlleles == null || edit.cpicAlleles.length == 0) {
          writer.println("# missing: " + edit.id);
          continue;
        }
        if (edit.cpicAlleles.length == 1) {
          if (!definitionFile.getChromosome().equals("chrX") && !definitionFile.getChromosome().equals("chrY")) {
            throw new IllegalStateException(edit.id + " only has one allele defined and is not on chrX or chrY");
          }
        } else if (edit.cpicAlleles.length > 2) {
          throw new IllegalStateException(edit.id + " has too many alleles (" +
              Arrays.toString(edit.cpicAlleles) + ")");
        }
        StringBuilder builder = new StringBuilder();
        for (String cpicAllele : edit.cpicAlleles) {
          if (builder.length() > 0) {
            builder.append(m_isPhased ? "|" : "/");
          }
          String vcfAllele = vl.getCpicToVcfAlleleMap().get(cpicAllele);
          if (vcfAllele == null) {
            throw new IllegalStateException(edit.id + ": cannot determine VCF allele for " + cpicAllele + ", expected one of " + vl.getCpicToVcfAlleleMap().keySet());
          }
          if (vcfAllele.equals(vl.getRef())) {
            builder.append("0");
          } else {
            int idx = vl.getAlts().indexOf(vcfAllele);
            if (idx == -1) {
              throw new IllegalStateException(edit.id + ":  VCF allele (" + vcfAllele + ") for " + cpicAllele +
                  " is neither ref not alt!");
            }
            builder.append((idx + 1));
          }
        }
        sample = builder.toString();
      }

      VcfHelper.printVcfLine(writer, definitionFile.getChromosome(), vl.getPosition(), vl.getRsid(), vl.getRef(),
          String.join(",", vl.getAlts()), info, sample);
    }

    if (editMap != null) {
      List<String> check = new ArrayList<>(editMap.keySet());
      matched.forEach(check::remove);
      if (check.size() > 0) {
        throw new IllegalStateException("Not involved in named alleles for " + definitionFile.getGeneSymbol() + ": " +
            String.join(", ", check));
      }
    }
  }



  private static class VcfEdit {
    private final String id;
    private final String[] cpicAlleles;

    private VcfEdit(String rsid, @Nullable String[] cpicAlleles) {
      this.id = rsid;
      this.cpicAlleles = cpicAlleles;
    }

    private VcfEdit(String chrom, long posiiton, @Nullable String[] cpicAlleles) {
      this.id = chrom + ":" + posiiton;
      this.cpicAlleles = cpicAlleles;
    }
  }
}
