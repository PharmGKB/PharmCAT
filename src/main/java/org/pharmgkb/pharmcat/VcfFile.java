package org.pharmgkb.pharmcat;

import java.io.BufferedReader;
import java.io.ByteArrayInputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.List;
import java.util.zip.GZIPInputStream;
import org.checkerframework.checker.nullness.qual.Nullable;
import org.pharmgkb.pharmcat.definition.DefinitionReader;
import org.pharmgkb.pharmcat.haplotype.VcfReader;
import org.pharmgkb.pharmcat.haplotype.VcfSampleReader;


/**
 * This class supports working with VCF files.
 * If file size is less than 20% of max available memory, it will be read into memory.
 *
 * @author Mark Woon
 */
public class VcfFile {
  private final Path m_vcfFile;
  private final boolean m_isGzipped;
  private final boolean m_readIntoMemory;
  private byte[] m_data;
  private List<String> m_samples = new ArrayList<>();


  public VcfFile(Path vcfFile) throws ReportableException, IOException {
    if (!isVcfFile(vcfFile)) {
      throw new ReportableException(vcfFile + " is not a VCF file");
    }
    m_vcfFile = vcfFile;
    m_isGzipped = isGzippedVcfFile(vcfFile);
    long maxMem = Runtime.getRuntime().maxMemory();
    m_readIntoMemory = Files.size(vcfFile) < (maxMem / 5);
  }

  public VcfFile(Path vcfFile, boolean readIntoMemory) throws ReportableException, IOException {
    if (!isVcfFile(vcfFile)) {
      throw new ReportableException(vcfFile + " is not a VCF file");
    }
    m_vcfFile = vcfFile;
    m_isGzipped = isGzippedVcfFile(vcfFile);
    m_readIntoMemory = readIntoMemory;
  }


  private BufferedReader open() throws IOException {
    if (m_readIntoMemory) {
      if (m_data == null) {
        m_data = Files.readAllBytes(m_vcfFile);
      }
      if (m_isGzipped) {
        return new BufferedReader(new InputStreamReader(new GZIPInputStream(new ByteArrayInputStream(m_data))));
      }
      return new BufferedReader(new InputStreamReader(new ByteArrayInputStream(m_data)));

    } else {
      if (m_isGzipped) {
        return new BufferedReader(new InputStreamReader(new GZIPInputStream(Files.newInputStream(m_vcfFile))));
      }
      return Files.newBufferedReader(m_vcfFile);
    }
  }


  public Path getFile() {
    return m_vcfFile;
  }

  public List<String> getSamples() throws IOException {
    if (m_samples.size() == 0) {
      try (BufferedReader reader = open()) {
        VcfSampleReader vcfSampleReader = new VcfSampleReader(reader);
        m_samples = vcfSampleReader.getSamples();
      }
    }
    return m_samples;
  }


  public VcfReader getReader(DefinitionReader definitionReader, @Nullable String sampleId)
      throws IOException {
    try (BufferedReader reader = open()) {
      return new VcfReader(definitionReader, reader, sampleId);
    }
  }


  public static boolean isGzippedVcfFile(Path vcfFile) {
    String filename = vcfFile.toString();
    return filename.endsWith(".vcf.bgz") || filename.endsWith(".vcf.gz");
  }

  public static boolean isVcfFile(Path vcfFile) {
    if (!Files.isRegularFile(vcfFile)) {
      return false;
    }
    String filename = vcfFile.toString();
    return filename.endsWith(".vcf") || isGzippedVcfFile(vcfFile);
  }
}
