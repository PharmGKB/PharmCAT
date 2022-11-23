package org.pharmgkb.pharmcat.haplotype;

import java.io.BufferedReader;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.List;
import com.google.common.base.Preconditions;
import org.pharmgkb.parser.vcf.VcfLineParser;
import org.pharmgkb.parser.vcf.VcfParser;
import org.pharmgkb.parser.vcf.model.ContigMetadata;
import org.pharmgkb.parser.vcf.model.VcfMetadata;
import org.pharmgkb.parser.vcf.model.VcfPosition;
import org.pharmgkb.parser.vcf.model.VcfSample;


/**
 * Simple VCF reader that only retrieves samples.
 *
 * @author Mark Woon
 */
public class VcfSampleReader implements VcfLineParser {
  private final List<String> m_samples = new ArrayList<>();


  public VcfSampleReader(Path vcfFile) throws IOException {
    Preconditions.checkNotNull(vcfFile);
    Preconditions.checkArgument(Files.isRegularFile(vcfFile), "%s is not a file", vcfFile);
    Preconditions.checkArgument(Files.isReadable(vcfFile), "%s is not readable", vcfFile);
    Preconditions.checkArgument(vcfFile.toString().endsWith(".vcf"), "%s is not a VCF file", vcfFile);

    // read VCF file
    try (BufferedReader reader = Files.newBufferedReader(vcfFile);
         VcfParser vcfParser = new VcfParser.Builder()
             .fromReader(reader)
             .parseWith(this)
             .build()) {
      String genomeBuild = null;
      VcfMetadata vcfMetadata = vcfParser.parseMetadata();
      for (ContigMetadata cm : vcfMetadata.getContigs().values()) {
        if (cm.getAssembly() != null) {
          if (genomeBuild == null) {
            genomeBuild = cm.getAssembly();
          } else if (!genomeBuild.equals(cm.getAssembly())) {
            throw new IllegalStateException("VCF file " + vcfFile.getFileName() + " uses different assemblies (" +
                genomeBuild + " vs " + cm.getAssembly() + " for contig)");
          }
        }
      }
      for (int x = 0; x < vcfMetadata.getNumSamples(); x += 1) {
        m_samples.add(vcfMetadata.getSampleName(x));
      }
    }
  }


  public List<String> getSamples() {
    return m_samples;
  }


  @Override
  public void parseLine(VcfMetadata metadata, VcfPosition position, List<VcfSample> sampleData) {
    // don't bother, we're not actually going to read data
  }
}
