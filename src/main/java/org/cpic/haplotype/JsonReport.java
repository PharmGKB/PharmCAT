package org.cpic.haplotype;

import java.io.BufferedWriter;
import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.Collection;
import java.util.Date;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import javax.annotation.Nonnull;
import com.google.common.base.Preconditions;
import com.google.common.collect.Lists;
import com.google.gson.Gson;
import com.google.gson.GsonBuilder;
import org.cpic.haplotype.model.json.DiplotypeCall;
import org.cpic.haplotype.model.json.HaplotypeCall;
import org.cpic.haplotype.model.json.HaplotyperResult;
import org.cpic.haplotype.model.json.Metadata;


/**
 * @author Mark Woon
 */
public class JsonReport {
  private DefinitionReader m_definitionReader;
  private HaplotyperResult m_root = new HaplotyperResult();
  private Path m_jsonFile;


  public JsonReport(@Nonnull DefinitionReader definitionReader) {
    Preconditions.checkNotNull(definitionReader);
    m_definitionReader = definitionReader;
  }


  public JsonReport forFile(@Nonnull Path vcfFile) {
    Preconditions.checkNotNull(vcfFile);
    Preconditions.checkArgument(vcfFile.toString().endsWith(".vcf"));
    Preconditions.checkArgument(Files.isRegularFile(vcfFile));

    Metadata md = new Metadata();
    //md.setCpicAnnotatorBuild();
    //md.setCpicDataBuild();
    md.setDate(new Date());
    //md.setGenomeAssembly();
    md.setInputFile(vcfFile.getName(vcfFile.getNameCount() - 1).toString());
    m_root.setMetadata(md);

    if (m_jsonFile == null) {
      String baseFilename = vcfFile.getName(vcfFile.getNameCount() - 1).toString();
      baseFilename = baseFilename.substring(0, baseFilename.length() - 4);
      m_jsonFile = vcfFile.getParent().resolve(baseFilename + ".json");
    }
    return this;
  }


  public JsonReport toFile(@Nonnull Path jsonFile) {
    Preconditions.checkNotNull(jsonFile);
    Preconditions.checkArgument(jsonFile.toString().endsWith(".json"));
    Preconditions.checkArgument(Files.isRegularFile(jsonFile));

    m_jsonFile = jsonFile;
    return this;
  }


  public JsonReport haplotype(@Nonnull String gene, @Nonnull List<List<HaplotypeMatch>> matches,
      Collection<SampleAllele> sampleAlleles) {

    Preconditions.checkNotNull(gene);
    Preconditions.checkNotNull(matches);

    DiplotypeCall diplotypeCall = new DiplotypeCall();
    diplotypeCall.setGene(gene);

    // get haplotype/diplotype info
    Map<String, HaplotypeCall> haplotypeMap = new HashMap<>();
    Set<String> diplotypes = new HashSet<>();
    for (List<HaplotypeMatch> pair : matches) {
      HaplotypeMatch hm1 = pair.get(0);
      HaplotypeMatch hm2 = pair.get(1);
      HaplotypeCall hap1 = haplotypeMap.get(hm1.getHaplotype().getCommonName());
      if (hap1 == null) {
        hap1 = new HaplotypeCall(hm1.getHaplotype().getCommonName(), hm1.getSequences());
        haplotypeMap.put(hap1.getName(), hap1);
      }
      HaplotypeCall hap2 = haplotypeMap.get(hm2.getHaplotype().getCommonName());
      if (hap2 == null) {
        hap2 = new HaplotypeCall(hm2.getHaplotype().getCommonName(), hm2.getSequences());
        haplotypeMap.put(hap2.getName(), hap2);
      }

      Lists.newArrayList(hap1, hap2);
    }
    diplotypeCall.setDiplotypes(diplotypes);
    diplotypeCall.setHaplotypes(new HashSet<>(haplotypeMap.values()));

    // get position info
    Map<Integer, SampleAllele> alleleMap = new HashMap<>();
    for (SampleAllele allele : sampleAlleles) {
      alleleMap.put(allele.getPosition(), allele);
    }
    for (Variant variant : m_definitionReader.getHaplotypePositions().get(gene)) {
      SampleAllele allele = alleleMap.get(variant.getPOS());
      String call;
      if (allele.isPhased()) {
        call = allele.getAllele1() + "|" + allele.getAllele2();
      } else {
        call = allele.getAllele1() + "/" + allele.getAllele2();
      }
      diplotypeCall.add(new org.cpic.haplotype.model.json.Variant(variant.getPOS(), variant.get_rsID(), call));
    }


    //diplotypeCall.setHaplotypesNotCalled();
    TSVfile tsvFile = m_definitionReader.getFiles().get(gene);
    diplotypeCall.setGeneVersion(tsvFile.getContentVersion() + " (" + tsvFile.getContentDate() + ")");
    diplotypeCall.setChromosome(tsvFile.getChromosome());
    m_root.addDiplotypeCall(diplotypeCall);

    return this;
  }


  public void print() throws IOException {

    Preconditions.checkState(m_jsonFile != null);
    Gson gson = new GsonBuilder().setPrettyPrinting().create();

    try (BufferedWriter writer = Files.newBufferedWriter(m_jsonFile, StandardCharsets.UTF_8)) {
      writer.write(gson.toJson(m_root));
    }
  }
}
