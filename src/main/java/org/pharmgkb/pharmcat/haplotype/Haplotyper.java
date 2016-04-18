package org.pharmgkb.pharmcat.haplotype;

import java.io.IOException;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.SortedSet;
import javax.annotation.Nonnull;
import javax.annotation.Nullable;
import com.google.common.base.Preconditions;
import com.google.common.collect.Lists;


/**
 * @author Mark Woon
 */
public class Haplotyper {
  private DefinitionReader m_definitionReader;
  private VcfReader m_vcfReader;

  public Haplotyper(@Nonnull DefinitionReader definitionReader) {

    Preconditions.checkNotNull(definitionReader);
    m_definitionReader = definitionReader;
    m_vcfReader = new VcfReader(calculateLocationsOfInterest(m_definitionReader));
  }

  protected static Set<String> calculateLocationsOfInterest(DefinitionReader definitionReader) {

    Set<String> locationsOfInterest = new HashSet<>();
    for (String gene : definitionReader.getHaplotypes().keySet()) {
      List<Variant> variants = definitionReader.getHaplotypePositions().get(gene);
      // map variant to chr:position
      for (Variant v : variants) {
        String chrPos = v.getCHROM() + ":" + v.getPOS();
        locationsOfInterest.add(chrPos);
      }
      // calculate permutations of all haplotypes
      for (Haplotype hap : definitionReader.getHaplotypes().get(gene)) {
        hap.calculatePermutations(variants);
      }
    }
    return locationsOfInterest;
  }


  public void call(@Nonnull Path vcfFile, @Nullable Path jsonFile) throws IOException {

    Map<String, SampleAllele> alleles = m_vcfReader.read(vcfFile);
    JsonReport report = new JsonReport(m_definitionReader)
        .forFile(vcfFile);
    if (jsonFile != null) {
      report.toFile(jsonFile);
    }
    // call haplotypes
    for (String gene : m_definitionReader.getHaplotypePositions().keySet()) {
      report.haplotype(gene, callHaplotype(alleles, gene), alleles.values());
    }
    report.print();
  }


  public List<List<HaplotypeMatch>> callHaplotype(Map<String, SampleAllele> alleleMap, String gene) {

    List<Variant> variants = m_definitionReader.getHaplotypePositions().get(gene);
    List<Haplotype> haplotypes = m_definitionReader.getHaplotypes().get(gene);

    List<SampleAllele> alleles = new ArrayList<>();
    for (Variant variant : variants) {
      String chrPos = variant.getCHROM() + ":" + variant.getPOS();
      SampleAllele allele = alleleMap.get(chrPos);
      if (allele == null) {
        throw new RuntimeException("Sample has no allele for " + chrPos + " (ref is " + variant.getREF() + ")");
      }
      alleles.add(allele);
    }

    // get sample permutations
    Set<String> permutations = CombinationUtil.generatePermutations(alleles);
    // compare sample permutations to haplotypes
    SortedSet<HaplotypeMatch> matches = ComparisonUtil.comparePermutations(permutations, haplotypes);

    if (permutations.size() == 1 && matches.size() == 1) {
      // sample is homozygous for all positions and it matches a single allele,
      // so we need to return that as a diplotype
      HaplotypeMatch hm = matches.first();
      List<List<HaplotypeMatch>> diplotype = new ArrayList<>();
      diplotype.add(Lists.newArrayList(hm, hm));
      return diplotype;
    }

    // find pair-wise matches
    List<List<Haplotype>> pairs = CombinationUtil.generatePerfectPairs(haplotypes);
    return ComparisonUtil.determinePairs(matches, pairs);
  }
}
