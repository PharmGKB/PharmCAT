package org.cpic.haplotype;

import java.io.IOException;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.SortedSet;
import javax.annotation.Nonnull;
import com.google.common.base.Preconditions;


/**
 * @author Mark Woon
 */
public class Haplotyper {
  private DefinitionReader m_definitionReader;
  private VcfReader m_vcfReader;

  public Haplotyper(@Nonnull DefinitionReader definitionReader) {

    Preconditions.checkNotNull(definitionReader);
    m_definitionReader = definitionReader;

    Set<String> locationsOfInterest = new HashSet<>();
    for (String gene : m_definitionReader.getHaplotypes().keySet()) {
      List<Variant> variants = m_definitionReader.getHaplotypePositions().get(gene);
      // map variant to chr:position
      for (Variant v : variants) {
        String chrPos = v.getCHROM() + ":" + v.getPOS();
        locationsOfInterest.add(chrPos);
      }
      // calculate permutations of all haplotypes
      for (Haplotype hap : m_definitionReader.getHaplotypes().get(gene)) {
        hap.calculatePermutations(variants);
      }
    }
    m_vcfReader = new VcfReader(locationsOfInterest);
  }


  public void call(@Nonnull Path vcfFile) throws IOException {

    Map<String, SampleAllele> alleles = m_vcfReader.read(vcfFile);
    // call haplotypes
    for (String gene : m_definitionReader.getHaplotypePositions().keySet()) {
      callHaplotype(alleles, gene);
    }
  }


  public void callHaplotype(Map<String, SampleAllele> alleleMap, String gene) {

    List<Variant> variants = m_definitionReader.getHaplotypePositions().get(gene);
    List<Haplotype> haplotypes = m_definitionReader.getHaplotypes().get(gene);

    List<SampleAllele> alleles = new ArrayList<>();
    for (Variant variant : variants) {
      String chrPos = variant.getCHROM() + ":" + variant.getPOS();
      SampleAllele allele = alleleMap.get(chrPos);
      if (allele == null) {
        throw new RuntimeException("Sample has no allele for " + chrPos);
      }
      alleles.add(allele);
    }

    // get sample permutations
    Set<String> permutations = CombinationUtil.generatePermutations(alleles);
    // compare sample permutations to haplotypes
    SortedSet<HaplotypeMatch> matches = ComparisonUtil.comparePermutations(permutations, haplotypes);

    // find pair-wise matches
    List<List<Haplotype>> pairs = CombinationUtil.generatePerfectPairs(haplotypes);
    ComparisonUtil.determinePairs(matches, pairs);
  }
}
