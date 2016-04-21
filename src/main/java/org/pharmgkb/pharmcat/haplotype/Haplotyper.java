package org.pharmgkb.pharmcat.haplotype;

import java.io.IOException;
import java.nio.file.Path;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.SortedMap;
import java.util.SortedSet;
import java.util.TreeMap;
import java.util.stream.Collectors;
import javax.annotation.Nonnull;
import com.google.common.base.Preconditions;
import com.google.common.collect.Lists;
import org.pharmgkb.pharmcat.haplotype.model.DiplotypeMatch;
import org.pharmgkb.pharmcat.haplotype.model.HaplotypeMatch;


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


  public Report call(@Nonnull Path vcfFile) throws IOException {

    SortedMap<String, SampleAllele> alleles = m_vcfReader.read(vcfFile);
    Report report = new Report(m_definitionReader)
        .forFile(vcfFile);
    // call haplotypes
    for (String gene : m_definitionReader.getHaplotypePositions().keySet()) {
      report.gene(gene, callHaplotype(alleles, gene), alleles.values());
    }
    return report;
  }


  public List<DiplotypeMatch> callHaplotype(SortedMap<String, SampleAllele> alleleMap, String gene) {

    List<Variant> variants = m_definitionReader.getHaplotypePositions().get(gene);
    List<Haplotype> haplotypes = m_definitionReader.getHaplotypes().get(gene);

    SortedMap<Integer, SampleAllele> haplotypeSampleMap = new TreeMap<>();
    for (Variant variant : variants) {
      String chrPos = variant.getCHROM() + ":" + variant.getPOS();
      SampleAllele allele = alleleMap.get(chrPos);
      if (allele == null) {
        throw new RuntimeException("Sample has no allele for " + chrPos + " (ref is " + variant.getREF() + ")");
      }
      haplotypeSampleMap.put(variant.getPOS(), allele);
    }

    // get sample permutations
    Set<String> permutations = CombinationUtil.generatePermutations(
        haplotypeSampleMap.values().stream()
            .sorted()
            .collect(Collectors.toList()));
    // compare sample permutations to haplotypes
    SortedSet<HaplotypeMatch> matches = ComparisonUtil.comparePermutations(permutations, haplotypes);

    if (permutations.size() == 1 && matches.size() == 1) {
      // sample is homozygous for all positions and it matches a single allele,
      // so we need to return that as a diplotype
      HaplotypeMatch hm = matches.first();
      DiplotypeMatch dm = new DiplotypeMatch(hm, hm);
      String seq = permutations.iterator().next();
      dm.addSequencePair(new String[] { seq, seq });
      return Lists.newArrayList(dm);
    }

    // find matched pairs
    return ComparisonUtil.determinePairs(haplotypeSampleMap, matches);
  }
}
