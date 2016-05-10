package org.pharmgkb.pharmcat.haplotype;

import java.lang.invoke.MethodHandles;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.SortedMap;
import java.util.TreeMap;
import java.util.stream.Collectors;
import javax.annotation.Nonnull;
import com.google.common.base.Preconditions;
import org.pharmgkb.pharmcat.definition.model.NamedAllele;
import org.pharmgkb.pharmcat.definition.model.VariantLocus;
import org.pharmgkb.pharmcat.haplotype.model.DiplotypeMatch;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;


/**
 * This is the data used to compute a {@link DiplotypeMatch}.
 *
 * @author Mark Woon
 */
public class MatchData {
  private static final Logger sf_logger = LoggerFactory.getLogger(MethodHandles.lookup().lookupClass());
  SortedMap<Integer, SampleAllele> geneSampleMap = new TreeMap<>();
  Set<VariantLocus> missingPositions = new HashSet<>();
  VariantLocus[] positions;
  List<NamedAllele> haplotypes;
  Set<String> permutations;


  /**
   * Organizes the {@link SampleAllele} data related to the {@code gene} of interest.
   *
   * @param alleleMap map of {@link SampleAllele}s from VCF
   */
  void marshallSampleData(@Nonnull SortedMap<String, SampleAllele> alleleMap, @Nonnull String chromosome,
      @Nonnull VariantLocus[] positions) {

    for (VariantLocus variant : positions) {
      String chrPos = chromosome + ":" + variant.getPosition();
      SampleAllele allele = alleleMap.get(chrPos);
      if (allele == null) {
        missingPositions.add(variant);
        sf_logger.info("Sample has no allele for {}", chrPos);
        continue;
      }
      allele = allele.forVariant(variant);
      geneSampleMap.put(variant.getPosition(), allele);
    }
  }

  /**
   * Organizes the {@link NamedAllele} data for analysis.
   * This will also reorganize haplotypes to deal with samples that have missing alleles.
   */
  void marshallHaplotypes(VariantLocus[] allPositions, List<NamedAllele> allHaplotypes) {

    if (missingPositions.isEmpty()) {
      positions = allPositions;
      haplotypes = allHaplotypes;

    } else {
      // handle missing positions by duplicating haplotype and eliminating missing positions
      VariantLocus[] availablePositions = new VariantLocus[allPositions.length - missingPositions.size()];
      for (int x = 0, y = 0; x < allPositions.length; x += 1) {
        if (!missingPositions.contains(allPositions[x])) {
          availablePositions[y] = allPositions[x];
          y += 1;
        }
      }
      positions = availablePositions;
      haplotypes = new ArrayList<>();
      for (NamedAllele hap : allHaplotypes) {
        // get alleles for positions we have data on
        String[] availableAlleles = new String[positions.length];
        for (int x = 0; x < positions.length; x += 1) {
          availableAlleles[x] = hap.getAllele(positions[x]);
        }
        NamedAllele newHap = new NamedAllele(hap.getId(), hap.getName(), availableAlleles);
        newHap.setFunction(hap.getFunction());
        newHap.setPopFreqMap(hap.getPopFreqMap());
        newHap.finalize(positions);
        if (newHap.getNumValidAlleles() > 0) {
          haplotypes.add(newHap);
        }
      }
    }
  }


  /**
   * Assumes that missing alleles in {@link NamedAllele}s should be the reference.
   */
  void defaultMissingAllelesToReference() {

    List<NamedAllele> updatedHaplotypes = new ArrayList<>();
    NamedAllele referenceHaplotype = null;
    for (NamedAllele hap : haplotypes) {
      if (referenceHaplotype == null) {
        referenceHaplotype = hap;
        updatedHaplotypes.add(hap);
        continue;
      }

      String[] refAlleles = referenceHaplotype.getAlleles();
      String[] curAlleles = hap.getAlleles();
      Preconditions.checkState(refAlleles.length == curAlleles.length);

      String[] newAlleles = new String[refAlleles.length];
      for (int x = 0; x < refAlleles.length; x += 1) {
        if (curAlleles[x] == null) {
          newAlleles[x] = refAlleles[x];
        } else {
          newAlleles[x] = curAlleles[x];
        }
      }

      NamedAllele fixedHap = new NamedAllele(hap.getId(), hap.getName(), newAlleles);
      fixedHap.setFunction(hap.getFunction());
      fixedHap.setPopFreqMap(hap.getPopFreqMap());
      fixedHap.finalize(positions);
      updatedHaplotypes.add(fixedHap);
    }

    haplotypes = updatedHaplotypes;
  }


  /**
   * Generate all permutations of sample at positions of interest.
   */
  void generateSamplePermutations() {

    permutations = CombinationUtil.generatePermutations(
        geneSampleMap.values().stream()
            .sorted()
            .collect(Collectors.toList())
    );
  }
}
