package org.pharmgkb.pharmcat.reporter.caller;

import java.util.Arrays;
import java.util.List;
import java.util.Objects;
import java.util.Optional;
import java.util.stream.Collectors;
import java.util.stream.Stream;
import com.google.common.base.Preconditions;
import com.google.common.collect.ImmutableList;
import org.pharmgkb.common.comparator.HaplotypeNameComparator;
import org.pharmgkb.pharmcat.haplotype.model.GeneCall;
import org.pharmgkb.pharmcat.reporter.DiplotypeFactory;
import org.pharmgkb.pharmcat.reporter.model.result.Diplotype;
import org.pharmgkb.pharmcat.reporter.model.result.GeneReport;
import org.pharmgkb.pharmcat.reporter.model.result.Haplotype;
import org.pharmgkb.pharmcat.reporter.model.result.Observation;
import org.pharmgkb.pharmcat.util.HaplotypeActivityComparator;


public class DpydCustomCaller {
  private static final String DPYD = "DPYD";

  public static boolean shouldBeUsedOn(GeneReport report) {
    return report != null && report.getGene().equals(DPYD);
  }

  public static Optional<Diplotype> makeLookupCalls(GeneReport report, DiplotypeFactory diplotypeFactory, GeneCall geneCall) {
    Preconditions.checkNotNull(report);
    Preconditions.checkArgument(report.getGene().equals(DPYD), "Can only be used on DPYD");

    if (report.isPhased()) {
      if (geneCall.getDiplotypes().size() > 0) {
        Optional<Haplotype> hap1 = diplotypeFactory.makeLeastFunctionHaplotypeByName(geneCall.getDiplotypes().stream()
            .map(d -> d.getHaplotype1().getHaplotype().getName())
            .flatMap(n -> Arrays.stream(n.split(" \\+ ")))
            .collect(Collectors.toList()));
        Optional<Haplotype> hap2 = diplotypeFactory.makeLeastFunctionHaplotypeByName(geneCall.getDiplotypes().stream()
            .map(d -> d.getHaplotype2().getHaplotype().getName())
            .flatMap(n -> Arrays.stream(n.split(" \\+ ")))
            .collect(Collectors.toList()));
        if (hap1.isPresent() && hap2.isPresent()) {
          List<String> diplotypeStrings = ImmutableList.of(hap1.get().getName() + "/" + hap2.get().getName());
          List<Diplotype> diplotypes = diplotypeFactory.makeDiplotypes(diplotypeStrings);
          return diplotypes.stream().findFirst();
        }
        else {
          return Optional.empty();
        }
      }
      else {
        return Optional.empty();
      }
    }
    // if unphased
    else {
      List<Haplotype> haplotypes = report.getMatcherDiplotypes().stream()
          .flatMap(d -> Stream.of(d.getAllele1(), d.getAllele2()))
          .filter(Objects::nonNull)
          .sorted(HaplotypeActivityComparator.getComparator())
          .toList();

      if (haplotypes.size() >= 2) {
        String diplotypeName = haplotypes.subList(0, 2).stream()
            .map(Haplotype::getName)
            .sorted(HaplotypeNameComparator.getComparator())
            .collect(Collectors.joining("/"));
        List<Diplotype> newDiplotypes = diplotypeFactory.makeDiplotypes(ImmutableList.of(diplotypeName));
        newDiplotypes.forEach(d -> d.setObserved(Observation.INFERRED));
        return Optional.of(newDiplotypes.get(0));
      } else {
        return Optional.empty();
      }
    }
  }
}
