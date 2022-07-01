package org.pharmgkb.pharmcat.reporter.caller;

import java.util.List;
import java.util.Objects;
import java.util.Optional;
import java.util.stream.Collectors;
import java.util.stream.Stream;
import com.google.common.base.Preconditions;
import com.google.common.collect.ImmutableList;
import org.pharmgkb.common.comparator.HaplotypeNameComparator;
import org.pharmgkb.pharmcat.reporter.DiplotypeFactory;
import org.pharmgkb.pharmcat.reporter.model.result.Diplotype;
import org.pharmgkb.pharmcat.reporter.model.result.GeneReport;
import org.pharmgkb.pharmcat.reporter.model.result.Haplotype;
import org.pharmgkb.pharmcat.util.HaplotypeActivityComparator;


public class DpydCustomCaller {
  private static final String DPYD = "DPYD";

  public static boolean shouldBeUsedOn(GeneReport report) {
    return report != null
        && report.getGene().equals(DPYD)
        && !report.isPhased();
  }

  public static Optional<Diplotype> makeLookupCalls(GeneReport report, DiplotypeFactory diplotypeFactory) {
    Preconditions.checkNotNull(report);
    Preconditions.checkArgument(report.getGene().equals(DPYD), "Can only be used on DPYD");

    if (report.isPhased()) {
      return Optional.empty();
    } else {
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
        return Optional.of(newDiplotypes.get(0));
      } else {
        return Optional.empty();
      }
    }
  }
}
