package org.pharmgkb.pharmcat.haplotype.benchmark;

import java.util.function.Consumer;


/**
 * Definition of one benchmark scenario: how to build the VCF and how to configure the matcher.
 */
record Scenario(
    String name,
    String gene,
    boolean findCombinations,
    boolean topCandidateOnly,
    boolean callCyp2d6,
    Consumer<BenchmarkVcfBuilder> setup) {
}
