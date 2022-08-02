package org.pharmgkb.pharmcat.haplotype;

import java.nio.file.Files;
import java.nio.file.Path;
import java.text.SimpleDateFormat;
import java.util.Collection;
import java.util.Date;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;
import com.google.common.base.Preconditions;
import org.pharmgkb.common.util.PathUtils;
import org.pharmgkb.pharmcat.definition.model.DefinitionExemption;
import org.pharmgkb.pharmcat.definition.model.DefinitionFile;
import org.pharmgkb.pharmcat.definition.model.NamedAllele;
import org.pharmgkb.pharmcat.definition.model.VariantLocus;
import org.pharmgkb.pharmcat.haplotype.model.DiplotypeMatch;
import org.pharmgkb.pharmcat.haplotype.model.GeneCall;
import org.pharmgkb.pharmcat.haplotype.model.HaplotypeMatch;
import org.pharmgkb.pharmcat.haplotype.model.Metadata;
import org.pharmgkb.pharmcat.haplotype.model.Result;
import org.pharmgkb.pharmcat.haplotype.model.Variant;


/**
 * This class collects results for the {@link NamedAlleleMatcher}.
 *
 * @author Mark Woon
 */
public class ResultBuilder {
  private final DefinitionReader m_definitionReader;
  private final Result m_result = new Result();
  private final SimpleDateFormat m_dateFormat = new SimpleDateFormat("MM/dd/yy");
  private final boolean m_topCandidatesOnly;
  private final boolean m_findCombinations;
  private final boolean m_callCyp2d6;


  public ResultBuilder(DefinitionReader definitionReader, boolean topCandidatesOnly, boolean findCombinations,
      boolean callCyp2d6) {
    Preconditions.checkNotNull(definitionReader);
    m_definitionReader = definitionReader;
    m_topCandidatesOnly = topCandidatesOnly;
    m_findCombinations = findCombinations;
    m_callCyp2d6 = callCyp2d6;
  }

  public Result build() {
    return m_result;
  }


  public ResultBuilder forFile(Path vcfFile, Map<String, Collection<String>> warnings) {
    Preconditions.checkNotNull(vcfFile);
    Preconditions.checkArgument(vcfFile.toString().endsWith(".vcf"));
    Preconditions.checkArgument(Files.isRegularFile(vcfFile));

    m_result.setMetadata(new Metadata(NamedAlleleMatcher.VERSION, m_definitionReader.getGenomeBuild(),
        PathUtils.getFilename(vcfFile), new Date(), m_topCandidatesOnly, m_findCombinations, m_callCyp2d6));
    if (warnings != null) {
      m_result.setVcfWarnings(warnings);
    }
    return this;
  }


  /**
   * Builds result for gene when VCF has with no samples for it.
   */
  protected ResultBuilder gene(String gene, MatchData matchData) {
    Preconditions.checkNotNull(gene);
    m_result.addGeneCall(initGeneCall(gene, matchData));
    return this;
  }


  /**
   * Adds diplotype results for specified gene.
   */
  protected ResultBuilder gene(String gene, MatchData matchData, List<DiplotypeMatch> matches) {
    Preconditions.checkNotNull(gene);

    GeneCall geneCall = initGeneCall(gene, matchData);
    // get haplotype/diplotype info
    for (DiplotypeMatch dm : matches) {
      geneCall.addDiplotype(dm);
    }

    m_result.addGeneCall(geneCall);
    return this;
  }


  /**
   * Add haplotype results for specified gene.
   * This should only be used for genes that we only want to get haplotype matches for (e.g. DPYD).
   */
  protected ResultBuilder gene(String gene, MatchData matchData, Set<HaplotypeMatch> matches) {
    Preconditions.checkNotNull(gene);

    GeneCall geneCall = initGeneCall(gene, matchData);
    if (matches != null) {
      geneCall.addAllHaplotypes(matches);
    }

    m_result.addGeneCall(geneCall);
    return this;
  }


  private GeneCall initGeneCall(String gene, MatchData matchData) {
    DefinitionFile definitionFile = m_definitionReader.getDefinitionFile(gene);
    String definitionVersion = m_dateFormat.format(definitionFile.getModificationDate());
    String chromosome = definitionFile.getChromosome();

    Set<String> matchableHaps = matchData.getHaplotypes().stream()
        .map(NamedAllele::getName)
        .collect(Collectors.toSet());
    Set<String> uncallableHaplotypes = m_definitionReader.getHaplotypes(gene).stream()
        .map(NamedAllele::getName)
        .filter(n -> !matchableHaps.contains(n))
        .collect(Collectors.toSet());
    Set<String> ignoredHaplotypes;
    DefinitionExemption exemption = m_definitionReader.getExemption(gene);
    if (exemption != null) {
      uncallableHaplotypes = uncallableHaplotypes.stream()
          .filter(h -> !exemption.shouldIgnoreAllele(h))
          .collect(Collectors.toSet());
      ignoredHaplotypes = exemption.getIgnoredAlleles().stream()
          .map(String::toUpperCase)
          .collect(Collectors.toSet());
    } else {
      ignoredHaplotypes = new HashSet<>();
    }

    GeneCall geneCall = new GeneCall(definitionVersion, chromosome, gene, matchData, uncallableHaplotypes,
        ignoredHaplotypes);

    // get position info
    for (VariantLocus variant : matchData.getPositions()) {
      geneCall.add(new Variant(variant, matchData.getSampleAllele(variant.getPosition())));
    }
    geneCall.finalizeVariants();

    return geneCall;
  }
}
