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
import org.pharmgkb.pharmcat.haplotype.model.Metadata;
import org.pharmgkb.pharmcat.haplotype.model.Result;
import org.pharmgkb.pharmcat.haplotype.model.Variant;


/**
 * This class collects results for the {@link NamedAlleleMatcher}.
 *
 * @author Mark Woon
 */
public class ResultBuilder {
  private DefinitionReader m_definitionReader;
  private Result m_result = new Result();
  private SimpleDateFormat m_dateFormat = new SimpleDateFormat("MM/dd/yy");


  public ResultBuilder(DefinitionReader definitionReader) {
    Preconditions.checkNotNull(definitionReader);
    m_definitionReader = definitionReader;
  }

  public Result build() {
    return m_result;
  }


  public ResultBuilder forFile(Path vcfFile, Map<String, Collection<String>> warnings) {
    Preconditions.checkNotNull(vcfFile);
    Preconditions.checkArgument(vcfFile.toString().endsWith(".vcf"));
    Preconditions.checkArgument(Files.isRegularFile(vcfFile));

    m_result.setMetadata(new Metadata(NamedAlleleMatcher.VERSION, m_definitionReader.getGenomeBuild(),
        PathUtils.getFilename(vcfFile), new Date()));
    if (warnings != null) {
      m_result.setVcfWarnings(warnings);
    }
    return this;
  }


  protected ResultBuilder gene(String gene, MatchData matchData, List<DiplotypeMatch> matches) {
    Preconditions.checkNotNull(gene);

    DefinitionFile tsvFile = m_definitionReader.getDefinitionFile(gene);
    String definitionVersion = m_dateFormat.format(tsvFile.getModificationDate());
    String chromosome = tsvFile.getChromosome();

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
          .filter(h -> !exemption.shouldIgnore(h))
          .collect(Collectors.toSet());
      ignoredHaplotypes = exemption.getIgnoredAlleles().stream()
          .map(String::toUpperCase)
          .collect(Collectors.toSet());
    } else {
       ignoredHaplotypes = new HashSet<>();
    }

    GeneCall geneCall = new GeneCall(definitionVersion, chromosome, gene, matchData, uncallableHaplotypes,
        ignoredHaplotypes);
    if (matches != null) {
      // get haplotype/diplotype info
      for (DiplotypeMatch dm : matches) {
        geneCall.addDiplotype(dm);
      }
    }

    // get position info
    for (VariantLocus variant : matchData.getPositions()) {
      geneCall.add(new Variant(variant, matchData.getSampleAllele(variant.getVcfPosition())));
    }

    m_result.addDiplotypeCall(geneCall);

    return this;
  }
}
