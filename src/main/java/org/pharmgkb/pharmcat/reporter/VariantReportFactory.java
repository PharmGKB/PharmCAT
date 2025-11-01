package org.pharmgkb.pharmcat.reporter;

import java.util.Comparator;
import java.util.Map;
import java.util.Objects;
import java.util.TreeMap;
import com.google.common.collect.Multimap;
import com.google.common.collect.TreeMultimap;
import org.jspecify.annotations.Nullable;
import org.pharmgkb.pharmcat.Env;
import org.pharmgkb.pharmcat.definition.model.DefinitionFile;
import org.pharmgkb.pharmcat.definition.model.NamedAllele;
import org.pharmgkb.pharmcat.definition.model.VariantLocus;
import org.pharmgkb.pharmcat.haplotype.model.Variant;
import org.pharmgkb.pharmcat.reporter.model.VariantReport;
import org.pharmgkb.pharmcat.util.HaplotypeNameComparator;


/**
 * Factory class to help create new {@link VariantReport} instances for a given gene
 *
 * @author Ryan Whaley
 */
public class VariantReportFactory {

  private final String m_gene;
  private final String m_chr;
  private final Multimap<Long, String> m_variantAlleleMap = TreeMultimap.create(
      Comparator.naturalOrder(),
      HaplotypeNameComparator.getComparator()
  );
  private final Map<Long, String> m_referenceAlleleMap = new TreeMap<>();

  /**
   * Create a new factory for the specified <code>gene</code> (HGNC symbol). This will gather all necessary information
   * from the definition files.
   *
   * @param gene a gene's HGNC symbol
   */
  public VariantReportFactory(String gene, String chr, Env env) {
    m_gene = gene;
    m_chr = chr;

    DefinitionFile definitionFile = env.getDefinitionReader().getDefinitionFile(gene);
    VariantLocus[] allVariants = definitionFile.getVariants();

    for (NamedAllele namedAllele : definitionFile.getNamedAlleles()) {
      if (namedAllele.isReference()) {
        namedAllele.initialize(definitionFile.getVariants());
        for (VariantLocus v : allVariants) {
          m_referenceAlleleMap.put(v.getPosition(), Objects.requireNonNull(namedAllele.getAllele(v)));
        }
      } else {
        @Nullable String[] definingAlleles = namedAllele.getAlleles();
        boolean hasSuballeles = !definitionFile.getSuballelesMap().isEmpty();
        for (int i = 0; i < definingAlleles.length; i++) {
          String alleleValue = definingAlleles[i];
          VariantLocus locus = allVariants[i];

          if (alleleValue != null) {
            String hapName = namedAllele.getName();
            if (hasSuballeles && definitionFile.getSuballelesMap().containsKey(hapName)) {
              hapName = definitionFile.getSuballelesMap().get(hapName);
            }
            m_variantAlleleMap.put(locus.getPosition(), hapName);
          }
        }
      }
    }
  }

  /**
   * Make a new {@link VariantReport} based on information found in the matcher's {@link Variant} class.
   * @param variant a {@link Variant} from the matcher
   * @return a {@link VariantReport} object that has additional data from the given {@link Variant}
   */
  public VariantReport make(Variant variant) {
    return initializeReport(new VariantReport(m_chr, m_gene, variant));
  }

  /**
   * Make a new {@link VariantReport} based on information found in the {@link VariantLocus} object
   * @param locus a {@link VariantLocus} object
   * @return a {@link VariantReport} object that has additional data from the given {@link Variant}
   */
  public VariantReport make(VariantLocus locus) {
    return initializeReport(new VariantReport(m_chr, m_gene, locus));
  }

  /**
   * Fill additional data gathered from the allele definition files
   */
  private VariantReport initializeReport(VariantReport variantReport) {
    variantReport.setReferenceAllele(m_referenceAlleleMap.get(variantReport.getPosition()));
    variantReport.setAlleles(m_variantAlleleMap.get(variantReport.getPosition()));
    return variantReport;
  }
}
