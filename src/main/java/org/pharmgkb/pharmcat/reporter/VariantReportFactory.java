package org.pharmgkb.pharmcat.reporter;

import java.util.Comparator;
import java.util.Map;
import java.util.TreeMap;
import com.google.common.collect.Multimap;
import com.google.common.collect.TreeMultimap;
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
          m_referenceAlleleMap.put(v.getPosition(), namedAllele.getAllele(v));
        }
      } else {
        String[] definingAlleles = namedAllele.getAlleles();
        for (int i = 0; i < definingAlleles.length; i++) {
          String alleleValue = definingAlleles[i];
          VariantLocus locus = allVariants[i];

          if (alleleValue != null) {
            m_variantAlleleMap.put(locus.getPosition(), namedAllele.getName());
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
    return initializeReport(new VariantReport(m_gene, variant));
  }

  /**
   * Make a new {@link VariantReport} based on information found in the {@link VariantLocus} object
   * @param locus a {@link VariantLocus} object
   * @return a {@link VariantReport} object that has additional data from the given {@link Variant}
   */
  public VariantReport make(VariantLocus locus) {
    return initializeReport(new VariantReport(m_gene, locus));
  }

  /**
   * Fill additional data gathered from the allele definition files
   */
  private VariantReport initializeReport(VariantReport variantReport) {
    assignNamedAlleles(variantReport);
    variantReport.setChr(m_chr);
    variantReport.setReferenceAllele(m_referenceAlleleMap.get(variantReport.getPosition()));
    return variantReport;
  }

  /**
   * Fills in the named alleles that the given {@link VariantReport} is used to define
   */
  private void assignNamedAlleles(VariantReport variantReport) {
    variantReport.setAlleles(m_variantAlleleMap.get(variantReport.getPosition()));
  }
}
