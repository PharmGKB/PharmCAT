package org.pharmgkb.pharmcat.reporter;

import java.io.IOException;
import java.util.Collection;
import java.util.Comparator;
import java.util.Iterator;
import java.util.Map;
import java.util.TreeMap;
import com.google.common.collect.Multimap;
import com.google.common.collect.TreeMultimap;
import org.pharmgkb.common.comparator.HaplotypeNameComparator;
import org.pharmgkb.pharmcat.Env;
import org.pharmgkb.pharmcat.definition.model.DefinitionFile;
import org.pharmgkb.pharmcat.definition.model.NamedAllele;
import org.pharmgkb.pharmcat.definition.model.VariantLocus;
import org.pharmgkb.pharmcat.haplotype.model.Variant;
import org.pharmgkb.pharmcat.reporter.model.VariantReport;


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
  private final Map<Long, String> m_wildAlleleMap = new TreeMap<>();

  /**
   * Create a new factory for the specified <code>gene</code> (HGNC symbol). This will gather all necessary information
   * from the definition files.
   * @param gene a gene's HGNC symbol
   * @throws IOException can occur from writing the JSON file
   */
  public VariantReportFactory(String gene, String chr, Env env) throws IOException {
    m_gene = gene;
    m_chr = chr;

    DefinitionFile definitionFile = env.getDefinitionReader().getDefinitionFile(gene);
    VariantLocus[] allVariants = definitionFile.getVariants();
    Iterator<NamedAllele> it = definitionFile.getNamedAlleles().iterator();

    NamedAllele wildNamedAllele = it.next();
    wildNamedAllele.initialize(definitionFile.getVariants());
    for (VariantLocus v : allVariants) {
      m_wildAlleleMap.put(v.getPosition(), wildNamedAllele.getAllele(v));
    }

    // purposely skip the first allele since it's "default"
    while (it.hasNext()) {
      NamedAllele namedAllele = it.next();

      String[] definingAlleles = namedAllele.getAlleles();
      for (int i = 0; i < definingAlleles.length; i++) {
        String alleleValue = definingAlleles[i];
        VariantLocus locus = allVariants[i];

        if (alleleValue != null) {
          // this is a special case for CYP2C19 since *1 is not the "reference" allele
          // this is here to prevent *1 from being shown repetitively in reports
          if (!gene.equals("CYP2C19") || !namedAllele.getName().equals("*1")) {
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
    variantReport.setWildtypeAllele(m_wildAlleleMap.get(variantReport.getPosition()));
    return variantReport;
  }

  /**
   * Fills in the named alleles that the given {@link VariantReport} is used to define
   */
  private void assignNamedAlleles(VariantReport variantReport) {
    Collection<String> alleles = m_variantAlleleMap.get(variantReport.getPosition());
    if (alleles != null) {
      variantReport.setAlleles(alleles);
    }
  }
}
