package org.pharmgkb.pharmcat.reporter;

import java.io.IOException;
import org.pharmgkb.pharmcat.definition.VariantAlleleMap;
import org.pharmgkb.pharmcat.definition.model.VariantLocus;
import org.pharmgkb.pharmcat.haplotype.model.Variant;
import org.pharmgkb.pharmcat.reporter.model.VariantReport;


/**
 * Factory class to help create new {@link VariantReport} instances for a given gene
 *
 * @author Ryan Whaley
 */
public class VariantReportFactory {

  private String m_gene;
  private VariantAlleleMap m_variantAlleleMap;

  public VariantReportFactory(String gene) throws IOException {
    m_gene = gene;
    m_variantAlleleMap = new VariantAlleleMap(gene);
  }

  public VariantReport make(Variant variant) {
    return new VariantReport(m_gene, variant).findAlleles(m_variantAlleleMap);
  }

  public VariantReport make(VariantLocus locus) {
    return new VariantReport(m_gene, locus).findAlleles(m_variantAlleleMap);
  }
}
