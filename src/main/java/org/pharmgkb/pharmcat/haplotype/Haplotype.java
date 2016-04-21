package org.pharmgkb.pharmcat.haplotype;

import java.util.ArrayList;
import java.util.List;
import java.util.regex.Pattern;
import javax.annotation.Nonnull;
import org.apache.commons.lang3.ObjectUtils;
import org.pharmgkb.pharmcat.util.HaplotypeNameComparator;


public class Haplotype implements Comparable<Haplotype> {

  private List<Variant> m_variants = new ArrayList<>();
  private String m_alleleId;
  private String m_name;
  private String m_functionStatus;
  private List<String> m_alleles = new ArrayList<>();
  private Pattern m_permutations;


  public Haplotype(String alleleId, String name) {
    m_alleleId = alleleId;
    m_name = name;
  }

  public Haplotype(String alleleId, String name, List<Variant> variants, String functionStatus, List<String> alleles) {
    m_variants = variants;
    m_alleleId = alleleId;
    m_name = name;
    m_functionStatus = functionStatus;
    m_alleles = alleles;
  }

  public void addVariant(Variant _Variant) {
    m_variants.add(_Variant);
  }

  public List<Variant> getVariants() {
    return m_variants;
  }

  public void setAlleleId(String _AlleleID) {
    m_alleleId = _AlleleID;
  }

  public String getAlleleId() {
    return m_alleleId;
  }

  public void setName(String _CommonName) {
    m_name = _CommonName;
  }

  public String getName() {
    return m_name;
  }

  public void setFunctionStatus(String _FunctionStatus) {
    m_functionStatus = _FunctionStatus;
  }

  public String getFunctionStatus() {
    return m_functionStatus;
  }

  public void addAllele(String _Allele) {
    m_alleles.add(_Allele);
  }

  public void addAlleles(List<String> _Alleles) {
    m_alleles = _Alleles;
  }

  public List<String> getAlleles() {
    return m_alleles;
  }


  @Override
  public String toString() {
    return m_name;
  }


  @Override
  public int compareTo(@Nonnull Haplotype o) {
    int rez = HaplotypeNameComparator.getComparator().compare(m_name, o.m_name);
    if (rez != 0) {
      return rez;
    }
    return ObjectUtils.compare(m_alleleId, o.m_alleleId);
  }


  private String resolveIupacCode(@Nonnull String allele) {
    if (allele.length() == 1) {
      return Iupac.lookup(allele).getRegex();
    }
    return allele;
  }

  public Pattern getPermutations() {
    return m_permutations;
  }

  public Pattern calculatePermutations(List<Variant> allVariants) {


    StringBuilder builder = new StringBuilder();
    int idx = 0;
    for (int x = 0; x < allVariants.size(); x++) {
      builder.append(allVariants.get(x).getPosition())
          .append(":");
      if (idx < m_variants.size() && m_variants.get(idx) == allVariants.get(x)) {
        builder.append(resolveIupacCode(m_alleles.get(idx)));
        idx += 1;
      } else {
        builder.append(".?");
      }
      builder.append(";");
    }
    m_permutations = Pattern.compile(builder.toString());
    return m_permutations;
  }
}
