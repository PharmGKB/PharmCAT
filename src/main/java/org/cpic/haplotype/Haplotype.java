package org.cpic.haplotype;

import org.apache.commons.lang3.ObjectUtils;
import org.cpic.util.HaplotypeNameComparator;

import javax.annotation.Nonnull;
import java.util.ArrayList;
import java.util.List;
import java.util.regex.Pattern;

public class Haplotype implements Comparable<Haplotype> {

	private List <Variant> Variants = new ArrayList<>();
	private String AlleleID;
	private String CommonName;
	private String FunctionStatus;
	private List <String> Alleles = new ArrayList<>();
  private Pattern m_permutations;


	public Haplotype(){
	}

	public Haplotype(ArrayList <Variant> _Variants, String _AlleleID, String _CommonName, String _FunctionStatus, ArrayList <String> _Alleles){
		Variants=_Variants;
		AlleleID = _AlleleID;
		CommonName = _CommonName;
		FunctionStatus = _FunctionStatus;
		Alleles=_Alleles;
	}

	public void addVariant(Variant _Variant){
		Variants.add(_Variant);
	}
	public List<Variant> getVariants(){
		return Variants;
	}
	public void setAlleleID(String _AlleleID){
		AlleleID = _AlleleID;
	}
	public String getAlleleID(){
		return AlleleID;
	}
	public void setCommonName(String _CommonName){
		CommonName = _CommonName;
	}
	public String getCommonName(){
		return CommonName;
	}
	public void setFunctionStatus(String _FunctionStatus){
		FunctionStatus = _FunctionStatus;
	}
	public String getFunctionStatus(){
		return FunctionStatus;
	}
	public void addAllele(String _Allele){
		Alleles.add(_Allele);
	}
	public void addAlleles(List <String> _Alleles){
		Alleles = _Alleles;
	}
	public List<String> getAlleles(){
		return Alleles;
	}


  @Override
  public String toString() {
    return CommonName;
  }


  @Override
  public int compareTo(Haplotype o) {
    int rez = HaplotypeNameComparator.getComparator().compare(CommonName, o.CommonName);
    if (rez != 0) {
      return rez;
    }
    return ObjectUtils.compare(AlleleID, o.AlleleID);
  }


  private String resolveIupacCode(@Nonnull String allele) {
    if (allele.length() == 1) {
      return Iupac.lookup(allele).getRegex();
    }
    return allele;
  }

  public Pattern calculatePermutations(List<Variant> allVariants) {


    StringBuilder builder = new StringBuilder();
    int idx = 0;
    for (int x = 0; x < allVariants.size(); x++) {
      builder.append(allVariants.get(x).getPOS())
          .append(":");
      if (idx < Variants.size() && Variants.get(idx) == allVariants.get(x)) {
        builder.append(resolveIupacCode(Alleles.get(idx)));
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
