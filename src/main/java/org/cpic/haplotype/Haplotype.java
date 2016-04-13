package org.cpic.haplotype;

import java.util.ArrayList;

import org.apache.commons.lang3.ObjectUtils;
import org.cpic.util.HaplotypeNameComparator;

public class Haplotype implements Comparable<Haplotype> {

	private ArrayList <Variant> Variants = new ArrayList<>();
	private String AlleleID;
	private String CommonName;
	private String FunctionStatus;
	private ArrayList <String> Alleles = new ArrayList<>();
	
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
	public ArrayList<Variant> getVaraints(){
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
	public void addAlleles(ArrayList <String> _Alleles){
		Alleles = _Alleles;
	}
	public ArrayList<String> getAlleles(){
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
}
