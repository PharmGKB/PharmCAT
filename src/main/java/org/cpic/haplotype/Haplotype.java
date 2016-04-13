package org.cpic.haplotype;

import java.util.ArrayList;

import org.apache.commons.lang3.ObjectUtils;
import org.cpic.util.HaplotypeNameComparator;

public class Haplotype implements Comparable<Haplotype> {

	private ArrayList <Variant> Variants = new ArrayList<>();
	private String AlleleID;
	private String CommonName;
	private String FunctionStatus;
	private ArrayList <String> NormalAlleles = new ArrayList<>();
	private ArrayList <String>  EffectAlleles = new ArrayList<>();
	
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
	public void addNormalAlleles(String _NormalAlleles){
		NormalAlleles.add(_NormalAlleles);
	}
	public ArrayList<String> getNormalAlleles(){
		return NormalAlleles;
	}
	public void addEffectAlleles(String _EffectAlleles){
		EffectAlleles.add(_EffectAlleles);
	}
	public ArrayList<String> getEffectAlleles(){
		return EffectAlleles;
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
