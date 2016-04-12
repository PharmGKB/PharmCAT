package org.cpic.haplotype;

import java.util.ArrayList;

import org.apache.commons.lang3.ObjectUtils;
import org.cpic.util.HaplotypeNameComparator;

public class Haplotype implements Comparable<Haplotype> {

	ArrayList <Variant> varaints;
	String AlleleID;
	String CommonName;
	String FunctionStatus;
	ArrayList <String> NormalAlleles;
	ArrayList <String>  EffectAlleles;


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
