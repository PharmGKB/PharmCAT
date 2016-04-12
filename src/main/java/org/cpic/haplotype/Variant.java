package org.cpic.haplotype;

import org.apache.commons.lang3.ObjectUtils;

public class Variant implements Comparable<Variant> {
	String CHROM;
	int POS;
	String cDNA;
	String ProteingEffect;
	String [] ALTs;
	String REF;
	String rsID;
	String GENE;


  @Override
  public int compareTo(Variant o) {

    int rez = ObjectUtils.compare(CHROM, o.CHROM);
    if (rez != 0) {
      return rez;
    }
    rez = ObjectUtils.compare(POS, o.POS);
    if (rez != 0) {
      return rez;
    }
    return ObjectUtils.compare(GENE, o.GENE);
  }
}
