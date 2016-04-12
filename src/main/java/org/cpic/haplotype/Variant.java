package org.cpic.haplotype;

import java.util.ArrayList;

import org.apache.commons.lang3.ObjectUtils;

public class Variant implements Comparable<Variant> {
	
	private String CHROM;
	private String GeneName;
	private int POS;
	private int GenePOS;
	private ArrayList<String> cDNA;
	private ArrayList<String> ProteingEffect;
	private ArrayList<String> ALTs;
	private String REF;
	private String rsID;
	
	
	public Variant(){
	}
	public Variant(String _CHROM, String _GeneName,String _cDNA){
		CHROM = _CHROM;
		GeneName = _GeneName;
		if (_cDNA.contains(";")){
			String [] fields = _cDNA.split(";");
			for (int i = 0; i < fields.length; i++){
				cDNA.add(fields[i].trim());
			}
		}
		cDNA.add(_cDNA);
	}
	public String getCHROM(){
		return CHROM;
	}
	public String getGeneName(){
		return GeneName;
	}
	public void set_cDNA(String _cDNA){
		if (_cDNA.contains(";")){
			String [] fields = _cDNA.split(";");
			for (int i = 0; i < fields.length; i++){
				cDNA.add(fields[i].trim());
			}
		}
		cDNA.add(_cDNA);
	}
	public void setProteingEffect(String _ProteingEffect){
		if (_ProteingEffect.contains(";")){
			String [] fields = _ProteingEffect.split(";");
			for (int i = 0; i < fields.length; i++){
				ProteingEffect.add(fields[i].trim());
			}
		}
			
		ProteingEffect.add(_ProteingEffect);
	}
	public void setPOS(String _POS){
		POS = Integer.parseInt(_POS);
	}
	public int getPOS(){
		return POS;
	}
	public void setGenePOS(String _GenePOS){
		GenePOS = Integer.parseInt(_GenePOS);
	}
	public void setREF(String _REF){
		REF = (_REF);
	}
	public void set_rsID(String _rsID){
		REF = (_rsID);
	}
	public void addALT(String _ALT){
		ALTs.add(_ALT);
	}
	public ArrayList<String> getALTs(){
		return ALTs;
	}





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
    return ObjectUtils.compare(GeneName, o.GeneName);
  }
}
