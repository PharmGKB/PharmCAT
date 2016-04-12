package org.cpic.haplotype;

import org.apache.commons.lang3.ObjectUtils;

import java.util.ArrayList;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

public class Variant implements Comparable<Variant> {

	private String CHROM;
	private String GeneName;
	private int POS;
	private int GenePOS;
	private ArrayList<String> HGVSg;
	private ArrayList<String> cDNA = new ArrayList<>();
	private ArrayList<String> ProteinEffect;
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
		}else{
			cDNA.add(_cDNA);
		}
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
		}else{
			cDNA.add(_cDNA);
		}
	}
	public void addProteingEffect(String _ProteingEffect){
		if (_ProteingEffect.contains(";")){
			String [] fields = _ProteingEffect.split(";");
			for (int i = 0; i < fields.length; i++){
				ProteinEffect.add(fields[i].trim());
			}
		}else{

			ProteinEffect.add(_ProteingEffect);
		}
	}
	public void addHGVSg(String _HGVSg){
		if (_HGVSg.contains(";")){
			String [] fields = _HGVSg.split(";");
			for (int i = 0; i < fields.length; i++){
				HGVSg.add(fields[i].trim());
			}
		}else{

			HGVSg.add(_HGVSg);
		}
	}
	private int getStartPOS(String _HGVSg){

		  Pattern p = Pattern.compile("\\d+");
		  Matcher m = p.matcher(_HGVSg);

		  return Integer.parseInt(m.group(1));

	}
	public boolean setStartPOS(){
		if (HGVSg.isEmpty()){
			return false;
		}else{
			POS = getStartPOS(HGVSg.get(0));
			return true;
		}
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
		rsID = (_rsID);
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
