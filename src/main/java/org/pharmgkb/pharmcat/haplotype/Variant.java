package org.pharmgkb.pharmcat.haplotype;

import java.util.ArrayList;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import javax.annotation.Nonnull;
import org.apache.commons.lang3.ObjectUtils;


/**
 * Object to hold information about varaints from the tsv file.
 *
 * @author nate
 *
 */
public class Variant implements Comparable<Variant> {

	private String CHROM;
	private String GeneName;
	private int POS;
	private String GenePOS;
	private ArrayList<String> HGVSg  = new ArrayList<>();
	private ArrayList<String> cDNA = new ArrayList<>();
	private ArrayList<String> ProteinEffect  = new ArrayList<>();
	private ArrayList<String> ALTs  = new ArrayList<>();
	private String REF;
	private String rsID;

	/**
	 * @param _CHROM
	 * @param _GeneName
	 * @param _cDNA
	 */
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
	/**
	 * @param _CHROM
	 * @param _GeneName
	 */
	public Variant(String _CHROM, String _GeneName){
		CHROM = _CHROM;
		GeneName = _GeneName;

	}
	/**
	 * @return
	 */
	public String getCHROM(){
		return CHROM;
	}
	/**
	 * @return
	 */
	public String getGeneName(){
		return GeneName;
	}
	/**
	 * @param _cDNA
	 */
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
	/**
	 * @param _ProteingEffect
	 */
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
	/**
	 * @param _HGVSg
	 */
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
	public String getHGVSg(){
		return HGVSg.toString();
	}
	/**
	 * @param _HGVSg
	 * @return
	 */
	private int getStartPOS(String _HGVSg){

		  Pattern p = Pattern.compile("(\\d+)");
		  Matcher m = p.matcher(_HGVSg);
		  //System.out.println(_HGVSg);
		  if(m.find()){
		  System.out.println(m.group(1));

		  return Integer.parseInt(m.group(1));
		  }
		  else {return -1;}

	}
	/**
	 * @return
	 */
	public boolean setStartPOS(){
		if (HGVSg.isEmpty()){
			return false;
		}else{
			POS = getStartPOS(HGVSg.get(0));
			return true;
		}
	}
	/**
	 * @param _POS
	 */
	public void setPOS(String _POS){
		POS = Integer.parseInt(_POS);
	}
	/**
	 * @return
	 */
	public int getPOS(){
		return POS;
	}
	/**
	 * @param _GenePOS
	 */
	public void setGenePOS(String _GenePOS){
		GenePOS = _GenePOS;
	}
	/**
	 * @return
	 */
	public String getGenePOS(){
		return GenePOS;
	}
	/**
	 * @param _REF
	 */
	public void setREF(String _REF){
		REF = (_REF);
	}
	/**
	 * @return
	 */
	public String getREF(){
		return REF;
	}
	/**
	 * @param _rsID
	 */
	public void set_rsID(String _rsID){
		rsID = (_rsID);
	}
	/**
	 * @return
	 */
	public String get_rsID(){
		return rsID;
	}
	/**
	 * @param _ALT
	 */
	public void addALT(String _ALT){
		ALTs.add(_ALT);
	}
	/**
	 * @return
	 */
	public ArrayList<String> getALTs(){
		return ALTs;
	}






  @Override
  public int compareTo(@Nonnull Variant o) {

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
