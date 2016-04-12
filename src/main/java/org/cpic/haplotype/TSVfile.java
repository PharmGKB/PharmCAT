package org.cpic.haplotype;

public class TSVfile {
	private String FileName;
	private int FormatVersion;
	private String GeneName;
	private String GeneID;
	private String ContentVersion;
	private String ContentDate;
	private String GenomeBuild;
	private String ChromosomeID;
	private String Chromosome;
	private String ProteinID;
	
	public TSVfile(String _FileName){
		FileName=_FileName;
	}
	
	public String getFileName(){
		return FileName;
	}
	

	public void setFormatVersion(String _FormatVersion){
		FormatVersion=Integer.parseInt(_FormatVersion);
	}
	
	public int getFormatVersion(){
		return FormatVersion;
	}
	
	public void setGeneName(String _GeneName){
		GeneName=_GeneName;
	}
	
	public String getGeneName(){
		return GeneName;
	}
	
	public void setChromosome(String _Chromosome){
		Chromosome=_Chromosome;
	}
	
	public String getChromosome(){
		return Chromosome;
	}
	
	public void setChromosomeID(String _ChromosomeID){
		Chromosome=_ChromosomeID;
	}
	
	public String getChromosomeID(){
		return ChromosomeID;
	}
	


}
