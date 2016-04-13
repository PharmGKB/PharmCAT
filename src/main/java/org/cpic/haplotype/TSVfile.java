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
	public void setGeneID(String _GeneID){
		GeneID=_GeneID;
	}
	
	public String getGeneID(){
		return GeneID;
	}
	
	public void setGeneName(String _GeneName){
		GeneName=_GeneName;
	}
	
	public String getGeneName(){
		return GeneName;
	}
	public void setContentVersion(String _ContentVersion){
		ContentVersion=_ContentVersion;
	}
	
	public String getContentVersion(){
		return ContentVersion;
	}
	public void setContentDate(String _ContentDate){
		ContentDate=_ContentDate;
	}
	
	public String getContentDate(){
		return ContentDate;
	}
	public void setGenomeBuild(String _GenomeBuild){
		GenomeBuild=_GenomeBuild;
	}
	
	public String getGenomeBuild(){
		return GenomeBuild;
	}
	public void setProteinID(String _ProteinID){
		ProteinID=_ProteinID;
	}
	
	public String getProteinID(){
		return ProteinID;
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
