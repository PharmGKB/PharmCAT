package org.cpic.haplotype;

/**
 * Object to hold haplotype metadata
 * @author nate
 *
 */

public class TSVfile {
	private String FileName;
	private int FormatVersion;
	private String GeneName;
	private String GeneID;
	private String GeneOrientation;
	private String ContentVersion;
	private String ContentDate;
	private String GenomeBuild;
	private String ChromosomeID;
	private String Chromosome;
	private String ProteinID;
	
	/**
	 * @param _FileName
	 */
	public TSVfile(String _FileName){
		FileName=_FileName;
	}
	
	public String getFileName(){
		return FileName;
	}
	

	/**
	 * @param _FormatVersion
	 */
	public void setFormatVersion(String _FormatVersion){
		FormatVersion=Integer.parseInt(_FormatVersion);
	}
	
	public int getFormatVersion(){
		return FormatVersion;
	}
	/**
	 * @param _GeneID
	 */
	public void setGeneID(String _GeneID){
		GeneID=_GeneID;
	}
	
	public String getGeneID(){
		return GeneID;
	}
	
	/**
	 * @param _GeneOrientation
	 */
	public void setGeneOrientation(String _GeneOrientation){
		GeneOrientation=_GeneOrientation;
	}
	
	/**
	 * @return
	 */
	public String getGeneOrientation(){
		return GeneOrientation;
	}
	/**
	 * @param _GeneName
	 */
	public void setGeneName(String _GeneName){
		GeneName=_GeneName;
	}
	
	/**
	 * @return
	 */
	public String getGeneName(){
		return GeneName;
	}
	/**
	 * @param _ContentVersion
	 */
	public void setContentVersion(String _ContentVersion){
		ContentVersion=_ContentVersion;
	}
	
	/**
	 * @return
	 */
	public String getContentVersion(){
		return ContentVersion;
	}
	
	/**
	 * @param _ContentDate
	 */
	public void setContentDate(String _ContentDate){
		ContentDate=_ContentDate;
	}
	
	/**
	 * @return
	 */
	public String getContentDate(){
		return ContentDate;
	}
	/**
	 * @param _GenomeBuild
	 */
	public void setGenomeBuild(String _GenomeBuild){
		GenomeBuild=_GenomeBuild;
	}
	
	/**
	 * @return
	 */
	public String getGenomeBuild(){
		return GenomeBuild;
	}
	/**
	 * @param _ProteinID
	 */
	public void setProteinID(String _ProteinID){
		ProteinID=_ProteinID;
	}
	
	/**
	 * @return
	 */
	public String getProteinID(){
		return ProteinID;
	}
	
	/**
	 * @param _Chromosome
	 */
	public void setChromosome(String _Chromosome){
		Chromosome=_Chromosome;
	}
	
	/**
	 * @return
	 */
	public String getChromosome(){
		return Chromosome;
	}
	
	/**
	 * @param _ChromosomeID
	 */
	public void setChromosomeID(String _ChromosomeID){
		Chromosome=_ChromosomeID;
	}
	
	/**
	 * @return
	 */
	public String getChromosomeID(){
		return ChromosomeID;
	}
	


}
