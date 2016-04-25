package org.pharmgkb.pharmcat.reporter.resultsJSON;

import java.util.List;
import org.pharmgkb.pharmcat.reporter.model.CPICinteraction;


public class MultiGeneInteraction extends Interaction{
  

	private List<Gene> relevantGene;
	private List<String> relevantDipCalls;
	
	public MultiGeneInteraction( ){
        super();
	}
	
	public MultiGeneInteraction( CPICinteraction inter){
	    	setName( inter.getName());
	        setRelatedChemicals( inter.getRelatedChemicals());
	        setSource( inter.getSource());
	        setSummaryHtml( inter.getSummaryHtml());
	        setTextHtml( inter.getTextHtml());        
	 }
	  
	public List<Gene> getRelevantGene() {
		return relevantGene;
	}
	public void setRelevantGene(List<Gene> relevantGene) {
		this.relevantGene = relevantGene;
	}
	public List<String> getRelevantDipCalls() {
		return relevantDipCalls;
	}
	public void setRelevantDipCalls(List<String> relevantDipCalls) {
		this.relevantDipCalls = relevantDipCalls;
	}
	

    
}
