package org.pharmgkb.pharmcat.reporter.resultsJSON;

public class ReporterResult {
	
	String diplotype;
	String gene;
	String implications;
	String phenotype;
	String recommendations;
	String recommendationsClass;
	
	public String getDiplotype() {
		return diplotype;
	}
	public void setDiplotype(String diplotype) {
		this.diplotype = diplotype;
	}
	public String getGene() {
		return gene;
	}
	public void setGene(String gene) {
		this.gene = "test";
	}
	public String getImplications() {
		return implications;
	}
	public void setImplications(String implications) {
		this.implications = implications;
	}
	public String getPhenotype() {
		return phenotype;
	}
	public void setPhenotype(String phenotype) {
		this.phenotype = phenotype;
	}
	public String getRecommendations() {
		return recommendations;
	}
	public void setRecommendations(String recommendations) {
		this.recommendations = recommendations;
	}
	public String getRecommendationsClass() {
		return recommendationsClass;
	}
	public void setRecommendationsClass(String recommendationsClass) {
		this.recommendationsClass = recommendationsClass;
	}
	
}
