package org.pharmgkb.pharmcat.reporter.model;

import java.util.ArrayList;

public class CalledType {
    String gene;
    String geneVersion;
    String chromosome;
    ArrayList<String> diplotypes;
    MetaData metaData;
    
    public String getGene() {
        return gene;
    }
    public void setGene(String gene) {
        this.gene = gene;
    }
    public String getGeneVersion() {
        return geneVersion;
    }
    public void setGeneVersion(String geneVersion) {
        this.geneVersion = geneVersion;
    }
    public String getChromosome() {
        return chromosome;
    }
    public void setChromosome(String chromosome) {
        this.chromosome = chromosome;
    }
    public ArrayList<String> getDiplotypes() {
        return diplotypes;
    }
    public void setDiplotypes(ArrayList<String> diplotypes) {
        this.diplotypes = diplotypes;
    }
    public MetaData getMetaData() {
        return metaData;
    }
    public void setMetaData(MetaData metaData) {
        this.metaData = metaData;
    }
    

}
