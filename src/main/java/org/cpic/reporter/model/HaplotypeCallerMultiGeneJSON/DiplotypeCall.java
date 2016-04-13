
package org.cpic.reporter.model.HaplotypeCallerMultiGeneJSON;

import java.util.ArrayList;
import java.util.List;
import javax.annotation.Generated;
import com.google.gson.annotations.Expose;
import com.google.gson.annotations.SerializedName;


public class DiplotypeCall {

    @SerializedName("gene")
    @Expose
    private String gene;
    @SerializedName("geneVersion")
    @Expose
    private String geneVersion;
    @SerializedName("chromosome")
    @Expose
    private String chromosome;
    @SerializedName("diplotypes")
    @Expose
    private List<String> diplotypes = new ArrayList<String>();
    @SerializedName("positionsUsed")
    @Expose
    private List<PositionsUsed> positionsUsed = new ArrayList<PositionsUsed>();
    @SerializedName("haplotypesNotCalled")
    @Expose
    private List<HaplotypesNotCalled> haplotypesNotCalled = new ArrayList<HaplotypesNotCalled>();
    @SerializedName("metaData")
    @Expose
    private MetaData metaData;

    /**
     * 
     * @return
     *     The gene
     */
    public String getGene() {
        return gene;
    }

    /**
     * 
     * @param gene
     *     The gene
     */
    public void setGene(String gene) {
        this.gene = gene;
    }

    /**
     * 
     * @return
     *     The geneVersion
     */
    public String getGeneVersion() {
        return geneVersion;
    }

    /**
     * 
     * @param geneVersion
     *     The geneVersion
     */
    public void setGeneVersion(String geneVersion) {
        this.geneVersion = geneVersion;
    }

    /**
     * 
     * @return
     *     The chromosome
     */
    public String getChromosome() {
        return chromosome;
    }

    /**
     * 
     * @param chromosome
     *     The chromosome
     */
    public void setChromosome(String chromosome) {
        this.chromosome = chromosome;
    }

    /**
     * 
     * @return
     *     The diplotypes
     */
    public List<String> getDiplotypes() {
        return diplotypes;
    }

    /**
     * 
     * @param diplotypes
     *     The diplotypes
     */
    public void setDiplotypes(List<String> diplotypes) {
        this.diplotypes = diplotypes;
    }

    /**
     * 
     * @return
     *     The positionsUsed
     */
    public List<PositionsUsed> getPositionsUsed() {
        return positionsUsed;
    }

    /**
     * 
     * @param positionsUsed
     *     The positionsUsed
     */
    public void setPositionsUsed(List<PositionsUsed> positionsUsed) {
        this.positionsUsed = positionsUsed;
    }

    /**
     * 
     * @return
     *     The haplotypesNotCalled
     */
    public List<HaplotypesNotCalled> getHaplotypesNotCalled() {
        return haplotypesNotCalled;
    }

    /**
     * 
     * @param haplotypesNotCalled
     *     The haplotypesNotCalled
     */
    public void setHaplotypesNotCalled(List<HaplotypesNotCalled> haplotypesNotCalled) {
        this.haplotypesNotCalled = haplotypesNotCalled;
    }

    /**
     * 
     * @return
     *     The metaData
     */
    public MetaData getMetaData() {
        return metaData;
    }

    /**
     * 
     * @param metaData
     *     The metaData
     */
    public void setMetaData(MetaData metaData) {
        this.metaData = metaData;
    }

}
