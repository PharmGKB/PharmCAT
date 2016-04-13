
package org.cpic.reporter.model.HaplotypeCallerMultiGeneJSON;

import javax.annotation.Generated;
import com.google.gson.annotations.Expose;
import com.google.gson.annotations.SerializedName;


public class MetaData {

    @SerializedName("cpicAnnotatorBuild")
    @Expose
    private String cpicAnnotatorBuild;
    @SerializedName("cpicDataBuild")
    @Expose
    private String cpicDataBuild;
    @SerializedName("genomeAssembly")
    @Expose
    private String genomeAssembly;
    @SerializedName("inputFile")
    @Expose
    private String inputFile;
    @SerializedName("dateTime")
    @Expose
    private String dateTime;

    /**
     * 
     * @return
     *     The cpicAnnotatorBuild
     */
    public String getCpicAnnotatorBuild() {
        return cpicAnnotatorBuild;
    }

    /**
     * 
     * @param cpicAnnotatorBuild
     *     The cpicAnnotatorBuild
     */
    public void setCpicAnnotatorBuild(String cpicAnnotatorBuild) {
        this.cpicAnnotatorBuild = cpicAnnotatorBuild;
    }

    /**
     * 
     * @return
     *     The cpicDataBuild
     */
    public String getCpicDataBuild() {
        return cpicDataBuild;
    }

    /**
     * 
     * @param cpicDataBuild
     *     The cpicDataBuild
     */
    public void setCpicDataBuild(String cpicDataBuild) {
        this.cpicDataBuild = cpicDataBuild;
    }

    /**
     * 
     * @return
     *     The genomeAssembly
     */
    public String getGenomeAssembly() {
        return genomeAssembly;
    }

    /**
     * 
     * @param genomeAssembly
     *     The genomeAssembly
     */
    public void setGenomeAssembly(String genomeAssembly) {
        this.genomeAssembly = genomeAssembly;
    }

    /**
     * 
     * @return
     *     The inputFile
     */
    public String getInputFile() {
        return inputFile;
    }

    /**
     * 
     * @param inputFile
     *     The inputFile
     */
    public void setInputFile(String inputFile) {
        this.inputFile = inputFile;
    }

    /**
     * 
     * @return
     *     The dateTime
     */
    public String getDateTime() {
        return dateTime;
    }

    /**
     * 
     * @param dateTime
     *     The dateTime
     */
    public void setDateTime(String dateTime) {
        this.dateTime = dateTime;
    }

}
