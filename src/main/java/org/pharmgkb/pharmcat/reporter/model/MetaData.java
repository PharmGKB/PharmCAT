package org.pharmgkb.pharmcat.reporter.model;

public class MetaData {
    String cpicAnnotatorBuild;
    String cpicDataBuild;
    String genomeAssembly;
    String inputFile;
    String dateTime;
    
    public String getCpicAnnotatorBuild() {
        return cpicAnnotatorBuild;
    }
    public void setCpicAnnotatorBuild(String cpicAnnotatorBuild) {
        this.cpicAnnotatorBuild = cpicAnnotatorBuild;
    }
    public String getCpicDataBuild() {
        return cpicDataBuild;
    }
    public void setCpicDataBuild(String cpicDataBuild) {
        this.cpicDataBuild = cpicDataBuild;
    }
    public String getGenomeAssembly() {
        return genomeAssembly;
    }
    public void setGenomeAssembly(String genomeAssembly) {
        this.genomeAssembly = genomeAssembly;
    }
    public String getInputFile() {
        return inputFile;
    }
    public void setInputFile(String inputFile) {
        this.inputFile = inputFile;
    }
    public String getDateTime() {
        return dateTime;
    }
    public void setDateTime(String dateTime) {
        this.dateTime = dateTime;
    }

}
