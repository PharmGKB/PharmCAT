
package org.pharmgkb.pharmcat.reporter.model;

import javax.annotation.Generated;
import com.google.gson.annotations.Expose;
import com.google.gson.annotations.SerializedName;

@Generated("org.jsonschema2pojo")
public class OntologyTerm {

    @SerializedName("src")
    @Expose
    private String src;
    @SerializedName("term")
    @Expose
    private String term;
    @SerializedName("termId")
    @Expose
    private String termId;


    /**
     * 
     * @return
     *     The src
     */
    public String getSrc() {
        return src;
    }

    /**
     * 
     * @param src
     *     The src
     */
    public void setSrc(String src) {
        this.src = src;
    }

    /**
     * 
     * @return
     *     The term
     */
    public String getTerm() {
        return term;
    }

    /**
     * 
     * @param term
     *     The term
     */
    public void setTerm(String term) {
        this.term = term;
    }

    /**
     * 
     * @return
     *     The termId
     */
    public String getTermId() {
        return termId;
    }

    /**
     * 
     * @param termId
     *     The termId
     */
    public void setTermId(String termId) {
        this.termId = termId;
    }

}
