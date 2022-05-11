
package org.pharmgkb.pharmcat.reporter.model.pgkb;

import com.google.gson.annotations.Expose;
import com.google.gson.annotations.SerializedName;


/**
 * PharmGKB Ontology Term Model
 */
public class OntologyTerm {

    @SerializedName("term")
    @Expose
    private String term;
    @SerializedName("termId")
    @Expose
    private String termId;


    public String getTerm() {
        return term;
    }

    public void setTerm(String term) {
        this.term = term;
    }

    public String getTermId() {
        return termId;
    }

    public void setTermId(String termId) {
        this.termId = termId;
    }
}
