
package org.cpic.reporter.model.HaplotypeCallerMultiGeneJSON;

import java.util.ArrayList;
import java.util.List;
import javax.annotation.Generated;
import com.google.gson.annotations.Expose;
import com.google.gson.annotations.SerializedName;


public class HaplotypeCallsMultiGene {

    @SerializedName("diplotypeCalls")
    @Expose
    private List<DiplotypeCall> diplotypeCalls = new ArrayList<DiplotypeCall>();

    /**
     * 
     * @return
     *     The diplotypeCalls
     */
    public List<DiplotypeCall> getDiplotypeCalls() {
        return diplotypeCalls;
    }

    /**
     * 
     * @param diplotypeCalls
     *     The diplotypeCalls
     */
    public void setDiplotypeCalls(List<DiplotypeCall> diplotypeCalls) {
        this.diplotypeCalls = diplotypeCalls;
    }

}
