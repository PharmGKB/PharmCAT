
package org.cpic.reporter.model.HaplotypeCallerJSON;

import java.util.ArrayList;
import java.util.List;
import javax.annotation.Generated;
import com.google.gson.annotations.Expose;
import com.google.gson.annotations.SerializedName;

@Generated("org.jsonschema2pojo")
public class HaplotypesNotCalled {

    @SerializedName("haplotype")
    @Expose
    private String haplotype;
    @SerializedName("positionsMissing")
    @Expose
    private List<PositionsMissing> positionsMissing = new ArrayList<PositionsMissing>();

    /**
     * 
     * @return
     *     The haplotype
     */
    public String getHaplotype() {
        return haplotype;
    }

    /**
     * 
     * @param haplotype
     *     The haplotype
     */
    public void setHaplotype(String haplotype) {
        this.haplotype = haplotype;
    }

    /**
     * 
     * @return
     *     The positionsMissing
     */
    public List<PositionsMissing> getPositionsMissing() {
        return positionsMissing;
    }

    /**
     * 
     * @param positionsMissing
     *     The positionsMissing
     */
    public void setPositionsMissing(List<PositionsMissing> positionsMissing) {
        this.positionsMissing = positionsMissing;
    }

}
