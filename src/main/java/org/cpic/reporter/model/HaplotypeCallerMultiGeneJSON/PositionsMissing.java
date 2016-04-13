
package org.cpic.reporter.model.HaplotypeCallerMultiGeneJSON;

import javax.annotation.Generated;
import com.google.gson.annotations.Expose;
import com.google.gson.annotations.SerializedName;


public class PositionsMissing {

    @SerializedName("position")
    @Expose
    private String position;
    @SerializedName("rsID")
    @Expose
    private String rsID;

    /**
     * 
     * @return
     *     The position
     */
    public String getPosition() {
        return position;
    }

    /**
     * 
     * @param position
     *     The position
     */
    public void setPosition(String position) {
        this.position = position;
    }

    /**
     * 
     * @return
     *     The rsID
     */
    public String getRsID() {
        return rsID;
    }

    /**
     * 
     * @param rsID
     *     The rsID
     */
    public void setRsID(String rsID) {
        this.rsID = rsID;
    }

}
