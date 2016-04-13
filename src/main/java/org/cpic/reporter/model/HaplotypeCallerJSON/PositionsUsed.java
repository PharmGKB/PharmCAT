
package org.cpic.reporter.model.HaplotypeCallerJSON;

import javax.annotation.Generated;
import com.google.gson.annotations.Expose;
import com.google.gson.annotations.SerializedName;

@Generated("org.jsonschema2pojo")
public class PositionsUsed {

    @SerializedName("position")
    @Expose
    private String position;
    @SerializedName("rsID")
    @Expose
    private String rsID;
    @SerializedName("VCFCall")
    @Expose
    private String VCFCall;

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

    /**
     * 
     * @return
     *     The VCFCall
     */
    public String getVCFCall() {
        return VCFCall;
    }

    /**
     * 
     * @param VCFCall
     *     The VCFCall
     */
    public void setVCFCall(String VCFCall) {
        this.VCFCall = VCFCall;
    }

}
