
package org.pharmgkb.pharmcat.reporter.model;

import javax.annotation.Generated;
import com.google.gson.annotations.Expose;
import com.google.gson.annotations.SerializedName;

@Generated("org.jsonschema2pojo")
public class RelatedGene {

    @SerializedName("objCls")
    @Expose
    private String objCls;
    @SerializedName("@id")
    @Expose
    private String Id;
    @SerializedName("@context")
    @Expose
    private String Context;
    @SerializedName("id")
    @Expose
    private String id;
    @SerializedName("symbol")
    @Expose
    private String symbol;
    @SerializedName("name")
    @Expose
    private String name;

    /**
     * 
     * @return
     *     The objCls
     *
    public String getObjCls() {
        return objCls;
    }

    /**
     * 
     * @param objCls
     *     The objCls
     *
    public void setObjCls(String objCls) {
        this.objCls = objCls;
    }

    /**
     * 
     * @return
     *     The Id
     *
    public String getId() {
        return Id;
    }

    /**
     * 
     * @param Id
     *     The @id
     *
    public void setId(String Id) {
        this.Id = Id;
    }

    /**
     * 
     * @return
     *     The Context
     */
    public String getContext() {
        return Context;
    }

    /**
     * 
     * @param Context
     *     The @context
     */
    public void setContext(String Context) {
        this.Context = Context;
    }

    /**
     * 
     * @return
     *     The id
     */
    public String getId() {
        return id;
    }

    /**
     * 
     * @param id
     *     The id
     */
    public void setId(String id) {
        this.id = id;
    }

    /**
     * 
     * @return
     *     The symbol
     */
    public String getSymbol() {
        return symbol;
    }

    /**
     * 
     * @param symbol
     *     The symbol
     */
    public void setSymbol(String symbol) {
        this.symbol = symbol;
    }

    /**
     * 
     * @return
     *     The name
     */
    public String getName() {
        return name;
    }

    /**
     * 
     * @param name
     *     The name
     */
    public void setName(String name) {
        this.name = name;
    }

}
