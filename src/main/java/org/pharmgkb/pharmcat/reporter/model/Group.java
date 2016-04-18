
package org.pharmgkb.pharmcat.reporter.model;

import java.util.ArrayList;
import java.util.List;
import javax.annotation.Generated;
import com.google.gson.annotations.Expose;
import com.google.gson.annotations.SerializedName;

@Generated("org.jsonschema2pojo")
public class Group {

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
    @SerializedName("name")
    @Expose
    private String name;
    @SerializedName("annotations")
    @Expose
    private List<Annotation> annotations = new ArrayList<Annotation>();
    @SerializedName("genotypes")
    @Expose
    private List<String> genotypes = new ArrayList<String>();
    @SerializedName("strength")
    @Expose
    private Strength strength;

    /**
     * 
     * @return
     *     The objCls
     */
    public String getObjCls() {
        return objCls;
    }

    /**
     * 
     * @param objCls
     *     The objCls
     */
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
     *
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

    /**
     * 
     * @return
     *     The annotations
     */
    public List<Annotation> getAnnotations() {
        return annotations;
    }

    /**
     * 
     * @param annotations
     *     The annotations
     */
    public void setAnnotations(List<Annotation> annotations) {
        this.annotations = annotations;
    }

    /**
     * 
     * @return
     *     The genotypes
     */
    public List<String> getGenotypes() {
        return genotypes;
    }

    /**
     * 
     * @param genotypes
     *     The genotypes
     */
    public void setGenotypes(List<String> genotypes) {
        this.genotypes = genotypes;
    }

    /**
     * 
     * @return
     *     The strength
     */
    public Strength getStrength() {
        return strength;
    }

    /**
     * 
     * @param strength
     *     The strength
     */
    public void setStrength(Strength strength) {
        this.strength = strength;
    }

}
