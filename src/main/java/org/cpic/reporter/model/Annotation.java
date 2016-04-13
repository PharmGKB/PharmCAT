
package org.cpic.reporter.model;

import javax.annotation.Generated;
import com.google.gson.annotations.Expose;
import com.google.gson.annotations.SerializedName;


public class Annotation {

    @SerializedName("id")
    @Expose
    private Integer id;
    @SerializedName("text")
    @Expose
    private String text;
    @SerializedName("textHtml")
    @Expose
    private String textHtml;
    @SerializedName("type")
    @Expose
    private Type type;

    /**
     * 
     * @return
     *     The id
     */
    public Integer getId() {
        return id;
    }

    /**
     * 
     * @param id
     *     The id
     */
    public void setId(Integer id) {
        this.id = id;
    }

    /**
     * 
     * @return
     *     The text
     */
    public String getText() {
        return text;
    }

    /**
     * 
     * @param text
     *     The text
     */
    public void setText(String text) {
        this.text = text;
    }

    /**
     * 
     * @return
     *     The textHtml
     */
    public String getTextHtml() {
        return textHtml;
    }

    /**
     * 
     * @param textHtml
     *     The textHtml
     */
    public void setTextHtml(String textHtml) {
        this.textHtml = textHtml;
    }

    /**
     * 
     * @return
     *     The type
     */
    public Type getType() {
        return type;
    }

    /**
     * 
     * @param type
     *     The type
     */
    public void setType(Type type) {
        this.type = type;
    }

}
