package org.pharmgkb.pharmcat.reporter.resultsJSON;

import java.util.List;
import org.pharmgkb.pharmcat.reporter.model.CPICinteraction;
import org.pharmgkb.pharmcat.reporter.model.Group;
import org.pharmgkb.pharmcat.reporter.model.RelatedChemical;


public class Interaction {
  
	private String name;
	
	private List<Group> groupList;
    private List<RelatedChemical> relatedChemicals;
    private String source;
    private String summaryHtml;
    private String textHtml;
    
    public Interaction(){
    }
    
    public Interaction( CPICinteraction inter){
    	this.name = inter.getName();
        this.relatedChemicals = inter.getRelatedChemicals();
        this.source = inter.getSource();
        this.summaryHtml = inter.getSummaryHtml();
        this.textHtml = inter.getTextHtml();        
    }

    public void addToGroup(Group gList ){
        groupList.add( gList );
    }

	public String getName() {
		return name;
	}

	public void setName(String name) {
		this.name = name;
	}

	public List<Group> getGroupList() {
		return groupList;
	}

	public void setGroupList(List<Group> groupList) {
		this.groupList = groupList;
	}

	public List<RelatedChemical> getRelatedChemicals() {
		return relatedChemicals;
	}

	public void setRelatedChemicals(List<RelatedChemical> relatedChemicals) {
		this.relatedChemicals = relatedChemicals;
	}

	public String getSource() {
		return source;
	}

	public void setSource(String source) {
		this.source = source;
	}

	public String getSummaryHtml() {
		return summaryHtml;
	}

	public void setSummaryHtml(String summaryHtml) {
		this.summaryHtml = summaryHtml;
	}

	public String getTextHtml() {
		return textHtml;
	}

	public void setTextHtml(String textHtml) {
		this.textHtml = textHtml;
	}
    
}
