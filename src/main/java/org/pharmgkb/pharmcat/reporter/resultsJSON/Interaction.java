package org.pharmgkb.pharmcat.reporter.resultsJSON;

import java.util.List;
import org.pharmgkb.pharmcat.reporter.model.CPICinteraction;
import org.pharmgkb.pharmcat.reporter.model.Group;
import org.pharmgkb.pharmcat.reporter.model.RelatedChemical;


public class Interaction {
    private List<Group> groupList;
    private List<RelatedChemical> relatedChemicals;
    private String source;
    private String summaryHtml;
    private String textHtml;
    
    public Interaction( CPICinteraction inter){
        this.relatedChemicals = inter.getRelatedChemicals();
        this.source = inter.getSource();
        this.summaryHtml = inter.getSummaryHtml();
        this.textHtml = inter.getTextHtml();        
    }

    public void addAllToGroup(List<Group> gList ){
        groupList.addAll( gList );
    }
}
