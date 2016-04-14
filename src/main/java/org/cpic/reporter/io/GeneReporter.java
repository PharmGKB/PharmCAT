package org.cpic.reporter.io;

import java.util.ArrayList;
import java.util.List;

import org.cpic.reporter.model.CPICException;
import org.cpic.reporter.model.HaplotypeCallerMultiGeneJSON.DiplotypeCall;

public class GeneReporter {
    
    List<CPICException> exceptList = new ArrayList<CPICException>();
    
    public void addToExceptList( CPICException except ){
        exceptList.add(except);
    }

}
