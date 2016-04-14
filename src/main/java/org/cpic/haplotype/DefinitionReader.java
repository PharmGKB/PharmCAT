package org.cpic.haplotype;

import java.io.BufferedReader;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Map;
import javax.annotation.Nonnull;
import com.google.common.base.Preconditions;
import com.google.common.collect.ArrayListMultimap;
import com.google.common.collect.ListMultimap;

/**
 * Read in haplotype definition files.
 *
 * @author Mark Woon
 */
public class DefinitionReader {
  private ListMultimap<String, Variant> m_haplotypePositions  = ArrayListMultimap.create();
  private ListMultimap<String, Haplotype> m_haplotypes = ArrayListMultimap.create();
  private Map<String, TSVfile> m_files = new HashMap<>();


  public Map<String, TSVfile> getFiles() {
	  return m_files;
  }


  public ListMultimap<String, Variant> getHaplotypePositions() {
    return m_haplotypePositions;
  }

  public ListMultimap<String, Haplotype> getHaplotypes() {
	    return m_haplotypes;
  }



  public void read(Path path) throws IOException {

    if (Files.isDirectory(path)) {
      Files.list(path)
          .filter(f -> f.toString().endsWith(".tsv"))
          .forEach(this::readFile);
    } else {
      readFile(path);
    }
  }






  private void readFile(@Nonnull Path file)  {

    Preconditions.checkNotNull(file);
    Preconditions.checkArgument(Files.isRegularFile(file));
    System.out.println(file);

    try (BufferedReader bufferedReader = Files.newBufferedReader(file)) {

    	TSVfile inProccessFile = new TSVfile(file.toString());
        ArrayList<Variant> variants = new ArrayList<>();
        ArrayList<Haplotype> haplotypes = new ArrayList<>();

        String line;
        boolean NormalFunctionAllele = true;

    	while((line = bufferedReader.readLine()) != null) {
    		//System.out.println(line);

            String[] fields = line.split("\t");
            for (int i = 0; i < fields.length; i++){
            	fields[i] = fields[i].trim();
            }

            if (fields[0].equals("FormatVersion")){
            	if (fields.length>0){
            		inProccessFile.setFormatVersion(fields[1]);
            	}else{
            		inProccessFile.setFormatVersion("");
            	}
            	
            }
            else if (fields[0].equals("GeneName")){
            	if (fields.length>0){
            		inProccessFile.setGeneName(fields[1]);
            	}else{
            		inProccessFile.setGeneName("");
            	}
            	inProccessFile.setGeneName(fields[1]);
            }
            else if (fields[0].equals("GeneRefSeq")){
            	inProccessFile.setGeneID(fields[1]);
            }
            else if (fields[0].equals("GeneOrientation")){
            	if (fields.length>0){
            		inProccessFile.setGeneOrientation(fields[1]);
            	}else{
            		inProccessFile.setGeneOrientation("");
            	}
            }
            else if (fields[0].equals("ContentDate")){
            	inProccessFile.setContentDate(fields[1]);
            }
            else if (fields[0].equals("ContentVersion")){
            	inProccessFile.setContentVersion(fields[1]);
            }
            else if (fields[0].equals("GenomeBuild")){
            	inProccessFile.setGenomeBuild(fields[1]);
            }
            else if (fields[0].equals("ChrName")){
            	inProccessFile.setChromosome(fields[1]);
            }
            else if (fields[0].equals("ChrRefSeq")){
            	inProccessFile.setChromosomeID(fields[1]);
            }
            else if (fields[0].equals("ProteinRefSeq")){
            	inProccessFile.setProteinID(fields[1]);
            }
            else if (fields[0].equals("ResourceNote")){
            	for (int i = 4; i < fields.length; i++){
                	Variant newVariant = new Variant(inProccessFile.getChromosome(),inProccessFile.getGeneName(),fields[i]);
                	variants.add(newVariant);
                }
            }
            else if (fields[0].equals("ProteinNote")){
            	for (int i = 4; i < fields.length; i++){
                	variants.get(i-4).addProteingEffect(fields[i]);
                }
            }
            else if (fields[0].equals("ChrPosition")){
            	for (int i = 4; i < fields.length; i++){
                	variants.get(i-4).addProteingEffect(fields[i]);
                	variants.get(i-4).setStartPOS();
                }
            }
            else if (fields[0].equals("GenePosition")){
            	for (int i = 4; i < fields.length; i++){
                	variants.get(i-4).addProteingEffect(fields[i]);
                }
            }
            else if (fields[0].equals("rsID")){
            	for (int i = 4; i < fields.length; i++){
                	variants.get(i-4).set_rsID(fields[i]);
                }

            }
            else if (fields[0].equals("Allele")){
            	ArrayList <String> alleles = new ArrayList<>();
            	ArrayList <Variant> hapVariants = new ArrayList<>();

            	int forLength = variants.size()+4;
            	if (forLength>fields.length){
            		forLength=fields.length;
            	}


            	for (int i = 4; i < forLength; i++){

                	if (NormalFunctionAllele){
                		variants.get(i-4).setREF(fields[i]);
                		alleles.add(fields[i]);
                		hapVariants.add(variants.get(i-4));

                	}
                	else{
                		if (!fields[i].trim().equals("")){
	                		variants.get(i-4).addALT(fields[i]);
	                		alleles.add(fields[i]);
	                		hapVariants.add(variants.get(i-4));
                		}
                	}

            	}
            	NormalFunctionAllele=false;
            	haplotypes.add(new Haplotype(hapVariants,fields[1],fields[2],fields[3],alleles));

            }
        }
      for (int i = 0; i < haplotypes.size(); i++){
        m_haplotypes.put(inProccessFile.getGeneName(),haplotypes.get(i));
      }
      for (int i = 0; i < variants.size(); i++){
        m_haplotypePositions.put(inProccessFile.getGeneName(),variants.get(i));
      }

	m_files.put(inProccessFile.getGeneName(), inProccessFile);



    } catch (Exception ex) {
	throw new RuntimeException("Failed to parse " + file, ex);
    }
  }

  public static void main(String[] args) {

    try {
      Path path = Paths.get(args[0]);
      DefinitionReader r = new DefinitionReader();
      r.read(path);
    } catch (Exception ex) {
      ex.printStackTrace();
    }
  }
}
