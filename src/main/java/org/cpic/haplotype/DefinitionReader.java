package org.cpic.haplotype;

import com.google.common.base.Preconditions;
import com.google.common.collect.ArrayListMultimap;
import com.google.common.collect.ListMultimap;

import java.io.BufferedReader;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.Set;

/**
 * Read in haplotype definition files.
 *
 * @author Mark Woon
 */
public class DefinitionReader {
  private ListMultimap<String, Variant> m_haplotypePositions  = ArrayListMultimap.create();
  private ListMultimap<String, Haplotype> m_haplotypes = ArrayListMultimap.create();

  private Set<TSVfile> s_files = new HashSet<>();


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






  private void readFile(Path file)  {

    Preconditions.checkArgument(Files.isRegularFile(file));
    System.out.println(file);

    

    try (BufferedReader bufferedReader = Files.newBufferedReader(file)) {
    	
    	TSVfile inProccessFile = new TSVfile(file.toString());
        ArrayList<Variant> variants = new ArrayList<>();
        ArrayList<Haplotype> haplotypes = new ArrayList<>();

        String line = null;
        boolean NormalFunctionAllele = true;

    	while((line = bufferedReader.readLine()) != null) {
    		//System.out.println(line);

            String[] fields = line.split("\t");
            for (int i = 0; i < fields.length; i++){
            	fields[i] = fields[i].trim();
            }

            if (fields[0].equals("FormatVersion")){
            	inProccessFile.setFormatVersion(fields[1]);
            }
            else if (fields[0].equals("Gene")){
            	inProccessFile.setGeneID(fields[0]);
            	inProccessFile.setGeneName(fields[1]);
            }
            else if (fields[0].equals("ContentVersion")){
            	inProccessFile.setContentVersion(fields[1]);
            	inProccessFile.setContentDate(fields[1]);
            }
            else if (fields[0].equals("GenomeBuild")){
            	inProccessFile.setGenomeBuild(fields[1]);
            }
            else if (fields[0].equals("Chomosome")){
            	inProccessFile.setChromosomeID(fields[1]);
            	inProccessFile.setChromosome(fields[2]);
            }
            else if (fields[0].equals("Protein")){
            	inProccessFile.setProteinID(fields[1]);
            }
            else if (fields[0].equals("cDNAChange")){
            	for (int i = 4; i < fields.length; i++){
                	Variant newVariant = new Variant(inProccessFile.getChromosome(),inProccessFile.getGeneName(),fields[i]);
                	variants.add(newVariant);
                }
            }
            else if (fields[0].equals("ProteinEffect")){
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
                		variants.get(i-4).addALT(fields[i]);
                		alleles.add(fields[i]);
                		hapVariants.add(variants.get(i-4));
                	}
                	
            	}
            	NormalFunctionAllele=false;
            	haplotypes.add(new Haplotype(hapVariants,fields[1],fields[2],fields[3],alleles));
            	
            }
        }
    	for (int i = 0; i < haplotypes.size(); i++){
    		m_haplotypes.put(inProccessFile.getGeneName(),haplotypes.get(i));    		
        }
    	for (int i = 0; i < haplotypes.size(); i++){
    		m_haplotypePositions.put(inProccessFile.getGeneName(),variants.get(i));
        }
    	
    	s_files.add(inProccessFile);



    } catch (Exception ex) {
    	throw new RuntimeException(ex);
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
