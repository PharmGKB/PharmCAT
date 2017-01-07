package org.pharmgkb.pharmcat.annotation;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.Collection;
import java.util.List;
import javax.annotation.Nonnull;
import com.google.common.base.Preconditions;
import com.google.common.base.Splitter;
import com.google.common.collect.HashMultimap;
import com.google.common.collect.Multimap;
import org.pharmgkb.pharmcat.ParseException;
import org.pharmgkb.pharmcat.annotation.model.RsidAnnotation;


/**
 * This parses annotation files.
 *
 * @author Mark Woon
 */
public class AnnotationReader {
  public static final String RSID_ANNOTATIONS_FILENAME = "rsid_annotations.tsv";
  public static final String NAMED_ALLELE_ANNOTATIONS_FILENAME = "namedAllele_annotations.tsv";
  private static final Splitter sf_tsvSplitter = Splitter.on("\t").trimResults();
  private Multimap<String, RsidAnnotation> m_rsidAnnotationByGene = HashMultimap.create();


  public @Nonnull Collection<RsidAnnotation> getRsidAnnotations(@Nonnull String gene) {
    return m_rsidAnnotationByGene.get(gene);
  }


  public void read(@Nonnull Path path) throws IOException {

    Preconditions.checkNotNull(path);
    if (!Files.exists(path)) {
      return;
    }

    if (Files.isDirectory(path)) {
      parseRsidAnnotations(path.resolve(RSID_ANNOTATIONS_FILENAME));
      parseNamedAlleleAnnotations(path.resolve(NAMED_ALLELE_ANNOTATIONS_FILENAME));

    } else {
      if (path.getFileName().toString().equals(RSID_ANNOTATIONS_FILENAME)) {
        parseRsidAnnotations(path);
      } else if (path.getFileName().toString().equals(NAMED_ALLELE_ANNOTATIONS_FILENAME)) {
        parseNamedAlleleAnnotations(path);
      }
    }
  }


  private void parseNamedAlleleAnnotations(@Nonnull Path file) {

    Preconditions.checkNotNull(file);
    if (!Files.exists(file)) {
      return;
    }
    Preconditions.checkArgument(Files.isRegularFile(file));

    // TODO: finish this
  }

  private void parseRsidAnnotations(@Nonnull Path file) {

    Preconditions.checkNotNull(file);
    if (!Files.exists(file)) {
      return;
    }
    Preconditions.checkArgument(Files.isRegularFile(file));

    try {
      for (String line : Files.readAllLines(file)) {
        List<String> data = sf_tsvSplitter.splitToList(line);
        String gene = data.get(1);
        String rsid = data.get(2);
        m_rsidAnnotationByGene.put(gene, new RsidAnnotation(gene, rsid));
      }
    } catch (IOException ex) {
      throw new ParseException("Couldn't read lines in " + file, ex);
    }
  }
}
