package org.pharmgkb.pharmcat.phenotype;

import java.io.BufferedReader;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.List;
import java.util.function.Predicate;
import com.google.common.base.Charsets;
import com.google.common.base.Preconditions;
import org.apache.commons.lang3.StringUtils;
import org.pharmgkb.pharmcat.phenotype.model.OutsideCall;


/**
 * Class to parse outside-call text files into {@link OutsideCall} objects.
 * This parser will ignore any lines starting with "#" or blank lines and interpret the rest as TSV formatted lines as:
 * <ol>
 *   <li>gene symbol</li>
 *   <li>diplotype</li>
 *   <li>phenotype</li>
 *   <li>activity score</li>
 * </ol>with
 * <p>
 * The diplotype call in the second "column" can optionally include the gene symbol as a prefix on the allele name, but
 * it will be ignored. The gene is read from the first "column". As an example, this means both
 * <code>CYP2D6*1/CYP2D6*2</code> and <code>*1/*2</code> values in the second "column" are equivalent.
 * <p>
 * The phenotype call in the third "column" is not required if the diplotype is supplied, the phenotyper will fill it
 * in. You can supply a phenotype without supplying a diplotype.
 * <p>
 * The diplotype call must include 2 alleles, and they must be separated by a "/".
 * <p>
 * Multiple diplotypes can be separated by " or ". For example, "*1/*2 or *2/*3".
 *
 * @author Ryan Whaley
 */
public class OutsideCallParser {
  private static final Predicate<String> sf_nonCommentLine = (l) -> StringUtils.isNotBlank(l) && !l.startsWith("#");

  public static List<OutsideCall> parse(Path filePath) throws IOException {
    Preconditions.checkNotNull(filePath);

    List<OutsideCall> calls = new ArrayList<>();
    try (BufferedReader reader = Files.newBufferedReader(filePath, Charsets.UTF_8)) {
      int x = 0;
      String line;
      while ((line = reader.readLine()) != null) {
        x += 1;
        if (sf_nonCommentLine.test(line)) {
          calls.add(new OutsideCall(line, x));
        }
      }
    }
    return calls;
  }

  public static List<OutsideCall> parse(String outsideCallData) {
    List<OutsideCall> calls = new ArrayList<>();
    String[] lines = StringUtils.stripToEmpty(outsideCallData).split("\n");
    for (int x = 0; x < lines.length; x += 1) {
      String line = lines[x];
      if (sf_nonCommentLine.test(line)) {
        calls.add(new OutsideCall(line, x + 1));
      }
    }
    return calls;
  }
}
