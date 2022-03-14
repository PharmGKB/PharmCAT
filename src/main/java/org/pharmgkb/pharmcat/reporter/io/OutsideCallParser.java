package org.pharmgkb.pharmcat.reporter.io;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.Arrays;
import java.util.List;
import java.util.function.Predicate;
import java.util.stream.Collectors;
import com.google.common.base.Preconditions;
import org.apache.commons.lang3.StringUtils;
import org.pharmgkb.pharmcat.reporter.model.OutsideCall;


/**
 * Class to parse outside-call text files into {@link OutsideCall} objects. This was created with Astrolabe in mind but 
 * can be used for any outside calls. This parser will ignore any lines starting with "#" or blank lines and interpret 
 * the rest as TSV formatted lines with the first field as gene symbol, the second field as diplotype call, and the
 * third field as the phenotype call.
 *
 * The diplotype call in the second "column" can optionally include the gene symbol as a prefix on the allele name, but
 * it will be ignored. The gene is read from the first "column". As an example, this means both
 * <code>CYP2D6*1/CYP2D6*2</code> and <code>*1/*2</code> values in the second "column" are equivalent.
 *
 * The phenotype call in the third "column" is not required if the diplotype is supplied, the phenotyper will fill it
 * in. You can supply a phenotype without supplying a diplotype.
 *
 * The diplotype call must include 2 alleles and they must be separated by a "/".
 *
 * Multiple diplotypes can be separated by " or ". For example, "*1/*2 or *2/*3".
 *
 * @author Ryan Whaley
 */
public class OutsideCallParser {
  private static final Predicate<String> sf_nonCommentLine = (l) -> !l.startsWith("#");

  public static List<OutsideCall> parse(Path filePath) throws IOException {
    Preconditions.checkNotNull(filePath);

    return Files.lines(filePath)
        .filter(sf_nonCommentLine)
        .map(OutsideCall::new)
        .collect(Collectors.toList());
  }

  public static List<OutsideCall> parse(String outsideCallData) {
    return Arrays.stream(StringUtils.stripToEmpty(outsideCallData).split("\n"))
        .filter(sf_nonCommentLine)
        .map(OutsideCall::new)
        .collect(Collectors.toList());
  }
}
