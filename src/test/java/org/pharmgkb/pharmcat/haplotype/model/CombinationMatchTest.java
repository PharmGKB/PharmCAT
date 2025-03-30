package org.pharmgkb.pharmcat.haplotype.model;

import java.util.List;
import java.util.SortedSet;
import java.util.TreeSet;
import java.util.stream.Collectors;
import org.junit.jupiter.api.Test;
import org.pharmgkb.pharmcat.definition.model.NamedAllele;
import org.pharmgkb.pharmcat.definition.model.VariantLocus;

import static org.junit.jupiter.api.Assertions.assertEquals;


/**
 * This is a JUnit test for {@link CombinationMatch}.
 *
 * @author Mark Woon
 */
class CombinationMatchTest {

  @Test
  void compareWithHaplotypeMatch() {

    VariantLocus[] refVariants = new VariantLocus[0];
    NamedAllele na1 = new NamedAllele("1", "*2", new String[0], new String[0], false);
    NamedAllele na2 = new NamedAllele("2", "*3", new String[0], new String[0], false);
    CombinationMatch cm1 = new CombinationMatch(refVariants, "123", List.of(na1, na2), null);
    HaplotypeMatch hapMatch = new HaplotypeMatch(na1);

    SortedSet<Object> set = new TreeSet<>();
    set.add(hapMatch);
    set.add(cm1);
    System.out.println(set.stream().map(Object::toString).collect(Collectors.joining(" / ")));
    assertEquals(2, set.size());
    assertEquals(hapMatch, set.first());

    set.clear();
    set.add(cm1);
    set.add(hapMatch);
    System.out.println(set.stream().map(Object::toString).collect(Collectors.joining(" / ")));
    assertEquals(2, set.size());
    assertEquals(hapMatch, set.first());


    CombinationMatch cm2 = new CombinationMatch(refVariants, "123", List.of(na1, na2), null);
    set.clear();
    set.add(cm1);
    set.add(cm2);
    System.out.println(set.stream().map(Object::toString).collect(Collectors.joining(" / ")));
    assertEquals(1, set.size());


    NamedAllele na3 = new NamedAllele("3", "*4", new String[0], new String[0], false);
    CombinationMatch cm3 = new CombinationMatch(refVariants, "123", List.of(na1, na2, na3), null);
    set.clear();
    set.add(cm1);
    set.add(cm3);
    System.out.println(set.stream().map(Object::toString).collect(Collectors.joining(" / ")));
    assertEquals(2, set.size());
    assertEquals(cm1, set.first());
  }
}
