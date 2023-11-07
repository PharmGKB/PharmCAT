package org.pharmgkb.pharmcat.phenotype.model;

import org.junit.jupiter.api.Test;
import org.pharmgkb.pharmcat.reporter.BadOutsideCallException;

import static org.junit.jupiter.api.Assertions.assertThrows;


/**
 * This is a JUnit test for {@link OutsideCall}.
 *
 * @author Mark Woon
 */
class OutsideCallTest {


  @Test
  void testBad() {
    assertThrows(BadOutsideCallException.class, () -> {
      try {
        new OutsideCall("CYP2D6", 1);
      } catch (BadOutsideCallException ex) {
        System.out.println(ex.getMessage());
        throw ex;
      }
    });

    assertThrows(BadOutsideCallException.class, () -> {
      try {
        new OutsideCall("CYP2D6\t", 2);
      } catch (BadOutsideCallException ex) {
        System.out.println(ex.getMessage());
        throw ex;
      }
    });

    assertThrows(BadOutsideCallException.class, () -> {
      try {
        new OutsideCall("CYP2D6\t   ", 3);
      } catch (BadOutsideCallException ex) {
        System.out.println(ex.getMessage());
        throw ex;
      }
    });

    assertThrows(BadOutsideCallException.class, () -> {
      try {
        new OutsideCall("CYP2D6\tCYP2D6 *1/CYP2D6 *2/CYP2D6 *3", 4);
      } catch (BadOutsideCallException ex) {
        System.out.println(ex.getMessage());
        throw ex;
      }
    });
  }

  @Test
  void testGood() {
    new OutsideCall("CYP2D6\t*1", 1);
    new OutsideCall("CYP2D6\t*1\t",2);
    new OutsideCall("CYP2D6\t*1\t\t",3);
    new OutsideCall("CYP2D6\t\tIM", 4);
    new OutsideCall("CYP2D6\t\t\t2.0", 5);
    new OutsideCall("CYP2D6\t*1\tIM\t", 6);
    new OutsideCall("CYP2D6\t*1\tIM\t2.0", 7);
    new OutsideCall("CYP2D6\t*1\t\t2.0", 8);
    new OutsideCall("CYP2D6\t\tIM\2.0", 9);
  }
}
