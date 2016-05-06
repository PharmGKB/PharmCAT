/*
 ----- BEGIN LICENSE BLOCK -----
 This Source Code Form is subject to the terms of the Mozilla Public License, v.2.0.
 If a copy of the MPL was not distributed with this file, You can obtain one at http://mozilla.org/MPL/2.0/.
 ----- END LICENSE BLOCK -----
 */
package org.pharmgkb.common.io;

import java.io.IOException;
import java.lang.invoke.MethodHandles;
import java.security.GeneralSecurityException;
import java.security.KeyFactory;
import java.security.PrivateKey;
import java.security.spec.PKCS8EncodedKeySpec;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import javax.annotation.Nonnull;
import javax.xml.bind.DatatypeConverter;
import com.google.api.client.googleapis.auth.oauth2.GoogleCredential;
import com.google.api.client.googleapis.javanet.GoogleNetHttpTransport;
import com.google.api.client.http.HttpTransport;
import com.google.api.client.json.JsonFactory;
import com.google.api.client.json.gson.GsonFactory;
import com.google.common.base.Preconditions;
import com.google.common.collect.Lists;
import org.slf4j.LoggerFactory;


/**
 * This helper class handles the prep for using Google API.
 *
 * @author Mark Woon
 */
public class GoogleApiHelper implements AutoCloseable {
  private static final Pattern sf_privateKeyPattern =
      Pattern.compile("-----BEGIN .*?PRIVATE KEY-----(.*)-----END .*?PRIVATE KEY-----");
  private HttpTransport m_httpTransport;
  private JsonFactory m_jsonFactory;
  private GoogleCredential m_credential;


  /**
   * Constructs a {@link GoogleCredential}.
   */
  public GoogleApiHelper(@Nonnull String userId, @Nonnull String privateKey, String... scopes)
      throws IOException, GeneralSecurityException {

    Preconditions.checkNotNull(userId);
    Preconditions.checkNotNull(privateKey);
    m_jsonFactory = GsonFactory.getDefaultInstance();
    m_httpTransport = GoogleNetHttpTransport.newTrustedTransport();
    m_credential = new GoogleCredential.Builder()
        .setTransport(m_httpTransport)
        .setJsonFactory(m_jsonFactory)
        .setServiceAccountId(userId)
        .setServiceAccountPrivateKey(getPrivateKey(privateKey))
        .setServiceAccountScopes(Lists.newArrayList(scopes))
        .build();
  }

  @Override
  public void close() {
    try {
      m_httpTransport.shutdown();
    } catch (IOException ex) {
      LoggerFactory.getLogger(MethodHandles.lookup().lookupClass()).error("Problem shutting down", ex);
    }
  }


  public JsonFactory getJsonFactory() {
    return m_jsonFactory;
  }


  public HttpTransport getHttpTransport() {
    return m_httpTransport;
  }


  public GoogleCredential getCredential() {
    return m_credential;
  }



  /**
   * This class converts a private key in PEM/PKCS#8 formatted String into a {@link PrivateKey}.
   */
  private static PrivateKey getPrivateKey(String key) throws GeneralSecurityException {

    key = stripKey(key);
    PKCS8EncodedKeySpec keySpec = new PKCS8EncodedKeySpec(DatatypeConverter.parseBase64Binary(key));
    KeyFactory keyFactory = KeyFactory.getInstance("RSA");
    return keyFactory.generatePrivate(keySpec);
  }


  /**
   * Strips header and footer from key and concatenates it to a single line.
   */
  static String stripKey(String key) {
    key = key.replaceAll("\\r|\\n", "");
    Matcher m = sf_privateKeyPattern.matcher(key);
    if (m.find()) {
      return m.group(1);
    }
    return key;
  }
}
