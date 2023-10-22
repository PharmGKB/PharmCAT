package org.pharmgkb.pharmcat.util;

import java.lang.reflect.Type;
import java.time.Instant;
import java.time.format.DateTimeParseException;
import java.util.Date;
import com.google.gson.JsonDeserializationContext;
import com.google.gson.JsonDeserializer;
import com.google.gson.JsonElement;
import com.google.gson.JsonParseException;
import com.google.gson.JsonPrimitive;
import com.google.gson.JsonSerializationContext;
import com.google.gson.JsonSerializer;


/**
 * This is a GSON adapter for deserializing/serializing ISO 8601 strings to/from {@link Date}s.
 * <p>
 * This is necessary until <a href="https://github.com/google/gson/issues/2472">GSON-2472</a> gets fixed.
 *
 * @author Mark Woon
 */
public class GsonDateAdapter implements JsonSerializer<Date>, JsonDeserializer<Date> {

  @Override
  public Date deserialize(JsonElement json, Type typeOfT, JsonDeserializationContext context) throws JsonParseException {
    try {
      return Date.from(Instant.parse(json.getAsString()));
    } catch (DateTimeParseException ex) {
      throw new JsonParseException(ex);
    }
  }

  @Override
  public JsonElement serialize(Date src, Type typeOfSrc, JsonSerializationContext context) {
    return new JsonPrimitive(src.toInstant().toString());
  }
}
