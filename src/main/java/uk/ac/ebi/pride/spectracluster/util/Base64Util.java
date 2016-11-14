/*
 $Id: Base64Util.java 3402 2009-08-26 07:26:13Z fredrik $

 Copyright (C) 2006, 2007 Gregory Vincic

 This file is part of Proteios.
 Available at http://www.proteios.org/

 Proteios is free software; you can redistribute it and/or modify it
 under the terms of the GNU General Public License as published by
 the Free Software Foundation; either version 2 of the License, or
 (at your option) any later version.

 Proteios is distributed in the hope that it will be useful, but
 WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with this program; if not, write to the Free Software
 Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
 02111-1307, USA.
 */
package uk.ac.ebi.pride.spectracluster.util;

import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.util.ArrayList;
import java.util.List;
import java.util.zip.DataFormatException;
import java.util.zip.Inflater;

/**
 * This class supports conversion to/from Base64-coded data.
 * Used e.g. to encode the data (peak) part of mass spectra
 * in an mzData file.
 *
 * @author Olle
 * @version 2.0
 */
public class Base64Util {
    /**
     * Decodes Base64 coded data String
     * and returns an ArrayList&lt;double&gt;
     * containing the extracted data.
     * Uses class uk.ac.ebi.pride.spectracluster.util.Base64Coder to decode data.
     *
     * @param doublePrecision boolean flag, true if precision == "64".
     * @param bigEndian       boolean flag indicating endian type.
     * @param dataString      String with Base64 coded data.
     * @return ArrayList&lt;double&gt; with extracted values.
     */
    public static List<Double> decode(boolean doublePrecision, boolean bigEndian, String dataString) {
        /*
         * Decode the data and put it into byte array
		 */
        boolean zLibCompression = false;
        /*
		 * Return result as list
		 */
        return decode(doublePrecision, bigEndian, zLibCompression, dataString);
    }


    /**
     * Decodes Base64 coded data String
     * and returns an ArrayList&lt;double&gt;
     * containing the extracted data.
     * Uses class uk.ac.ebi.pride.spectracluster.util.Base64Coder to decode data.
     *
     * @param doublePrecision boolean flag, true if precision == "64".
     * @param bigEndian       boolean flag indicating endian type.
     * @param zLibCompression boolean flag indicating zLib compression.
     * @param dataString      String with Base64 coded data.
     * @return ArrayList&lt;double&gt; with extracted values.
     */
    public static List<Double> decode(boolean doublePrecision, boolean bigEndian, boolean zLibCompression, String dataString) {
		/*
		 * Decode the data and put it into byte array
		 */
        byte[] dataByteArray = Base64Coder.decode(dataString.toCharArray());
        if (zLibCompression) {
            // Decompress the bytes
            Inflater decompresser = new Inflater();
            decompresser.setInput(dataByteArray, 0, dataByteArray.length);
            byte[] uncompressedByteArray = new byte[3 * dataByteArray.length];
            int uncompressedArrayLength = 0;
            try {
                uncompressedArrayLength = decompresser.inflate(uncompressedByteArray);
            } catch (DataFormatException e) {
                throw new IllegalStateException("DataFormatException when decompressing byte array", e);
            }
            decompresser.end();
            // Use uncompressed byte array
            dataByteArray = new byte[uncompressedArrayLength];
            System.arraycopy(uncompressedByteArray, 0, dataByteArray, 0, uncompressedArrayLength);
        }
        ByteBuffer dataByteBuffer = ByteBuffer.wrap(dataByteArray);
		/*
		 * Get size of result list from
		 * byte array size and precision.
		 */
        int arraySize = dataByteArray.length;
        int bytesPerValue = 4;
        if (doublePrecision) {
			/*
			 * Precision == "64"
			 */
            bytesPerValue = 8;
        }
        int dataLength = arraySize / bytesPerValue;
		/*
		 * Create result arrayList
		 */
        List<Double> resultArray = new ArrayList<Double>(1);
		/*
		 * Java works with big endian by default.
		 * If the data is little endian, the
		 * byte order is changed.
		 *
		 * This configuration influences how data
		 * is retrieved using methods like
		 * getFloat() and getDouble(),
		 */
        if (!bigEndian) {
            dataByteBuffer.order(ByteOrder.LITTLE_ENDIAN);
        }
		/*
		 * Put the data into the result list
		 */
        if (!doublePrecision) {
			/*
			 * Precision == "32"
			 */
            for (int i = 0; i < dataLength; i++) {
                float value = dataByteBuffer.getFloat();
                resultArray.add(i, new Float(value).doubleValue());
            }
        } else {
			/*
			 * Precision == "64"
			 */
            for (int i = 0; i < dataLength; i++) {
                double value = dataByteBuffer.getDouble();
                resultArray.add(i, value);
            }
        }
		/*
		 * Return result as list
		 */
        return resultArray;
    }


    /**
     * Encodes Base64 coded data List&lt;double&gt;
     * and returns a String
     * containing the encoded data.
     * Uses class uk.ac.ebi.pride.spectracluster.util.Base64Coder to encode data.
     *
     * @param doublePrecision boolean flag, true if precision == "64".
     * @param bigEndian       boolean flag indicating endian type.
     * @param inDataList      ArrayList with double or integer (32bit) data values.
     * @return String with encoded values.
     */
    public static String encode(boolean doublePrecision, boolean bigEndian, List<? extends Number> inDataList) {
		/*
		 * Get size of byte array from
		 * input list size and precision.
		 */
        int dataLength = inDataList.size();
        int bytesPerValue = 4;
        if (doublePrecision) {
			/*
			 * Precision == "64"
			 */
            bytesPerValue = 8;
        }
        int arraySize = bytesPerValue * dataLength;
		/*
		 * Create dataByteBuffer ByteBuffer
		 */
        byte[] dataByteArray = new byte[arraySize];
        ByteBuffer dataByteBuffer = ByteBuffer.wrap(dataByteArray);
		/*
		 * Java works with big endian by default.
		 * If the data is little endian, the
		 * byte order is changed.
		 *
		 * Note! This configuration must be set before
		 * the ByteBuffer is loaded with data using
		 * methods like putFloat() and putDouble(),
		 * as it influences how the data will be stored.
		 */
        if (!bigEndian) {
            dataByteBuffer.order(ByteOrder.LITTLE_ENDIAN);
        }
		/*
		 * Put the data into the byte buffer
		 */
        if (!doublePrecision) {
			/*
			 * Precision == "32"
			 */
            for (Number anInDataList : inDataList) {
                if (anInDataList instanceof Double) {
                    Double doubleVal = (Double) anInDataList;
                    dataByteBuffer.putFloat(doubleVal.floatValue());
                } else if (anInDataList instanceof Integer) {
                    Integer intVal = (Integer) anInDataList;
                    dataByteBuffer.putInt(intVal);
                }
            }
        } else {
			/*
			 * Precision == "64"
			 */
            for (Number anInDataList : inDataList) {
                Double doubleVal = (Double) anInDataList;
                dataByteBuffer.putDouble(doubleVal);
            }
        }
		/*
		 * Encode the data and put it into char array
		 */
        char[] dataCharArray = Base64Coder.encode(dataByteBuffer.array());
		/*
		 * Return result as String
		 */
        return new String(dataCharArray);
    }

}
