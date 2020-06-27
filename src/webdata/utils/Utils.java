package webdata.utils;

import java.math.BigInteger;
import java.nio.ByteBuffer;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashSet;

/**
 * Static class for utility methods
 */
public final class Utils {

    /**
     * Empty and private constructor to make this class static.
     */
    private Utils() {}

    /**
     * Returns a Byte ArrayList of exactly 4 bytes.
     * @param arr Byte ArrayList to pad
     * @return Byte ArrayList of size 4
     */
    public static ArrayList<Byte> padByte(ArrayList<Byte> arr) {
        int pad = 4 - arr.size();
        int i = 0;
        Byte zero = 0;
        while (i < pad) {
            arr.add(0, zero);
            ++i;
        }
        return arr;
    }

    /**
     * Returns a byte array of exactly 4 bytes.
     * @param arr byte array to pad
     * @return byte array of size 4
     */
    public static byte[] padByte(byte[] arr) {
        byte[] newArr = new byte[4];
        int pad = 4 - arr.length;
        for (int i = 0; i < arr.length; ++i) {
            newArr[i + pad] = arr[i];
        }
        return newArr;
    }

    /**
     * Convert the given Integer to it's byte representation
     * @param i Integer to convert
     * @return Byte ArrayList of the Integer's representation
     */
    public static ArrayList<Byte> intToByte(final Integer i) {
        BigInteger bi = BigInteger.valueOf(i);
        ArrayList<Byte> bigByte = new ArrayList<>();
        for (byte b: bi.toByteArray()) {
            bigByte.add(b);
        }
        return bigByte;
    }

    /**
     * Convert a byte array representing an integer to it's corresponding int
     * @param intBytes The byte array
     * @return The corresponding int
     */
    public static int byteArrayToInt(byte[] intBytes){
        ByteBuffer byteBuffer = ByteBuffer.wrap(intBytes);
        return byteBuffer.getInt();
    }

    /**
     * Convert an ArrayList of Short to short array
     * @param list ArrayList of Short
     * @param arr Array to populate
     */
    public static void toPrimitiveArray(ArrayList<Short> list, short[] arr) {
        for (int i = 0; i < list.size(); ++i) {
            arr[i] = list.get(i);
        }
    }

    /**
     * Convert an ArrayList of Byte to byte array
     * @param list ArrayList of Byte
     * @param arr Array to populate
     */
    public static void toPrimitiveArray(ArrayList<Byte> list, byte[] arr) {
        for (int i = 0; i < list.size(); ++i) {
            arr[i] = list.get(i);
        }
    }

    public static double[] integerCollectionToDoubleArray(Collection<Integer> values) {
        double[] result = new double[values.size()];
        int i = 0;
        for (Integer val: values) {
            result[i] = (double) val;
            ++i;
        }
        return result;
    }

    public static double[] doubleCollectionToDoubleArray(Collection<Double> values) {
        double[] result = new double[values.size()];
        int i = 0;
        for (double val: values) {
            result[i] = val;
            ++i;
        }
        return result;
    }
}
