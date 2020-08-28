package webdata.utils;

import java.math.BigInteger;
import java.nio.ByteBuffer;
import java.util.*;

/**
 * Static class for utilities
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


    ////////////////////////////////////////////////////////////////////////////////////////////////////////
    // **                        Damerau-Levenshtein Edit Distance Algorithm                           ** //
    ////////////////////////////////////////////////////////////////////////////////////////////////////////

    /**
     * Compute the Damerau-Levenshtein distance between the specified source
     * string and the specified target string.
     * @param s1 First String, Considered as the "correct" word
     * @param s2 Second String, Considered as the "misspelled" word
     * @return
     */
    public static DLDistance DLD(String s1, String s2) {
        if (s1 == null || s2 == null) {  // Invalid input
            return null;
        }

        int[][] distances = new int[s1.length() + 1][s2.length() + 1];

        // If there exists a distance to calculate
        if (!s1.equals(s2)) {
            for (int i = 0; i < s1.length() + 1; ++i) {
                distances[i][0] = i;
            }
            for (int j = 0; j < s2.length() + 1; ++j) {
                distances[0][j] = j;
            }

            for (int i = 0; i < s1.length(); ++i) {
                for (int j = 0; j < s2.length(); ++j) {
                    int cost = (s1.charAt(i) == s2.charAt(j)) ? 0 : 1;

                    distances[i + 1][j + 1] = min(distances[i][j + 1] + 1,  // Insertion
                            distances[i + 1][j] + 1,  // Deletion
                            distances[i][j] + cost  // Substitution
                    );
                    // Transposition
                    if (i > 0 && j > 0 && s1.charAt(i) == s2.charAt(j - 1) && s1.charAt(i - 1) == s2.charAt(j)) {
                        distances[i + 1][j + 1] = min(distances[i+1][j+1], distances[i - 1][j - 1] + cost);
                    }
                }
            }
        }

        return new DLDistance(s2, s1, distances);
    }

    /**
     * Get the minimum of given int values
     * @param args Integers to compare
     * @return The minimum value
     */
    public static int min(int... args) {
        int min = Integer.MAX_VALUE;
        for (int arg : args) {
            if (arg < min) {
                min = arg;
            }
        }
        return min;
    }

    /**
     * Find all N-grams for the given term.
     * @param N -gram
     * @param edgeMark Marker for the edge of the word.
     * @param term The term to find n-grams for
     * @return
     */
    public static String[] findNGrams(int N, String edgeMark, String term) {
        // Create the n-grams array with size as the possible amount of n-grams for this term length.
        String[] ngrams = new String[term.length() + 2 - (N - 1)];

        // Make the term of the form $term$ (if edgeMark = "$").
        String newTerm = edgeMark.concat(term).concat(edgeMark);
        for (int j = 0; j < newTerm.length() - N + 1; ++j) {
            ngrams[j] = newTerm.substring(j, j + N);
        }

        return ngrams;
    }

    /**
     * Return the union size of 2 String arrays
     * @param a first array
     * @param b second array
     * @return Array union size
     */
    public static int unionSizeStringArrays(String[] a, String[] b) {
        Set<String> setA = new HashSet<>(Arrays.asList(a));
        Set<String> setB = new HashSet<>(Arrays.asList(b));
        setA.addAll(setB);
        return setA.size();
    }

    /**
     * Return the intersection size of 2 String arrays
     * @param a first array
     * @param b second array
     * @return Array intersection size
     */
    public static int intersectSizeStringArrays(String[] a, String[] b) {
        Set<String> setA = new HashSet<>(Arrays.asList(a));
        Set<String> setB = new HashSet<>(Arrays.asList(b));
        setA.retainAll(setB);
        return setA.size();
    }
}
