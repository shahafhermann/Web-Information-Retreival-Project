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

        if (s1.equals(s2)) {  // No distance to compute
            int[][] empty = {};
            return new DLDistance(s2, s1, empty);
        }

        // The max possible distance
        int inf = s1.length() + s2.length();

        // Create and initialize the character array indices
        HashMap<Character, Integer> da = new HashMap<>();
        for (int i = 0; i < s1.length(); ++i) {
            da.put(s1.charAt(i), 0);
        }
        for (int j = 0; j < s2.length(); ++j) {
            da.put(s2.charAt(j), 0);
        }

        // Create the distance matrix H[0 .. s1.length+1][0 .. s2.length+1]
        int[][] distances = new int[s1.length() + 2][s2.length() + 2];

        // initialize the left and top edges of H
        for (int i = 0; i <= s1.length(); ++i) {
            distances[i + 1][0] = inf;
            distances[i + 1][1] = i;
        }

        for (int j = 0; j <= s2.length(); ++j) {
            distances[0][j + 1] = inf;
            distances[1][j + 1] = j;

        }

        // fill in the distance matrix H
        // look at each character in s1
        for (int i = 1; i <= s1.length(); ++i) {
            int db = 0;

            // look at each character in s2
            for (int j = 1; j <= s2.length(); ++j) {
                int i1 = da.get(s2.charAt(j - 1));
                int j1 = db;

                int cost = 1;
                if (s1.charAt(i - 1) == s2.charAt(j - 1)) {
                    cost = 0;
                    db = j;
                }

                distances[i + 1][j + 1] = min(
                        distances[i][j] + cost, // substitution
                        distances[i + 1][j] + 1, // insertion
                        distances[i][j + 1] + 1, // deletion
                        distances[i1][j1] + (i - i1 - 1) + 1 + (j - j1 - 1));
            }

            da.put(s1.charAt(i - 1), i);
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
        for (int j = 0; j < newTerm.length() - N; ++j) {
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
