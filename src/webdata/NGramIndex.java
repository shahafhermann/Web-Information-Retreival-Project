package webdata;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.TreeMap;

/**
 * A class for building an n-gram index of the given terms.
 */
public class NGramIndex implements Serializable {

    private int N = 2;  // Default to Bi-Gram
    public static String EDGE_MARK = "$";

    private TreeMap<String, ArrayList<Integer>> ngrams;

    /**
     * Constructor.
     * @param terms The terms to index
     */
    public NGramIndex(ArrayList<String> terms) {
        ngrams = new TreeMap<>();
        buildIndex(terms);
    }

    /**
     * Constructor, allows to change the default N. This is a set-up for future improvements.
     * @param terms The terms to index
     * @param n -gram.
     */
    public NGramIndex(ArrayList<String> terms, int n) {
        this.N = n;
        ngrams = new TreeMap<>();
        buildIndex(terms);
    }

    /**
     * Build the index
     * @param terms The terms to index
     */
    private void buildIndex(ArrayList<String> terms) {
        for (int i = 0; i < terms.size(); ++i) {
            String newTerm = EDGE_MARK.concat(terms.get(i)).concat(EDGE_MARK);
            for (int j = 0; j < newTerm.length() - N; ++j) {
                String gram = newTerm.substring(j, j + N);
                ArrayList<Integer> curList = (ngrams.containsKey(gram)) ? ngrams.get(gram) : new ArrayList<>();
                curList.add(i);
                ngrams.put(gram, curList);
            }
        }
    }

    /**
     * Get N
     */
    public int getN() { return N; }

    /**
     * Get the terms with the given n-gram
     * @param gram The n-gram to get
     * @return The terms as an ArrayList of Integers, representing the term's number.
     */
    public ArrayList<Integer> getNGrams(String gram) { return ngrams.get(gram); }
}
