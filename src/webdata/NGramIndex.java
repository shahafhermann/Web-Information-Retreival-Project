package webdata;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.TreeMap;

public class NGramIndex implements Serializable {

    private int N = 2;  // Default to Bi-Gram
    public static String EDGE_MARK = "$";

    private TreeMap<String, ArrayList<Integer>> ngrams;

    public NGramIndex(ArrayList<String> terms) {
        ngrams = new TreeMap<>();
        buildIndex(terms);
    }

    public NGramIndex(ArrayList<String> terms, int n) {
        this.N = n;
        ngrams = new TreeMap<>();
        buildIndex(terms);
    }

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

    public int getN() { return N; }

    public ArrayList<Integer> getNGrams(String term) { return ngrams.get(term); }
}
