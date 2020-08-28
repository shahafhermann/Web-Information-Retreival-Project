package webdata.utils;

/**
 * Wrapper class for an edit in Damerau-Levenshtein edit distance algorithm.
 */
public class Edit {

    private String type;
    private int index;

    /**
     * Constructor
     * @param type Error type
     * @param index Index in the word where the error happened
     */
    public Edit(String type, int index) {
        this.type = type;
        this.index = index;
    }

    /**
     * Get the error type
     */
    public String getType() { return type; }

    /**
     * Get the index in the word where the error happened
     */
    public int getIndex() { return index; }
}
