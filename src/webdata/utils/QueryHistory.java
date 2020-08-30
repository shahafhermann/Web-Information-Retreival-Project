package webdata.utils;

import java.io.Serializable;
import java.util.TreeMap;

/**
 * Singleton class, recording all queries entered by the user (as long as the index is not erased).
 * This class only records single terms entered in the queries,
 * rather than pairs or more (following the unigram language model).
 */
public class QueryHistory implements Serializable {

    private TreeMap<String, Integer> history = new TreeMap<>();

    /**
     * Private and empty constructor
     */
    public QueryHistory() {}

    /**
     * Save the given term to the history, incrementing it's count if already exists or adding as new if not.
     * If the term is null, do nothing.
     * @param term The term to save
     */
    public void save(String term) {
        if (term == null) return;

        int count = history.getOrDefault(term, 0);
        history.put(term, count + 1);
    }

    public TreeMap<String, Integer> getHistory() {
        return history;
    }
}
