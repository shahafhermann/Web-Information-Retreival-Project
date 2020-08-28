package webdata;

import webdata.utils.LetterProbability;
import webdata.utils.Utils;

import java.io.*;
import java.util.*;

public class IndexReader {

    Dictionary tokenDict;
    Dictionary productDict;
    ReviewData rd;
    NGramIndex tokenNGI;
    NGramIndex productNGI;
    LetterProbability tokenLp;
    LetterProbability productLp;

    /**
     * Creates an IndexReader which will read from the given directory
     * @param dir The directory to read from.
     */
    public IndexReader(String dir) {
        tokenDict = (Dictionary) readObject(dir, IndexWriter.tokenDictFileName);
        productDict = (Dictionary) readObject(dir, IndexWriter.productDictFileName);
        rd = (ReviewData) readObject(dir, IndexWriter.reviewDataFileName);
        tokenNGI = (NGramIndex) readObject(dir, IndexWriter.tokenNGIFileName);
        productNGI = (NGramIndex) readObject(dir, IndexWriter.productNGIFileName);
        tokenLp = (LetterProbability) readObject(dir, IndexWriter.tokenLetterProbabilityFileName);
        productLp = (LetterProbability) readObject(dir, IndexWriter.productLetterProbabilityFileName);
    }

    private Object readObject(String dir, String fileName) {
        try {
            ObjectInputStream reader = new ObjectInputStream(new FileInputStream(dir + File.separator + fileName));
            Object o = reader.readObject();
            reader.close();
            return o;

        } catch(IOException|ClassNotFoundException e) {
            System.err.println(e.getMessage());
            System.exit(1);
        }
        return null;  // We'll never reach this
    }

    // ------------------------------------------------------------ //

    /**
     * Get the the n-grams for the given term
     * @param term Term to find n-grams of
     * @param isProduct Indicate if this method should refer to the product or token index.
     * @return The n-grams as a String array
     */
    public String[] getNGrams(String term, boolean isProduct) {
        return (isProduct) ?
                Utils.findNGrams(productNGI.getN(), NGramIndex.EDGE_MARK, term) :
                Utils.findNGrams(tokenNGI.getN(), NGramIndex.EDGE_MARK, term);
    }

    /**
     * Find all term numbers that have ngrams in common with term.
     * @param ngrams String array of ngrams to find common terms with.
     * @param isProduct Indicate if this method should refer to the product or token index.
     * @return A set of the term numbers with common ngrams.
     */
    public Set<Integer> findTermsWithCommonNgrams(String[] ngrams, boolean isProduct) {
        Set<Integer> allTerms = new HashSet<>();

        // For each "gram", get the list of other terms that contain this "gram"
        for (String gram : ngrams) {
            ArrayList<Integer> commonTerms = (isProduct) ?
                                                            productNGI.getNGrams(gram) :
                                                            tokenNGI.getNGrams(gram);
            allTerms.addAll(commonTerms);
        }

        return allTerms;
    }

    /**
     * Find the term corresponding with the given id.
     * @param id The term's number
     * @param isProduct Indicate if this method should refer to the product or token dictionary.
     * @return The term corresponding with the id.
     */
    public String getTermById(int id, boolean isProduct) {
        return (isProduct) ? productDict.searchTerm(id) : tokenDict.searchTerm(id);
    }

    // ------------------------------------------------------------ //

    /**
     * @param reviewId The review to get the product id for.
     * @return The product identifier for the given review.
     *         Returns null if there is no review with the given identifier.
     */
    public String getProductId(int reviewId) {
        return ((1 <= reviewId) && (reviewId <= rd.getNumOfReviews())) ? rd.getReviewProductId(reviewId - 1) : null;
    }

    /**
     * @param reviewId The review to get the score for.
     * @return The score for a given review.
     *         Returns -1 if there is no review with the given identifier.
     */
    public int getReviewScore(int reviewId) {
        return ((1 <= reviewId) && (reviewId <= rd.getNumOfReviews())) ?
                rd.getScore(reviewId - 1) : -1;
    }

    /**
     * @param reviewId The review to get the numerator for the helpfulness for.
     * @return The numerator for the helpfulness of a given review
     *         Returns -1 if there is no review with the given identifier
     */
    public int getReviewHelpfulnessNumerator(int reviewId) {
        return ((1 <= reviewId) && (reviewId <= rd.getNumOfReviews())) ?
                rd.getHelpfulnessNumerator(reviewId - 1) : -1;
    }

    /**
     * @param reviewId The review to get the denominator for the helpfulness for.
     * @return The denominator for the helpfulness of a given review
     *         Returns -1 if there is no review with the given identifier
     */
    public int getReviewHelpfulnessDenominator(int reviewId) {
        return ((1 <= reviewId) && (reviewId <= rd.getNumOfReviews())) ?
                rd.getHelpfulnessDenominator(reviewId - 1) : -1;
    }

    /**
     * @param reviewId The review to get the number of tokens for.
     * @return The number of tokens in a given review
     *         Returns -1 if there is no review with the given identifier
     */
    public int getReviewLength(int reviewId) {
        return ((1 <= reviewId) && (reviewId <= rd.getNumOfReviews())) ?
                rd.getTokensPerReview(reviewId - 1) : -1;
    }


    // ------------------------------------------------------------ //


    /**
     * @param token The token to check.
     * @return The number of reviews containing a given token (i.e., word)
     *         Returns 0 if there are no reviews containing this token
     */
    public int getTokenFrequency(String token) {
        int i = tokenDict.searchTerm(token.toLowerCase());
        if (i < 0 || i >= tokenDict.getNumOfTerms()) {
            return 0;
        }
        long pos = tokenDict.getPostingPtr(i);
        return tokenDict.readLength(pos);
    }

    /**
     * @param token The token to check.
     * @return The number of times that a given token (i.e., word) appears in the reviews indexed
     *         Returns 0 if there are no reviews containing this token
     */
    public int getTokenCollectionFrequency(String token) {
        int i = tokenDict.searchTerm(token.toLowerCase());
        if (i < 0 || i >= tokenDict.getNumOfTerms()) {
            return 0;
        }
        return tokenDict.getFrequency(i);
    }

    /**
     * @param token The token to check.
     * @return A series of integers of the form id-1, freq-1, id-2, freq-2, ... such that
     *         id-n is the n-th review containing the given token and freq-n is the number of times that the token
     *         appears in review id-n.
     *         Note that the integers should be sorted by id.
     *         Returns an empty Enumeration if there are no reviews containing this token.
     */
     public Enumeration<Integer> getReviewsWithToken(String token) {
         return enumHelper(tokenDict, token.toLowerCase());
     }


     // --------------------------------------------------------- //


     /**
     * @return The number of product reviews available in the system.
     */
    public int getNumberOfReviews() {
        return rd.getNumOfReviews();
    }

    /**
     * @return The number of tokens in the system (Tokens should be counted as many times as they appear).
     */
    public int getTokenSizeOfReviews() {
        int tokenCount = 0;
        for (int i = 1; i <= getNumberOfReviews(); ++i) {
            tokenCount += getReviewLength(i);
        }
        return tokenCount;
    }


    // ---------------------------------------------------------- //


    /**
     * @param productId The id of the product to check.
     * @return The ids of the reviews for a given product identifier.
     *         Note that the integers returned should be sorted by id.
     *         Returns an empty Enumeration if there are no reviews for this product.
     */
    public Enumeration<Integer> getProductReviews(String productId) {
        return enumHelper(productDict, productId);
    }


    // ---------------------------------------------------------- //


    /**
     * Get the Enumaration list for the given Dictionary and term.
     * @param dict Dictionary
     * @param term Term
     * @return Enumaration
     */
    private Enumeration<Integer> enumHelper(Dictionary dict, String term) {
        int i = dict.searchTerm(term);
        if (i < 0 || i >= tokenDict.getNumOfTerms()) {
            return new Vector<Integer>().elements();
        }
        long pos = dict.getPostingPtr(i);
        long nextPos = (i + 1 < dict.getNumOfTerms()) ? dict.getPostingPtr(i + 1) : -1;
        Integer[] list = dict.read(pos, nextPos);

        Vector<Integer> reviewsWithToken = new Vector<>(Arrays.asList(list));
        return reviewsWithToken.elements();
    }
}
