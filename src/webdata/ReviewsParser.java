package webdata;

import java.io.*;
import java.util.ArrayList;
import java.util.HashSet;

/**
 * A parser for a file of reviews.
 */
public class ReviewsParser {

    /* Data */
    private HashSet<String> tokenSet = new HashSet<>();
    private HashSet<String> productIdSet = new HashSet<>();

    private ArrayList<Byte> reviewScore = new ArrayList<>();
    private ArrayList<Short> reviewHelpfulnessNumerator = new ArrayList<>();
    private ArrayList<Short> reviewHelpfulnessDenominator = new ArrayList<>();
    private ArrayList<Short> tokensPerReview = new ArrayList<>();
    private int numOfReviews = 0;
    private StringBuilder productIds = new StringBuilder();

    /* String constants */
    private final String SPLIT_TOKENS_REGEX = "[^A-Za-z0-9]+";

    /**
     * Empty the data structures stored in this instance.
     */
    void clear() {
        this.reviewScore = new ArrayList<>();
        this.reviewHelpfulnessNumerator = new ArrayList<>();
        this.reviewHelpfulnessDenominator = new ArrayList<>();
        this.tokensPerReview = new ArrayList<>();
        numOfReviews = 0;
        productIds = new StringBuilder();
    }

    /**
     * Return the tokens HashSet
     */
    HashSet<String> getTokenSet() { return tokenSet; }

    /**
     * Return the productIds HashSet
     */
    HashSet<String> getProductIdSet() { return productIdSet; }

    /**
     * Return the review scores as an ArrayList of Strings
     */
    ArrayList<Byte> getReviewScore() {
        return reviewScore;
    }

    /**
     * Return the review helpfulness numerator as an ArrayList of Strings
     */
    ArrayList<Short> getReviewHelpfulnessNumerator() {
        return reviewHelpfulnessNumerator;
    }

    /**
     * Return the review helpfulness denominator as an ArrayList of Strings
     */
    ArrayList<Short> getReviewHelpfulnessDenominator() {
        return reviewHelpfulnessDenominator;
    }

    /**
     * Return the review product IDs as an ArrayList of Strings
     */
    String getProductIds() {
        return productIds.toString();
    }

    /**
     * Return the number of token per review as an ArrayList of Strings
     */
    ArrayList<Short> getTokensPerReview() {
        return tokensPerReview;
    }

    /**
     * Return the number of reviews
     */
    int getNumOfReviews() { return numOfReviews;}

    /**
     * Return the number unique tokens
     */
    int getNumOfTokens() { return tokenSet.size(); }

    /**
     * Return the number of unique products
     */
    int getNumOfproducts() { return productIdSet.size(); }

    /**
     * Break a text to all it's tokens (alphanumeric).
     * @param text The text to break
     */
    private void breakText(String text) {
        String[] tokens = text.split(SPLIT_TOKENS_REGEX);
        int tokenCounter = 0;
        for (String token: tokens) {
            if (!token.isEmpty()) {
                tokenSet.add(token);

                ++tokenCounter;
            }
        }
        tokensPerReview.add((short) tokenCounter);
    }

    /**
     * Parse a string resembling a review helpfulness to it's numerator and denominator.
     * @param term The review helpfulness as String
     */
    private void writeReviewHelpfulness(String term) {
        String[] split = term.split("/");
        reviewHelpfulnessNumerator.add(Short.parseShort(split[0]));
        reviewHelpfulnessDenominator.add(Short.parseShort(split[1]));
    }

    /**
     * Parse a string resembling a review score.
     * @param term The review helpfulness as String
     */
    private void writeReviewScore(String term) {
        reviewScore.add(Byte.parseByte(term.split("\\.")[0]));
    }

    /**
     * Parse the file
     * @param inputFile The file to parse
     */
    void parseFile(String inputFile) {
        try (BufferedReader reader = new BufferedReader(new FileReader(new File(inputFile)), (int)Math.pow(2, 20))){
            String line = reader.readLine();
            String textBuffer = "";
            boolean textFlag = false;
            String term;
            while (line != null){

                if (textFlag && !line.startsWith("product/productId: ")) {
                    textBuffer = textBuffer.concat(" ").concat(line);
                    line = reader.readLine();
                    continue;
                }

                if (line.startsWith("product/productId: ")) {
                    textFlag = false;
                    if (!textBuffer.isEmpty()) {
                        breakText(textBuffer.toLowerCase());
                    }
                    ++numOfReviews;
                    term = line.substring(19);
                    productIds.append(term);
                    productIdSet.add(term);
                    line = reader.readLine();
                    continue;
                }

                if (line.startsWith("review/helpfulness: ")) {
                    writeReviewHelpfulness(line.substring(20));
                    line = reader.readLine();
                    continue;
                }

                if (line.startsWith("review/score: ")) {
                    writeReviewScore(line.substring(14));
                    line = reader.readLine();
                    continue;
                }

                if (line.startsWith("review/text:")) {
                    textFlag = true;
                    textBuffer = line.substring(12);
                    line = reader.readLine();
                    continue;
                }

                line = reader.readLine();
            }

            if (!textBuffer.isEmpty()) {
                breakText(textBuffer.toLowerCase());
            }
        } catch (IOException e) {
            System.err.println(e.getMessage());
            System.exit(1);
        }
    }
}
