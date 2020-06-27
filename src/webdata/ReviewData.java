package webdata;

import webdata.utils.Utils;
import java.io.Serializable;
import java.util.ArrayList;

/**
 * An object representing the data for the reviews
 */
class ReviewData implements Serializable {
    /* A String with all product IDs concatenated */
    private String productId;

    /* The length of a single product ID */
    private byte productIdLen;

    /* Array holding the numerator part of the review's helpfulness */
    private short[] reviewHelpfulnessNumerator;

    /* Array holding the denominator part of the review's helpfulness */
    private short[] reviewHelpfulnessDenominator;

    /* Array holding the review scores */
    private byte[] reviewScore;

    /* Array holding the number of tokens per review */
    private short[] tokensPerReview;

    /* The total number of reviews */
    private int numOfReviews;

    /**
     * Construct the review data object
     * @param productId String representing all product IDs concatenated
     * @param reviewHelpfulnessN ArrayList of Shorts representing all helpfulness Numerator data
     * @param reviewHelpfulnessD ArrayList of Shorts representing all helpfulness Denominator data
     * @param reviewScore ArrayList of Strings representing all review scores
     * @param tokensPerReview ArrayList of Strings representing the number of tokens per review
     * @param numOfReviews The total number of reviews
     */
    ReviewData (String productId, ArrayList<Short> reviewHelpfulnessN, ArrayList<Short> reviewHelpfulnessD,
                ArrayList<Byte> reviewScore, ArrayList<Short> tokensPerReview, int numOfReviews) {
        productIdLen = (byte) (productId.length() / numOfReviews);
        this.productId = productId;

        this.reviewHelpfulnessNumerator = new short[numOfReviews];
        this.reviewHelpfulnessDenominator = new short[numOfReviews];
        Utils.toPrimitiveArray(reviewHelpfulnessN, this.reviewHelpfulnessNumerator);
        Utils.toPrimitiveArray(reviewHelpfulnessD, this.reviewHelpfulnessDenominator);

        this.reviewScore = new byte[numOfReviews];
        Utils.toPrimitiveArray(reviewScore, this.reviewScore);

        this.tokensPerReview = new short[numOfReviews];
        Utils.toPrimitiveArray(tokensPerReview, this.tokensPerReview);
        this.numOfReviews = numOfReviews;
    }

    /**
     * Empty the data structures stored in this instance.
     */
    void clear() {
        this.productId = "";
        this.productIdLen = 0;
        this.reviewHelpfulnessNumerator = null;
        this.reviewHelpfulnessDenominator = null;
        this.reviewScore = null;
        this.tokensPerReview = null;
        this.numOfReviews = 0;
    }

    /**
     * Return the score for the requested review i
     */
    byte getScore(int i) { return reviewScore[i]; }

    /**
     * Return the helpfulness numerator for the requested review i
     */
    short getHelpfulnessNumerator(int i) { return reviewHelpfulnessNumerator[i]; }

    /**
     * Return the helpfulness denominator for the requested review i
     */
    short getHelpfulnessDenominator(int i) { return reviewHelpfulnessDenominator[i]; }

    /**
     * Return review i's productID
     */
    String getReviewProductId(int i) {
        return productId.substring(i * productIdLen, (i * productIdLen) + productIdLen);
    }

    /**
     * Return the number of tokens in review i
     */
    short getTokensPerReview(int i) { return tokensPerReview[i]; }

    /**
     * Return the number of reviews
     */
    int getNumOfReviews() { return numOfReviews; }
}
