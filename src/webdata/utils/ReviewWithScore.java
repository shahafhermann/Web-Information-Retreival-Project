package webdata.utils;

public class ReviewWithScore implements Comparable<ReviewWithScore> {

    private int reviewNumber;
    private double score;


    /**
     * Constructor
     */
    public ReviewWithScore(int reviewNumber, double score) {
        this.reviewNumber = reviewNumber;
        this.score = score;
    }

    public int getReviewNumber() {
        return reviewNumber;
    }

    public double getScore() {
        return score;
    }

    /**
     * Returns an int representing the order between this line and the given other line.
     * @param o The other line to compare to.
     * @return Negative number if this line is smaller than other, 0 if equal, and positive if it is greater.
     */
    @Override
    public int compareTo(ReviewWithScore o) {
        double res = o.score - this.score;
        if (res > 0) return 1;
        else if (res == 0) return this.reviewNumber - o.reviewNumber;
        else return -1;
    }
}
