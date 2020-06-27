package webdata.utils;

public class ProductWithScore implements Comparable<ProductWithScore> {

    private String productId;
    private double score;

    /**
     * Constructor
     */
    public ProductWithScore(String productId, double score) {
        this.productId = productId;
        this.score = score;
    }

    public String getProductId() {
        return productId;
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
    public int compareTo(ProductWithScore o) {
        double res = o.score - this.score;
        if (res > 0) return 1;
        else if (res == 0) return this.productId.compareTo(o.productId);
        else return -1;
    }
}
