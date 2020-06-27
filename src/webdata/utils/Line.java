package webdata.utils;

/**
 * This class represents a line in the relation files with idX and two additional columns separated by spaces.
 */
public class Line implements Comparable<Line>{
    private static final String COLUMN_DELIMINATOR = "#";
    private static final int TERM_INDEX = 0;
    private static final int REVIEW_ID_INDEX = 1;
    private static final int FREQUENCY_INDEX = 2;

    // Data members
    private int term;
    private int reviewId;
    private int frequency;

    /**
     * Initializes the line with the given line string.
     * @param line Line string.
     */
    public Line(String line){
        this.initLineParts(line);
    }

    /**
     * Initializes the different line parts.
     */
    private void initLineParts(String line){
        String[] lineParts = line.split(COLUMN_DELIMINATOR);
        this.term = Integer.parseInt(lineParts[TERM_INDEX]);
        this.reviewId = Integer.parseInt(lineParts[REVIEW_ID_INDEX]);
        this.frequency = Integer.parseInt(lineParts[FREQUENCY_INDEX]);
    }

    /**
     * Get the line's term
     */
    public int getTerm() {
        return term;
    }

    /**
     * Get the line's reviewId
     */
    public int getReviewId() {
        return reviewId;
    }

    /**
     * Get the line's reviewId
     */
    public int getFrequency() { return frequency; }

    /**
     * Set the line's frequency member
     * @param freq The frequency to set
     */
    public void setFrequency(int freq) { this.frequency = freq; }

    /**
     * Returns an int representing the order between this line and the given other line.
     * @param o The other line to compare to.
     * @return Negative number if this line is smaller than other, 0 if equal, and positive if it is greater.
     */
    @Override
    public int compareTo(Line o) {
        if (this.term == o.term){
            return this.reviewId - o.reviewId;
        }
        return this.term - o.term;
    }

    /**
     * Returns a string representing the line.
     * @return The line string representation.
     */
    @Override
    public String toString() {
        return this.term + "#" + this.reviewId  + "#" + this.frequency;
    }

    @Override
    public int hashCode() {
        int code = this.toString().hashCode();
        return code;
    }

    @Override
    public boolean equals(Object o) {

        // If the object is compared with itself then return true
        if (o == this) {
            return true;
        }

        /* Check if o is an instance of Line or not
          "null instanceof [type]" also returns false */
        if (!(o instanceof Line)) {
            return false;
        }

        // typecast o to Line so that we can compare data members
        Line l = (Line) o;

        // Compare the data members and return accordingly
        return compareTo(l) == 0;
    }
}
