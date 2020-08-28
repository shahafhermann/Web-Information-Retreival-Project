package webdata.utils;

import java.util.ArrayList;

/**
 * Wrapper for Damerau-Levenshtein minimum edit distance and the list of required edits to correct a word.
 */
public class DLDistance {

    private String wrong;
    private String correct;
    int[][] distanceMatrix;

    /**
     * Constructor
     * @param wrong The "misspelled" word
     * @param correct The "correct" word
     * @param distanceMatrix The distance matrix calculated with the Damerau-Levenshtein edit distance algorithm
     */
    public DLDistance(String wrong, String correct, int[][] distanceMatrix) {
        this.wrong = wrong;
        this.correct = correct;
        this.distanceMatrix = distanceMatrix;
    }

    /**
     * Get the misspelled word
     */
    public String getWrong() {
        return wrong;
    }

    /**
     * Get the correct word
     */
    public String getCorrect() {
        return correct;
    }

    /**
     * Get the minimum edit distance
     */
    public int getDistance() {
        return distanceMatrix[correct.length() + 1][wrong.length() + 1];
    }

    /**
     * Get a list of the edits required to turn 'wrong' into 'correct'.
     * The list is constructed using back tracing for dynamic programming.
     * @return The list of edits as an ArrayList of Edit.
     */
    public ArrayList<Edit> getEdits() {
        ArrayList<Edit> edits = new ArrayList<>();

        int i = correct.length();
        int j = wrong.length();
        while (i != 0 && j != 0) {
            if (i > 1 && j > 1 && correct.charAt(i - 1) == wrong.charAt(j - 2)
                    && correct.charAt(i - 2) == wrong.charAt(j - 1)) {
                if (distanceMatrix[i - 2][j - 2] < distanceMatrix[i][j]) {
                    edits.add(0, new Edit("trans", i - 2));
                    i -= 2;
                    j -= 2;
                    continue;
                }
            }

            int sub = distanceMatrix[i][j];
            int ins = distanceMatrix[i + 1][j];
            int del = distanceMatrix[i][j + 1];
            int min = Utils.min(sub, ins, del);

            if (min == sub) {
                if (distanceMatrix[i + 1][j + 1] > distanceMatrix[i][j]) {
                    edits.add(0, new Edit("sub", i - 1));
                }
                i--;
                j--;
            } else if (min == ins) {
                edits.add(0, new Edit("ins", i - 1));
                j--;
            } else if (min == del) {
                edits.add(0, new Edit("del", i - 1));
                i--;
            }
        }

        return edits;
    }
}
