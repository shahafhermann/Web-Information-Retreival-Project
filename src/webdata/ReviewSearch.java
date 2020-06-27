package webdata;

import com.sun.source.tree.Tree;
import webdata.utils.ProductWithScore;
import webdata.utils.ReviewWithScore;
import webdata.utils.Utils;

import java.util.*;

public class ReviewSearch {

    private static final int C = 30;
    private IndexReader ir;

    private final int deleteCost = 1;
    private final int insertCost = 1;
    private final int replaceCost = 1;

    /**
     * Constructor
     */
    public ReviewSearch(IndexReader iReader) {
        ir = iReader;
    }

    /* -------------------------------- Vector Space Search ------------------------------------ */

    private double[] calcLtf(double[] termFrequencies) {
        Arrays.parallelSetAll(termFrequencies,
                i -> ((termFrequencies[i] == 0) ? 0 : (1 + Math.log10(termFrequencies[i]))));
        return termFrequencies;
    }

    private double[] calcTidf(double[] termFrequencies) {
        int numOfReviews = ir.getNumberOfReviews();
        Arrays.parallelSetAll(termFrequencies,
                i -> ((termFrequencies[i] == 0) ? 0 : (Math.log10(numOfReviews / termFrequencies[i]))));
        return termFrequencies;
    }

    private TreeMap<String, Integer> histogramQuery(Enumeration<String> query) {
        List<String> queryList = Collections.list(query);
        TreeMap<String, Integer> hist = new TreeMap<>();
        for (String term: queryList) {
            Integer freq = (hist.keySet().contains(term.toLowerCase())) ? hist.get(term.toLowerCase()) + 1 : 1;
            hist.put(term.toLowerCase(), freq);
        }
        return hist;
    }

    private TreeMap<Integer, TreeMap<String, Integer>> getReviews(Set<String> queryTerms) {
        TreeMap<Integer, TreeMap<String, Integer>> reviews = new TreeMap<>();
        for (String term: queryTerms) {
            Enumeration<Integer> termReviewAndFrequencies = ir.getReviewsWithToken(term);
            while (termReviewAndFrequencies.hasMoreElements()) {
                int reviewNumber = termReviewAndFrequencies.nextElement();
                int frequencyInReview = termReviewAndFrequencies.nextElement();
                TreeMap<String, Integer> curMap;
                if (reviews.keySet().contains(reviewNumber)) { // If the current review number already exists in the reviews map
                    curMap = reviews.get(reviewNumber); // Get the currently existing "vector" of this review
                } else {
                    curMap = new TreeMap<>(); // Create a "vector" of all available unique query terms
                    for (String t: queryTerms) {
                        curMap.put(t, 0);
                    }
                }
                curMap.replace(term, frequencyInReview);
                reviews.put(reviewNumber, curMap);
            }
        }
        return reviews;
    }

    private double[] computeLTCOfQuery(TreeMap<String, Integer> queryHist) {
        double[] ltf = calcLtf(webdata.utils.Utils.integerCollectionToDoubleArray(queryHist.values()));
        double[] termFrequencies = new double[ltf.length];
        int i = 0;
        for (String token: queryHist.keySet()) {
            termFrequencies[i] = ir.getTokenFrequency(token);
            ++i;
        }
        termFrequencies = calcTidf(termFrequencies);

        double[] ltc = new double[ltf.length];
        double cosNormalization = 0;
        for (i = 0; i < ltf.length; ++i) {
            ltc[i] = ltf[i] * termFrequencies[i];
            cosNormalization += Math.pow(ltc[i], 2);
        }
        cosNormalization = Math.sqrt(cosNormalization);
        final double cosNormalizationConst = cosNormalization;
        if (cosNormalizationConst != 0) {
            Arrays.parallelSetAll(ltc,
                    j -> (ltc[j] / cosNormalizationConst));
        }
        return ltc;
    }

    private double calcReviewScore(double[] freqOfTerms, double[] qqq) {
        freqOfTerms = calcLtf(freqOfTerms);
        double score = 0;
        for (int i = 0; i < qqq.length; ++i) {
            score += freqOfTerms[i] * qqq[i];
        }
        return score;
    }

    /**
     * Returns a list of the id-s of the k most highly ranked reviews for the
     * given query, using the vector space ranking function lnn.ltc (using the
     * SMART notation)
     * The list should be sorted by the ranking
     */
    public Enumeration<Integer> vectorSpaceSearch(Enumeration<String> query, int k) {
        // Compute qqq:
        TreeMap<String, Integer> queryHist = histogramQuery(query);
        double[] ltc = computeLTCOfQuery(queryHist);
        TreeMap<String, Double> queryVec = new TreeMap<>();
        int i = 0;
        for (String term: queryHist.keySet()) {
            if (ltc[i] != 0) {
                queryVec.put(term, ltc[i]);
            }
            ++i;
        }

        // Compute ddd:
        TreeMap<Integer, TreeMap<String, Integer>> reviews = getReviews(queryVec.keySet());
        ReviewWithScore[] reviewWithScores = new ReviewWithScore[reviews.size()];
        i = 0;
        for (Integer review: reviews.keySet()) {
            double score = calcReviewScore(Utils.integerCollectionToDoubleArray(reviews.get(review).values()),
                                           Utils.doubleCollectionToDoubleArray(queryVec.values()));
            ReviewWithScore rws = new ReviewWithScore(review, score);
            reviewWithScores[i] = rws;
            ++i;
        }
        return getBestReviews(k, reviewWithScores);
    }


    /* -------------------------------- Language Model Search ------------------------------------ */


    private double calcMcProb(String token) {
        return ((double) ir.getTokenCollectionFrequency(token)) / ((double) ir.getTokenSizeOfReviews());
    }

    private double calcMdProb(String token, int reviewId) {
        Enumeration<Integer> reviewsWithFrequency = ir.getReviewsWithToken(token);
        int curReview;
        double curFrequency;
        while (reviewsWithFrequency.hasMoreElements()) {
            curReview = reviewsWithFrequency.nextElement();
            curFrequency = reviewsWithFrequency.nextElement();
            if (curReview == reviewId) {
                return curFrequency / (double) ir.getReviewLength(reviewId);
            }
        }
        return 0;
    }

    private double mixtureModelPerReview(List<String> query, double lambda, int reviewId) {
        double score = 1;
        double probabilty;
        for (String term: query) {
            probabilty = (lambda * calcMdProb(term, reviewId)) + ((1 - lambda) * calcMcProb(term));
            score *= probabilty;
        }
        return score;
    }

    private Set<Integer> getRelevantReviews(Set<String> querySet) {
        Set<Integer> relevantReviews = new HashSet<>();
        int curReview;
        for (String term: querySet) {
            Enumeration<Integer> termPostingList = ir.getReviewsWithToken(term);
            while (termPostingList.hasMoreElements()) {
                curReview = termPostingList.nextElement();
                termPostingList.nextElement();
                relevantReviews.add(curReview);
            }
        }
        return relevantReviews;
    }

    /**
     * Returns a list of the id-s of the k most highly ranked reviews for the
     * given query, using the language model ranking function, smoothed using a
     * mixture model with the given value of lambda
     * The list should be sorted by the ranking
     */
    public Enumeration<Integer> languageModelSearch(Enumeration<String> query, double lambda, int k) {
        ArrayList<String> queryList = new ArrayList<>();
        while (query.hasMoreElements()) {
            queryList.add(query.nextElement().toLowerCase());
        }
        Set<String> querySet = new HashSet<>(queryList);
        Set<Integer> relevantReviews = getRelevantReviews(querySet);

        ReviewWithScore[] reviewWithScores = new ReviewWithScore[relevantReviews.size()];
        int i = 0;
        double score;
        for (int reviewId: relevantReviews) {
            score = mixtureModelPerReview(queryList, lambda, reviewId);
            ReviewWithScore rws = new ReviewWithScore(reviewId, score);
            reviewWithScores[i] = rws;
            ++i;
        }
        return getBestReviews(k, reviewWithScores);

    }

    private Enumeration<Integer> getBestReviews(int k, ReviewWithScore[] reviewWithScores) {
        Arrays.sort(reviewWithScores);
        int numOfBestResults = Math.min(k, reviewWithScores.length);
        Integer[] bestResults = new Integer[numOfBestResults];
        for (int i = 0; i < bestResults.length; ++i) {
            bestResults[i] = reviewWithScores[i].getReviewNumber();
        }

        Vector<Integer> bestReviews = new Vector<>(Arrays.asList(bestResults));
        return bestReviews.elements();
    }


    /* ------------------------------------ Product Search ---------------------------------------- */


    /**
     * Returns a list of the id-s of the k most highly ranked productIds for the
     * given query using a function of your choice
     * The list should be sorted by the ranking
     *
     * Algorithm:
     * 1.	Find the top C reviews that match the query using VectorSpaceSearch (‘c’ is an integer greater than 1).
     * 2.	Assign each review with a weight corresponding it's ranking, where sum(weights)=1,
     *      and the higher the ranking – the bigger the weight.
     * 3.	Extract Product IDs of each of the above reviews.
     *      If productId is mentioned in more than one review, sum their weights.
     * 4.	For each productId get the posting list and find ALL the reviews that it appears in.
     * 5.	For each review of each productId, get the helpfulness and score.
     *      If helpfulness is not in range [0,1] then discard it.
     *      Then save (helpfulness*score) as the ‘new_review_score’.
     * 6.	Calculate the average of: average(‘all_new_review_scores_of_pid’)
     *      and median(‘all_new_review_scores_of_pid’) for each product.
     * 7.	Normalize the above result for each product by the sum of all results.
     * 8.	Combine the weights calculated in steps 3 and 7 and normalize.
     * 9.	Return top k.
     */
    public Collection<String> productSearch(Enumeration<String> query, int k) {
        // Find all relevant reviews according to the query
        Enumeration<Integer> allRelevantReviews = vectorSpaceSearch(query, C);

        List<Integer> allReviewsList = Collections.list(allRelevantReviews); // Convert to list

        // Assign each review with a weight corresponding it's position (sum(weights)=1)
        int numOfRelevantReviews = allReviewsList.size();
        double[] weights = new double[numOfRelevantReviews];
        int denominator = (numOfRelevantReviews * (numOfRelevantReviews + 1)) / 2;
        for (int i = 0; i < numOfRelevantReviews; ++i) {
            weights[i] = ((double) (numOfRelevantReviews - i)) / denominator;
        }

        // Extract Product IDs of each review.
        HashMap<String, Double> productWeightMap = new HashMap<>();
        for (int i = 0; i < numOfRelevantReviews; ++i) {
            String productId = ir.getProductId(allReviewsList.get(i));
            if (productWeightMap.containsKey(productId)) {
                double weight = productWeightMap.get(productId);
                productWeightMap.put(productId, weight + weights[i]);
            } else {
                productWeightMap.put(productId, weights[i]);
            }
        }

        // Iterate over all reviews of all products
        HashMap<String, Double> productNewScores = new HashMap<>();
        double sumOfScores = 0;
        for (String productId: productWeightMap.keySet()) {
            Enumeration<Integer> productReviews = ir.getProductReviews(productId);

            ArrayList<Double> newReviewScores = new ArrayList<>();
            while (productReviews.hasMoreElements()) {
                int reviewId = productReviews.nextElement();
                int score = ir.getReviewScore(reviewId);
                int helpfulnessNumerator = ir.getReviewHelpfulnessNumerator(reviewId);
                int helpfulnessDenominator = ir.getReviewHelpfulnessDenominator(reviewId);
                if (helpfulnessNumerator > helpfulnessDenominator) {
                    continue;
                }

                double helpfulness = (helpfulnessDenominator == 0) ?
                        1 : ((double) helpfulnessNumerator) / helpfulnessDenominator;
                newReviewScores.add(score * helpfulness);
            }

            // Get median and average and calculate a score for this product
            Collections.sort(newReviewScores);
            int numOfReviews = newReviewScores.size();
            double median;
            if (numOfReviews % 2 == 0) {
                median = (newReviewScores.get(numOfReviews / 2) + newReviewScores.get((numOfReviews / 2) - 1)) / 2;
            } else {
                median = newReviewScores.get(numOfReviews / 2);
            }
            double avg = 0;
            for (double newScore: newReviewScores) {
                avg += newScore;
            }
            avg /= numOfReviews;

            double newScore = (avg + median) / 2;
            productNewScores.put(productId, newScore);
            sumOfScores += newScore;
        }

        // Calculate final weight and normalize
        ProductWithScore[] productWithScores = new ProductWithScore[productWeightMap.size()];
        int i = 0;
        for (String productId: productWeightMap.keySet()) {
            double normalizedScore = productNewScores.get(productId) / sumOfScores;
            double curWeight = productWeightMap.get(productId);
            productWithScores[i] = new ProductWithScore(productId, (curWeight + normalizedScore) / 2);
            ++i;
        }

        // Sort by values and return top k
        Arrays.sort(productWithScores);
        int numOfBestResults = Math.min(k, productWithScores.length);
        ArrayList<String> bestResults = new ArrayList<>();
        for (i = 0; i < numOfBestResults; ++i) {
            bestResults.add(productWithScores[i].getProductId());
        }

        return bestResults;
    }


    ////////////////////////////////////////////////////////////////////////////////////////////////////////
    // **                        Damerau-Levenshtein Edit Distance Algorithm                           ** //
    ////////////////////////////////////////////////////////////////////////////////////////////////////////

    private int threeWayMin(int a, int b, int c) {
        return Math.min(Math.min(a, b), c);
    }

    /**
     * Compute the Damerau-Levenshtein distance between the specified source
     * string and the specified target string.
     */
    public int DLD(String source, String target) {
        int sourceLength = source.length();
        int targetLength = target.length();
        if (sourceLength == 0) return targetLength;
        if (targetLength == 0) return sourceLength;

        // Declare and init the distances matrix
        int[][] distances = new int[targetLength + 1][sourceLength + 1];
        for (int i = 0; i <= targetLength; ++i) {
            distances[i][0] = i;
        }
        for (int j = 0; j <= sourceLength; ++j) {
            distances[0][j] = j;
        }

        // Calculate weights
        int cost;
        for (int i = 1; i <= targetLength; ++i) {
            for (int j = 1; j <= sourceLength; ++j) {
                // todo: instead of 1, use weights[i][j]
                cost = (source.charAt(j - 1) == target.charAt(i - 1)) ? 0 : 1;
                distances[i][j] = threeWayMin((distances[i - 1][j] + deleteCost),
                        (distances[i - 1][j] + insertCost),
                        (distances[i - 1][j - 1] + cost));
            }
        }

        return distances[targetLength][sourceLength];








//        int deleteCost = 1, insertCost = 1, replaceCost = 1, swapCost = 1;
//
//        // Handle the cases where one of the input strings is empty
//        if (source.length() == 0) {
//            return target.length() * insertCost;
//        }
//        if (target.length() == 0) {
//            return source.length() * deleteCost;
//        }
//
//        // Init the distance matrix
//        int[][] distances = new int[source.length()][target.length()];
//        Map<Character, Integer> sourceIndexByCharacter = new HashMap<>();
//        if (source.charAt(0) != target.charAt(0)) {
//            distances[0][0] = Math.min(replaceCost, deleteCost + insertCost);
//        }
//        sourceIndexByCharacter.put(source.charAt(0), 0);
//
//        for (int i = 1; i < source.length(); i++) {
//            int deleteDistance = distances[i - 1][0] + deleteCost;
//            int insertDistance = (i + 1) * deleteCost + insertCost;
//            int matchDistance = i * deleteCost
//                    + (source.charAt(i) == target.charAt(0) ? 0 : replaceCost);
//            distances[i][0] = Math.min(Math.min(deleteDistance, insertDistance),
//                    matchDistance);
//        }
//        for (int j = 1; j < target.length(); j++) {
//            int deleteDistance = (j + 1) * insertCost + deleteCost;
//            int insertDistance = distances[0][j - 1] + insertCost;
//            int matchDistance = j * insertCost
//                    + (source.charAt(0) == target.charAt(j) ? 0 : replaceCost);
//            distances[0][j] = Math.min(Math.min(deleteDistance, insertDistance),
//                    matchDistance);
//        }
//        for (int i = 1; i < source.length(); i++) {
//            int maxSourceLetterMatchIndex = source.charAt(i) == target.charAt(0) ? 0
//                    : -1;
//            for (int j = 1; j < target.length(); j++) {
//                Integer candidateSwapIndex = sourceIndexByCharacter.get(target
//                        .charAt(j));
//                int jSwap = maxSourceLetterMatchIndex;
//                int deleteDistance = distances[i - 1][j] + deleteCost;
//                int insertDistance = distances[i][j - 1] + insertCost;
//                int matchDistance = distances[i - 1][j - 1];
//                if (source.charAt(i) != target.charAt(j)) {
//                    matchDistance += replaceCost;
//                } else {
//                    maxSourceLetterMatchIndex = j;
//                }
//                int swapDistance;
//                if (candidateSwapIndex != null && jSwap != -1) {
//                    int iSwap = candidateSwapIndex;
//                    int preSwapCost;
//                    if (iSwap == 0 && jSwap == 0) {
//                        preSwapCost = 0;
//                    } else {
//                        preSwapCost = distances[Math.max(0, iSwap - 1)][Math.max(0, jSwap - 1)];
//                    }
//                    swapDistance = preSwapCost + (i - iSwap - 1) * deleteCost
//                            + (j - jSwap - 1) * insertCost + swapCost;
//                } else {
//                    swapDistance = Integer.MAX_VALUE;
//                }
//                distances[i][j] = Math.min(Math.min(Math
//                        .min(deleteDistance, insertDistance), matchDistance), swapDistance);
//            }
//            sourceIndexByCharacter.put(source.charAt(i), i);
//        }
//        return distances[source.length() - 1][target.length() - 1];
    }
}
