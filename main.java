package webdata;

import java.text.SimpleDateFormat;
import java.util.*;

public class main {
    public static void main(String[] args) {
        String dir = "/Users/shahaf/Documents/UNI/אחזור מידע באינטרנט/ex3/indexFiles";
//        String file = "/Users/shahaf/Documents/UNI/אחזור מידע באינטרנט/datasets/Movies_&_TV.txt";
        String file = "/Users/shahaf/Documents/UNI/אחזור מידע באינטרנט/ex3/100.txt";
//        String file = "/Users/shahaf/Documents/UNI/אחזור מידע באינטרנט/datasets/Musical_Instruments.txt";
//        String file = "/Users/shahaf/Documents/UNI/אחזור מידע באינטרנט/datasets/Electronics.txt";


        IndexWriter siw = new IndexWriter();
        siw.write(file, dir);
        IndexReader ir = new IndexReader(dir);

//        test1(ir);
        test2(ir);

        siw.removeIndex(dir);

    }

    private static void takeTime(String msg) {
        String timeStamp =  new SimpleDateFormat("HH.mm.ss.SS").format(new java.util.Date());
        System.err.println(msg + " at " + timeStamp);
    }

    private static void test2(IndexReader ir) {
//        String[] queryTerms = {"best", "ever", "best", "query", "ever", "ever", "shahaf", "shahaf", "zohar", "zohar"};
//        String[] queryTerms = {"Good", "dog", "FoOd", "for", "sensitive", "dogs"};
        String[] queryTerms = {"love", "candy"};
        Vector<String> query = new Vector<>(Arrays.asList(queryTerms));
        ReviewSearch rs = new ReviewSearch(ir);
        Enumeration<Integer> res;

//        res = rs.vectorSpaceSearch(query.elements(), 5);
//        while (res.hasMoreElements()) {
//            System.out.println(res.nextElement());
//        }
//
//        res = rs.languageModelSearch(query.elements(), 0.25, 5);
//        while (res.hasMoreElements()) {
//            System.out.println(res.nextElement());
//        }

        Collection<String> goodProducts = rs.productSearch(query.elements(), 10);
        for (String product: goodProducts) {
            System.out.println(product);
        }

    }

    private static void test1(IndexReader ir) {
        System.err.println("\n\n--------- Started TEST ---------\n\n");
//        int[] vals = {-1, 0, 1, 1000, 99, 10001, 2000000, 8000000, 7850072, 7850073, 7850071};
//        for (int val: vals) {
//            System.out.println("Value: " + val);
//            System.out.println(ir.getProductId(val));
//            System.out.println(ir.getReviewScore(val));
//            System.out.println(ir.getReviewHelpfulnessNumerator(val));
//            System.out.println(ir.getReviewHelpfulnessDenominator(val));
//            System.out.println(ir.getReviewLength(val));
//        }

        /* ****/
        String[] svals = {"spill" , "wistful", "peaceful" , "cough" ,	"crooked" , "special", 	"used" ,
                "opposite" , "imperfect" , "white",	"level", "alert", "coat", "fancy", "sleepy", "risk",
                "bye", "feeling", "manage", "story", "fixed", "rude", "modern", "ticket", "beginner",
                "summer", "observation", "respect", "dynamic", "delicious", "quiet", "bridge", "strong",
                "addicted", "reflective", "unlock", "silent", "attraction", "way", "assorted", "spoon",
                "scratch", "baseball", "lying", "unwieldy", "mug", "tall", "serious", "lumber", "drum",
                "living", "unable", "health", "whip", "tempt", "testy", "window", "part", "return", "value",
                "spray", "adhesive", "grip", "vague", "roasted", "smelly", "fold", "purple", "cloudy",
                "prevent", "mysterious", "replace", "bucket", "moor", "certain", "pleasant", "introduce",
                "cool", "jolly", "dreary", "caring", "voice", "nauseating", "pancake", "glossy", "cushion",
                "true", "exciting", "desk", "present", "try", "fabulous", "tiger", "tiny", "cry",
                "inexpensive", "straw", "ugly", "circle", "holiday"};
        System.out.println();
        String tokenFrequencyStart =  new SimpleDateFormat("HH.mm.ss.SS").format(new java.util.Date());
        for (String s: svals) {
            int i = ir.getTokenFrequency(s);
        }
        String tokenFrequencyEnd =  new SimpleDateFormat("HH.mm.ss.SS").format(new java.util.Date());
        System.err.println("getTokenFrequency: start: " + tokenFrequencyStart + " , end: " + tokenFrequencyEnd);


        tokenFrequencyStart =  new SimpleDateFormat("HH.mm.ss.SS").format(new java.util.Date());
        for (String s: svals) {
            Enumeration<Integer> e = ir.getReviewsWithToken(s);
        }
        tokenFrequencyEnd =  new SimpleDateFormat("HH.mm.ss.SS").format(new java.util.Date());
        System.err.println("getReviewsWithToken: start: " + tokenFrequencyStart + " , end: " + tokenFrequencyEnd);

//        String[] svals2 = {"tall", "health", "apple"};
//        for (String s: svals2) {
//            System.out.println(ir.getTokenFrequency(s));
//            System.out.println(ir.getTokenCollectionFrequency(s));
//            Enumeration<Integer> e = ir.getReviewsWithToken(s);
//            while (e.hasMoreElements()) {
//                System.out.println(e.nextElement());
//            }
//        }
//        /* ****/
//
//        System.out.println();
//        System.out.println("Parameterless methods");
//        System.out.println(ir.getNumberOfReviews());
//        System.out.println(ir.getTokenSizeOfReviews());
//        String[] pvals = {"B0001Z3TLQ", "B0019CW0HE", "B0019CW0HF"};
//        for (String s: pvals) {
//            Enumeration<Integer> e = ir.getProductReviews(s);
//            while (e.hasMoreElements()) {
//                System.out.println(e.nextElement());
//            }
//        }
    }


}

//import org.omg.PortableInterceptor.SYSTEM_EXCEPTION;
//import webdata.IndexReader;
//import webdata.SlowIndexWriter;
//
//import java.io.*;
//import java.util.Enumeration;
//
//public class main {
//
//    static final String INDICES_DIR_NAME = "indices";
//    static final String REVIEWS_FILE_NAME_1000 = "1000.txt";
//
//    public static void main(String[] args) {
//        final String dir = System.getProperty("user.dir");
//        final String indicesDir = dir + File.separatorChar + INDICES_DIR_NAME;
//
//        buildIndex(indicesDir);
//        queryWordIndex(indicesDir);
//        queryMetaData(indicesDir);
//        queryReviewMetaData(indicesDir);
//        queryProductIndex(indicesDir);
//
//        deleteIndex(indicesDir);
//
//    }
//
//    private static void queryReviewMetaData(String indicesDir) {
//        IndexReader indexReader = new IndexReader(indicesDir);
//        int[] ridsTestCases = {1, 11, 12, 999, 1001, 0, -2, 10, 32, 522};
//        for (int rid : ridsTestCases) {
//            System.out.println("rid: " + rid + " " +
//                    indexReader.getProductId(rid) + " " +
//                    indexReader.getReviewScore(rid) + " " +
//                    indexReader.getReviewHelpfulnessNumerator(rid) + " " +
//                    indexReader.getReviewHelpfulnessDenominator(rid) + " " +
//                    indexReader.getReviewLength(rid));
//        }
//    }
//
//    private static void queryMetaData(String indicesDir) {
//        IndexReader indexReader = new IndexReader(indicesDir);
//        System.out.println(indexReader.getNumberOfReviews());
//        System.out.println(indexReader.getTokenSizeOfReviews());
//    }
//
//
//    private static void testGetReviewsWithToken(IndexReader indexReader,
//                                                String[] wordTestCases) {
//        System.out.println("Checking getReviewsWithToken...");
//        for (String word : wordTestCases) {
//            Enumeration<Integer> res = indexReader.getReviewsWithToken(word);
//            System.out.print(word + ": " + System.lineSeparator());
//            while (res.hasMoreElements()) {
//                System.out.print(res.nextElement().toString() + " ");
//            }
//            System.out.println();
//        }
//    }
//
//    private static void testGetTokenFrequency(IndexReader indexReader,
//                                              String[] wordTestCases) {
//        System.out.println("Checking getTokenFrequency...");
//        for (String word : wordTestCases) {
//            int numOfReviews = indexReader.getTokenFrequency(word);
//            System.out.println(word + ": " + "numOfReviews: " + numOfReviews);
//        }
//    }
//
//    private static void testGetTokenCollectionFrequency(IndexReader indexReader,
//                                                        String[] wordTestCases) {
//        System.out.println("Checking getTokenCollectionFrequency...");
//        for (String word : wordTestCases) {
//            int numOfMentions = indexReader.getTokenCollectionFrequency(word);
//            System.out.println(word + ": " + "numOfMentions: " + numOfMentions);
//        }
//    }
//
//
//    private static void queryWordIndex(String indicesDir) {
//        IndexReader indexReader = new IndexReader(indicesDir);
//        String[] wordTestCases = {"0", "bulba", "zzz", "1", "9oz", "a", "crunchy", "how", "laxative",
//                "prefer", "storebought", "zucchini", "the"};
//        testGetReviewsWithToken(indexReader, wordTestCases);
//        testGetTokenFrequency(indexReader, wordTestCases);
//        testGetTokenCollectionFrequency(indexReader, wordTestCases);
//
//    }
//
//
//    private static void queryProductIndex(String indicesDir) {
//        IndexReader indexReader = new IndexReader(indicesDir);
//        String[] productTestCases = {"A009ASDF5", "B099ASDF5", "B0001PB9FE", "B0002567IW", "B000ER6YO0",
//                "B000G6RYNE", "B006F2NYI2", "B009HINRX8", "B001E4KFG0"};
//        for (String pid : productTestCases) {
//            Enumeration<Integer> res = indexReader.getProductReviews(pid);
//            System.out.print(pid + ": " + System.lineSeparator());
//            while (res.hasMoreElements()) {
//                System.out.print(res.nextElement().toString() + " ");
//            }
//            System.out.println();
//        }
//    }
//
//
//    private static void deleteIndex(String indicesDir) {
//        SlowIndexWriter slowIndexWriter = new SlowIndexWriter();
//        slowIndexWriter.removeIndex(indicesDir);
//    }
//
//    private static void buildIndex(String indicesDir) {
//
//        SlowIndexWriter slowIndexWriter = new SlowIndexWriter();
//        slowIndexWriter.slowWrite(REVIEWS_FILE_NAME_1000, indicesDir);
//    }
//}
