package webdata;

import webdata.utils.LetterProbability;
import webdata.utils.QueryHistory;

import java.io.*;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.ArrayList;

/**
 *
 */
public class IndexWriter{

    static final String tokenDictFileName = "tokenDict";
    static final String productDictFileName = "productDict";
    static final String reviewDataFileName = "reviewData";
    static final String tokenNGIFileName = "tokenNGI";
    static final String productNGIFileName = "productNGI";
    static final String productPostingListFileName = "productPostingList";
    static final String tokenPostingListFileName = "tokenPostingList";
    static final String tokenLetterProbabilityFileName = "tokenLetterProbabilities";
    static final String productLetterProbabilityFileName = "productLetterProbabilities";
    static final String historyFileName = "queryHistory";

    private final String tokensFileName = "tokenFile";
    private final String productsFileName = "productFile";
    private final String sortedIndicator = "_sorted";

    /**
     * Given product review data, creates an on disk index
     * inputFile is the path to the file containing the review data
     * dir is the directory in which all index files will be created
     * if the directory does not exist, it should be created
     * @param inputFile The path to the file containing the review data.
     * @param dir the directory in which all index files will be created if the directory does not exist, it should be
     *            created.
     */
    public void write(String inputFile, String dir) {
        createDir(dir);
        String historyDir = dir + File.separator + "history";
        createDir(historyDir);

        String sortedTokensFilePath = dir + File.separator + tokensFileName + sortedIndicator;
        String sortedProductsFilePath = dir + File.separator + productsFileName + sortedIndicator;

        ReviewsParser parser = new ReviewsParser();
        parser.parseFile(inputFile);

        ReviewData rd = new ReviewData(parser.getProductIds(), parser.getReviewHelpfulnessNumerator(),
                parser.getReviewHelpfulnessDenominator(), parser.getReviewScore(),
                parser.getTokensPerReview(), parser.getNumOfReviews());

        try (ObjectOutputStream reviewDataWriter = new ObjectOutputStream(
                new FileOutputStream(dir + File.separator + reviewDataFileName))) {
            reviewDataWriter.writeObject(rd);
        } catch(IOException e) {
            System.err.println(e.getMessage());
            System.exit(1);
        }
        rd.clear();
        parser.clear();

        String tmpDirName = createTempDir(dir);
        Sorter sorter = new Sorter(new ArrayList<>(parser.getTokenSet()),
                                   new ArrayList<>(parser.getProductIdSet()),
                                   tmpDirName);
        sorter.sort(inputFile, sortedTokensFilePath, sortedProductsFilePath);
        removeIndex(tmpDirName);


        Dictionary tokenDict = buildDictionary(parser.getNumOfTokens(), sortedTokensFilePath,
                false, dir, sorter.getTokensArray());
        Dictionary productDict = buildDictionary(parser.getNumOfproducts(), sortedProductsFilePath,
                true, dir, sorter.getProductIdsArray());

        writeObject(dir, tokenDictFileName, tokenDict);
        writeObject(dir, productDictFileName, productDict);

        NGramIndex tokenNGI = new NGramIndex(sorter.getTokensArray());
        NGramIndex productNGI = new NGramIndex(sorter.getProductIdsArray());

        writeObject(dir, tokenNGIFileName, tokenNGI);
        writeObject(dir, productNGIFileName, productNGI);

        LetterProbability tokenLp = new LetterProbability(sorter.getTokensArray());
        LetterProbability productLp = new LetterProbability(sorter.getProductIdsArray());
        writeObject(dir, tokenLetterProbabilityFileName, tokenLp);
        writeObject(dir, productLetterProbabilityFileName, productLp);

        QueryHistory qh = new QueryHistory();
        writeObject(historyDir, historyFileName, qh);
    }

    /**
     * Write a Serializable object.
     * @param dir The directory to write to
     * @param fileName The file name to write to
     * @param o The object to write
     */
    public static void writeObject(String dir, String fileName, Object o) {
        try {
            ObjectOutputStream writer = new ObjectOutputStream(new FileOutputStream(dir + File.separator + fileName));
            writer.writeObject(o);
            writer.close();

        } catch(IOException e) {
            System.err.println(e.getMessage());
            System.exit(1);
        }
    }

    /**
     * Create a new directory. If already exists, delete it's contents.
     * @param dir The directory to create.
     */
    private void createDir(String dir) {
        File dirFile = new File(dir);
        if (dirFile.exists()) {
            purgeDirectory(dirFile);
        } else {  // Create it
            try{
                dirFile.mkdir();
            }
            catch(SecurityException e){
                System.err.println(e.getMessage());
                System.exit(1);
            }
        }
    }

    /**
     * Recursively delete all files and all subdirectories of the given directory
     * @param dir
     */
    private void purgeDirectory(File dir) {
        File[] entries = dir.listFiles();
        if (entries != null) {
            for(File file: entries){
                if (file.isDirectory())
                    purgeDirectory(file);
                file.delete();
            }
        }
    }

    /**
     * Create temporary directory for sorting
     * @param dir Directory  to create
     * @return The directory's name
     */
    private String createTempDir(String dir) {
        String tmpDirName = dir + File.separator + "tmp";
        removeIndex(tmpDirName);
        File tmpDir = new File(tmpDirName);
        if (!tmpDir.exists()) {
            try{
                tmpDir.mkdir();
            }
            catch(SecurityException e){
                System.err.println(e.getMessage());
                System.exit(1);
            }
        }
        return tmpDirName;
    }

    /**
     * Delete all index files by removing the given directory.
     * @param dir The directory to remove the index from.
     */
    public void removeIndex(String dir) {
        File dirFile = new File(dir);
        if (dirFile.exists()) {
            purgeDirectory(dirFile);
            dirFile.delete();
        }
    }

    /**
     * Build a dictionary object
     * @param numOfTerms Number of terms in the file
     * @param out The sorted file of terms
     * @param isProduct Indicates if the term is productId or token
     * @param dir The directory in which the dictionary is saved
     * @param mapping A map of a number to term (i is mapped to the string at index i)
     * @return The built dictionary
     */
    private Dictionary buildDictionary(int numOfTerms, String out, Boolean isProduct, String dir,
                                       ArrayList<String> mapping) {
        Dictionary dict = new Dictionary(numOfTerms, out, isProduct, dir, mapping);
        /* Delete sorted */
        try {
            Files.deleteIfExists(Paths.get(out));
        } catch (IOException e) {
            e.printStackTrace();
        }
        return dict;
    }
}
