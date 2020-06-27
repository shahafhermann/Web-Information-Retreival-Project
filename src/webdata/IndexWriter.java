package webdata;

import java.io.*;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.text.SimpleDateFormat;
import java.util.ArrayList;

/**
 *
 */
public class IndexWriter{

    static final String tokenDictFileName = "tokenDict";
    static final String productDictFileName = "productDict";
    static final String reviewDataFileName = "reviewData";
    static final String productPostingListFileName = "productPostingList";
    static final String tokenPostingListFileName = "tokenPostingList";
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
        File dirFile = new File(dir);
        if (dirFile.exists()) {
            removeFiles(dir);
        } else {  // Create it
            try{
                dirFile.mkdir();
            }
            catch(SecurityException e){
                System.err.println(e.getMessage());
                System.exit(1);
            }
        }

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

        try {
            /* Write the new files */
            ObjectOutputStream tokenDictWriter = new ObjectOutputStream(new FileOutputStream(dir + File.separator + tokenDictFileName));
            tokenDictWriter.writeObject(tokenDict);
            tokenDictWriter.close();

            ObjectOutputStream productDictWriter = new ObjectOutputStream(new FileOutputStream(dir + File.separator + productDictFileName));
            productDictWriter.writeObject(productDict);
            productDictWriter.close();
        } catch(IOException e) {
            System.err.println(e.getMessage());
            System.exit(1);
        }
    }

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
            String[] entries = dirFile.list();
            if (entries != null) {
                for(String s: entries){
                    File currentFile = new File(dirFile, s);
                    currentFile.delete();
                }
            }
            dirFile.delete();
        }
    }

    /**
     * Delete all index files (and only the files).
     * @param dir The directory to remove the index from.
     */
    private void removeFiles(String dir) {
        deleteFile(dir, tokenDictFileName);
        deleteFile(dir, productDictFileName);
        deleteFile(dir, reviewDataFileName);
        deleteFile(dir, productPostingListFileName);
        deleteFile(dir, tokenPostingListFileName);
    }

    /**
     * Delete a single file
     * @param dir The directory to delete from
     * @param fileName The file name
     */
    private void deleteFile(String dir, String fileName) {
        File f = new File(dir + File.separator + fileName);
        if (f.exists()) {
            f.delete();
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
