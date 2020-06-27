package webdata;

import webdata.utils.Line;
import webdata.utils.ReaderWrapper;

import java.io.*;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.*;

/**
 * A class for sorting the dataset
 */
public class Sorter {
    /* Data */
    private ArrayList<String> tokensArray;
    private ArrayList<String> productIdsArray;
    private HashMap<String, Integer> tokensMap;
    private HashMap<String, Integer> productIdMap;
    private int numOfReviews = 0;

    private ArrayList<webdata.utils.Line> tokenLines = new ArrayList<>();
    private ArrayList<webdata.utils.Line> productIdLines = new ArrayList<>();

    /* String constants */
    private static final String SPLIT_TOKENS_REGEX = "[^A-Za-z0-9]+";

    /* File paths to save the terms lists */
    private String tmpDir;
    private static final int NUM_OF_REVIEWS_PER_FILE = 10000;
    private static final String SORT_TEMP_TOKEN_FILE_NAME = "t_%d_%d.txt";
    private static final String SORT_TEMP_PRODUCT_FILE_NAME = "p_%d_%d.txt";
    private int numOfTempFiles = 0;
    private final int M = 1000;

    /**
     * Constructor
     * @param tokensArray ArrayList of token Strings
     * @param productIdsArray ArrayList of product id Strings
     * @param tmpDir Directory of temp files
     */
    Sorter(ArrayList<String> tokensArray, ArrayList<String> productIdsArray, String tmpDir) {
        Collections.sort(tokensArray);
        this.tokensArray = tokensArray;
        this.tokensMap = buildHashMap(tokensArray);

        Collections.sort(productIdsArray);
        this.productIdsArray = productIdsArray;
        this.productIdMap = buildHashMap(productIdsArray);

        this.tmpDir = tmpDir;
    }

    /**
     * Clear data members
     */
    void clear() {
        tokensMap = new HashMap<>();
        productIdMap = new HashMap<>();
    }

    /**
     * Build a hash table from the given ArrayList, where the String at index i is mapped to i.
     * @param toConvert The ArrayList to build from
     * @return The newly created HashMap
     */
    private static HashMap<String, Integer> buildHashMap(ArrayList<String> toConvert) {
        HashMap<String, Integer> hmap = new HashMap<>();
        for (int i = 0; i < toConvert.size(); ++i) {
            hmap.put(toConvert.get(i), i);
        }
        return hmap;
    }

    /**
     * Get the ArrayList of token Strings
     */
    ArrayList<String> getTokensArray() { return tokensArray; }

    /**
     * Get the ArrayList of product id Strings
     */
    ArrayList<String> getProductIdsArray() { return productIdsArray; }

    /**
     * Break a text to all it's tokens (alphanumeric).
     * @param text The text to break
     */
    private void breakText(String text) {
        ArrayList<String> tokens = new ArrayList<>(Arrays.asList(text.split(SPLIT_TOKENS_REGEX)));
        Collections.sort(tokens);
        String prevToken = "";
        int freq = 1;
        for (String token: tokens) {
            if (!token.isEmpty()) {
                if (!token.equals(prevToken)) {
                    if (!prevToken.isEmpty()) {
                        tokenLines.add(createLine(tokensMap.get(prevToken), freq));
                    }
                    prevToken = token;
                    freq = 1;
                } else {
                    ++freq;
                }
            }
        }

        if (!prevToken.isEmpty()) {
            tokenLines.add(createLine(tokensMap.get(prevToken), freq));
        }
    }

    /**
     * Create a new line object
     * @param term term of the line
     * @param freq frequency of the term
     * @return the new line object
     */
    private webdata.utils.Line createLine(int term, int freq){
        return new webdata.utils.Line(term + "#" + numOfReviews + "#" + freq);
    }

    /**
     * Parse the file
     * @param inputFile The file to parse
     */
    void firstPhase(String inputFile) {
        try (BufferedReader reader = new BufferedReader(new FileReader(new File(inputFile)))){
            String line = reader.readLine();
            String textBuffer = "";
            boolean textFlag = false;
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
                    if (numOfReviews % (NUM_OF_REVIEWS_PER_FILE + 1)  == 0) {
                        createTempFiles();
                    }
                    productIdLines.add(createLine(productIdMap.get(line.substring(19)), 1));
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
            createTempFiles();
        } catch (IOException e) {
            System.err.println(e.getMessage());
            System.exit(1);
        }
    }

    /**
     * Create temp files for the sorting process
     */
    private void createTempFiles() {
        Collections.sort(tokenLines);
        Collections.sort(productIdLines);
        writeMBlocks(tokenLines, SORT_TEMP_TOKEN_FILE_NAME);
        writeMBlocks(productIdLines, SORT_TEMP_PRODUCT_FILE_NAME);
        ++numOfTempFiles;
        tokenLines = new ArrayList<>();
        productIdLines = new ArrayList<>();
    }

    /**
     * Writes the blocks to a temp file with number tempNumber.
     * @param blocks The blocks to write.
     * @param fileName The file name to be used as the temp file
     */
    private void writeMBlocks(ArrayList<webdata.utils.Line> blocks, String fileName){
        try (BufferedWriter writer = new BufferedWriter(new FileWriter(
                new File(Paths.get(tmpDir, String.format(fileName, 0, numOfTempFiles)).toString())))){
            // Writes the block lines to the temp file
            for (webdata.utils.Line line:blocks) {
                writer.write(line.toString());
                writer.newLine();
            }
        } catch (IOException e) {
            System.err.println(e.getMessage());
        }
    }

    /**
     * This method performs merge sort on the first column of the table in the given in file, R(A,B,C) or S(A,D,E),
     * and writes the output to the out file.
     * @param in The pathname of the file to read from.
     * @param outTokens The pathname of the token file to write to.
     * @param outProducts The pathname of the product file to write to.
     */
    public void sort(String in, String outTokens, String outProducts) {
        firstPhase(in);
        clear();
        secondPhase(outTokens, tmpDir, numOfTempFiles, SORT_TEMP_TOKEN_FILE_NAME);
        secondPhase(outProducts, tmpDir, numOfTempFiles, SORT_TEMP_PRODUCT_FILE_NAME);
    }

    private void copyOut(File in, String out) {
        try (BufferedReader reader = new BufferedReader(new FileReader(in))) {
            try (BufferedWriter writer = new BufferedWriter(new FileWriter(new File(out)))) {
                String line = reader.readLine();
                while (line != null) {
                    writer.write(line);
                    writer.newLine();
                    line = reader.readLine();
                }
            }
        } catch (IOException e) {
            System.err.println(e.getMessage());
        }
    }

    /**
     * This method performs the second phase of the two phase sort algorithm. It reads all lines from
     * numberOfTempFiles temp files into the output file.
     * @param out The output file name.
     * @param tmpPath The path of the temp files.
     * @param numberOfTempFiles The number of temp files.
     */
    private void secondPhase(String out, String tmpPath, int numberOfTempFiles, String fileName){
        double mergeSteps = Math.ceil(Math.log(numberOfTempFiles) / Math.log(M));
        double numOfFiles = numberOfTempFiles;
        String outputFileName = out;
        if (mergeSteps == 0) {
            String in = String.format(fileName, 0, 0);
            copyOut(Paths.get(tmpDir, in).toFile(), out);
            try {
                Files.deleteIfExists(Paths.get(tmpDir, in));
            } catch (IOException e) {
                System.err.println(e.getMessage());
            }
            return;
        }
        for ( int currStep = 1; currStep <= mergeSteps; currStep++) {
            numOfFiles = Math.ceil(numOfFiles / M);

            int start = 0;
            int end = start + M;
            for (int outputFileIndex = 0; outputFileIndex < numOfFiles; outputFileIndex++){
                if ((currStep == mergeSteps) && (outputFileIndex == (numOfFiles - 1))) {
                    outputFileName = out;
                } else {
                    outputFileName = Paths.get(tmpDir, String.format(fileName, currStep, outputFileIndex)).toString();
                }
                mergeOnce(outputFileName, tmpPath, start, end , fileName, currStep - 1);
                start = end;
                end += M;
            }
        }
    }

    /**
     * Merge one chunk of files, starting at 'start' and ending at 'end', to a single merged file.
     * @param out A single sorted file for the chunk
     * @param tmpPath the path of the temp files directory
     * @param start first file to sort
     * @param end last file to sort
     * @param fileName The final sorted file name
     * @param prevStep
     */
    private void mergeOnce(String out, String tmpPath, int start, int end, String fileName, int prevStep) {
        PriorityQueue<ReaderWrapper> heapOfReaders = new PriorityQueue<>();

        try (BufferedWriter writer = new BufferedWriter(new FileWriter(new File(out)))) {
            this.initializeReaders(heapOfReaders, tmpPath, start, end, fileName, prevStep);

            // While there are more lines left, write the next minimal line to the output
            while (!heapOfReaders.isEmpty()) {
                ReaderWrapper minReader = heapOfReaders.poll();
                writer.write(minReader.getCurrLine().toString());
                writer.newLine();
                if (minReader.advancePtr()){
                    heapOfReaders.add(minReader);
                }
            }
            writer.flush();
            deleteTempFiles(start, end, fileName, prevStep);
        } catch (IOException e) {
            System.err.println(e.getMessage());
        }
    }

    /**
     * Initialize an array of readers to read from the files.
     * @param heapOfReaders An object wrapping readers from file as a minimum heap
     * @param tmpPath the path of the temp files directory
     * @param endingFileIndex last file to read from
     * @param startingFileIndex first file to read from
     * @param fileName The final sorted file name
     * @throws IOException
     */
    private void initializeReaders(PriorityQueue<ReaderWrapper> heapOfReaders, String tmpPath,
                                   int startingFileIndex, int endingFileIndex, String fileName, int prevStep)
            throws IOException {
        for (int i = startingFileIndex; i < endingFileIndex; i++) {
            Path filePath = Paths.get(tmpPath, String.format(fileName, prevStep, i));
            if (Files.exists(filePath)){
                BufferedReader br = new BufferedReader(new FileReader(filePath.toFile()));
                String line = br.readLine();
                if (line != null){
                    heapOfReaders.add(new ReaderWrapper(br, new Line(line)));
                }
            }else{
                break;
            }
        }
    }

    /**
     * Delete the temp files created in the sorting process
     * @param start first file to delete
     * @param end last file to delete
     * @param fileName The file name to delete. this is a template
     * @param prevStep this is required for the file name's template
     */
    private void deleteTempFiles(int start, int end, String fileName, int prevStep) {
        try {
            for (;start < end; ++start) {
                Files.deleteIfExists(Paths.get(tmpDir, String.format(fileName, prevStep, start)));
            }
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

}