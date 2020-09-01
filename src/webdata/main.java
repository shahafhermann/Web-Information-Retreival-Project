package webdata;
import java.text.SimpleDateFormat;
import java.util.*;

public class main {
    public static void main(String[] args) {
        String dir = "/Users/shahaf/Documents/UNI/אחזור מידע באינטרנט/Project/indexFiles";
        String file = "/Users/shahaf/Documents/UNI/אחזור מידע באינטרנט/Project/1000.txt";


        IndexWriter siw = new IndexWriter();
        siw.write(file, dir);
        IndexReader ir = new IndexReader(dir);
        ReviewSearch rs = new ReviewSearch(ir);

        ArrayList<String[]> queryTerms = new ArrayList<>();
        String[] queryTerms1 = {"treats", "that", "dogs", "love"};
        String[] queryTerms2 = {"lo", "dogs"};
        queryTerms.add(queryTerms1);
        queryTerms.add(queryTerms2);

        ArrayList<Vector<String>> queries = new ArrayList<>();
        for (String[] queryTerm : queryTerms) {
            Vector<String> query = new Vector<>(Arrays.asList(queryTerm));
            queries.add(query);
        }

        for (Vector<String> query : queries) {
            Enumeration<Integer> res;

            System.out.println("---- History Search ----");
            res = rs.vectorSpaceSearch(query.elements(), 5, true);
            System.out.println("Found: ");
            while (res.hasMoreElements()) {
                System.out.println(res.nextElement());
            }
        }

        System.out.println();

        for (Vector<String> query : queries) {
            Enumeration<Integer> res;

            System.out.println("---- Noisy Chanel with Bayes Rule ----");
            res = rs.vectorSpaceSearch(query.elements(), 5, false);
            System.out.println("Found: ");
            while (res.hasMoreElements()) {
                System.out.println(res.nextElement());
            }
        }

        siw.removeIndex(dir);

    }
}