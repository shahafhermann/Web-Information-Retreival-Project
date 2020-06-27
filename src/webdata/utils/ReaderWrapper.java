package webdata.utils;

import java.io.*;

/**
 * A class representing a minimum heap of BufferedReaders
 */
public class ReaderWrapper implements Comparable<ReaderWrapper> {

	BufferedReader br;
	webdata.utils.Line currLine;

	/**
	 * Constructor
	 * @param br BufferedReader
	 * @param currLine The current line being read
	 */
	public ReaderWrapper(BufferedReader br, webdata.utils.Line currLine){
		this.br = br;
		this.currLine = currLine;
	}

	/**
	 * returns the current line
	 */
	public Line getCurrLine(){
		return currLine;
	}

	/**
	 * Advance the reader
	 * @return True if advanced, False otherwise
	 * @throws IOException
	 */
	public boolean advancePtr() throws IOException{
		String line = br.readLine();
		if (line != null){
			currLine = new Line(line);
			return true;
		}
		currLine = null;
		br.close();
		return false;
	}

	@Override
	public int compareTo(ReaderWrapper o) {
		return currLine.compareTo(o.currLine);
	}
}
