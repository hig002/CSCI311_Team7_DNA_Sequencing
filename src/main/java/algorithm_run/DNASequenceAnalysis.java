package algorithm_run;

import java.io.*;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.Arrays;
import java.util.LinkedHashMap;
import java.util.Map;
import java.util.Scanner;

public class DNASequenceAnalysis {

	// Method to read a file into a String
	private static String readFileAsString(String fileName) throws IOException {
		// reads all bytes from a file and returns them as a string
		return new String(Files.readAllBytes(Paths.get(fileName)));
	}

	static Map<String, String> readFasta(String filename) {
		// Reads a FASTA file and returns a map of sequence headers to sequences
		// filename is the path to the FASTA file

		Map<String, String> sequencesMap = new LinkedHashMap<>();
		try {
			BufferedReader reader = new BufferedReader(new FileReader(filename));
			String line;
			StringBuilder sequence = new StringBuilder();
			String header = null;
			while ((line = reader.readLine()) != null) {
				if (line.startsWith(">")) {
					// If it's a header line, store previous sequence if any
					if (header != null && sequence.length() > 0) {
						sequencesMap.put(header, sequence.toString());
						sequence.setLength(0); // Clear the StringBuilder
					}
					header = line.substring(1); // Extract header (remove '>')
				} else {
					// If it's a sequence line, append to StringBuilder
					sequence.append(line);
				}
			}
			// Store the last sequence
			if (header != null && sequence.length() > 0) {
				sequencesMap.put(header, sequence.toString());
			}
			reader.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
		return sequencesMap;
	}

	/*
	 ********************************************************************************************************************
	 * ALGORITHM IMPLEMENTATIONS BELOW
	 * *****************************************************************************
	 * *************************************
	 */

	// Needleman-Wunsch Algorithm
	private static int nwAlgo(String query, String data, int matchScore, int mismatchScore, int gapScore) {

		// Implementation of the Needleman-Wunsch algorithm
		// query is the query sequence
		// data is the data sequence
		// matchScore is the score for a match
		// mismatchScore is the score for a mismatch
		// gapScore is the score for a gap

		// Construct the alignment matrix
		int[][] matrix = new int[query.length() + 1][data.length() + 1];
		// Fill in the base cases
		// Fill in the top row
		for (int col = 1; col < matrix[0].length; col++) {
			matrix[0][col] = col * gapScore;
		}
		// Fill in the left column
		for (int row = 1; row < matrix.length; row++) {
			matrix[row][0] = row * gapScore;
		}

		// Fill out the matrix table
		for (int row = 1; row < matrix.length; row++) {
			for (int col = 1; col < matrix[row].length; col++) {
				// Compare Sequences
				// Calculate Diagonal Score
				int diagScore = matrix[row - 1][col - 1];
				if (query.charAt(row - 1) == data.charAt(col - 1)) {
					diagScore += matchScore;
				} else {
					diagScore += mismatchScore;
				}
				// Calculate Top Score
				int topScore = matrix[row - 1][col] + gapScore;
				// Calculate Side Score
				int sideScore = matrix[row][col - 1] + gapScore;

				// Compare and find the largest number to use as the score for the matrix
				int finalScore = Math.max(topScore, sideScore);
				finalScore = Math.max(finalScore, diagScore);
				matrix[row][col] = finalScore;
			}
		}

		return matrix[matrix.length - 1][matrix[0].length - 1];
	}

	// Longest Common Substring
	private static String LCSubStr(char[] X, char[] Y, int m, int n) {
		// finds the longest common substring between two strings
		// X and Y are the two strings
		// m and n are the lengths of X and Y respectively

		int LCStuff[][] = new int[m + 1][n + 1];
		int result = 0; // To store length of the LCS
		int endPos = -1; // To store the ending position of the LCS in X

		for (int i = 0; i <= m; i++) {
			for (int j = 0; j <= n; j++) {
				if (i == 0 || j == 0)
					LCStuff[i][j] = 0;
				else if (X[i - 1] == Y[j - 1]) {
					LCStuff[i][j] = LCStuff[i - 1][j - 1] + 1;
					if (LCStuff[i][j] > result) {
						result = LCStuff[i][j];
						endPos = i; // Remember the end position of LCS in X
					}
				} else
					LCStuff[i][j] = 0;
			}
		}

		// If there's no common substring, return an appropriate message.
		if (endPos == -1) {
			return "No Common Substring Found!";
		}

		// Extracting the longest common substring from X
		String lcs = String.valueOf(Arrays.copyOfRange(X, endPos - result, endPos));
		return lcs;
	}

	// Levenshtein Distance
	private static int levenshteinDistance(String s1, String s2) {
		// finds the minimum number of single-character edits required to change one
		// word into another
		// s1 and s2 are the two strings
		s1 = s1.toLowerCase();
		s2 = s2.toLowerCase();

		int[][] save = new int[s1.length() + 1][s2.length() + 1];

		// initialize the first row and column
		for (int i = 0; i <= s1.length(); i++) {
			save[i][0] = i;
		}

		for (int j = 0; j <= s2.length(); j++) {
			save[0][j] = j;
		}

		// dynamic programming solution for the rest of the table
		for (int i = 1; i <= s1.length(); i++) {
			for (int j = 1; j <= s2.length(); j++) {
				int cost = 0;
				if (s1.charAt(i - 1) != s2.charAt(j - 1)) {
					cost = 1;
				}
				save[i][j] = Math.min(Math.min(save[i - 1][j] + 1, save[i][j - 1] + 1), save[i - 1][j - 1] + cost);
			}
		}

		return save[s1.length()][s2.length()];
	}

	// LCS with Dynamic Programming
	private static int lcsDP(String sequence, String query, int m, int n, int[][] dp) {
		// A Top-Down DP implementation of LCS problem

		// Returns length of LCS for X[0..m-1], Y[0..n-1] static int
		if (m == 0 || n == 0)
			return 0;

		if (dp[m][n] != -1)
			return dp[m][n];

		if (sequence.charAt(m - 1) == query.charAt(n - 1)) {
			dp[m][n] = 1 + lcsDP(sequence, query, m - 1, n - 1, dp);
			return dp[m][n];
		}

		dp[m][n] = Math.max(lcsDP(sequence, query, m, n - 1, dp),
				lcsDP(sequence, query, m - 1, n, dp));

		return dp[m][n];
	}

	public static void main(String[] args) {
		// program begins execution here
		try {
			Scanner scanner = new Scanner(System.in);
			String continueChoice = "";

			do {
				System.out.println("Enter the file path for the query sequence:");
				String queryFilePath = scanner.nextLine();
				String querySequence = readFileAsString("src/main/resources/" + queryFilePath);

				System.out.println("Enter the file path for the FASTA file:");
				String fastaFilePath = scanner.nextLine();
				Map<String, String> fastaSequences = readFasta("src/main/resources/" + fastaFilePath);

				System.out.println(
						"Choose the algorithm:\n1. Needleman-Wunsch\n2. Longest Common Substring\n3. Levenshtein Distance\n4. LCS with Dynamic Programming");
				int choice = scanner.nextInt();
				scanner.nextLine();

				runAlgorithm(choice, fastaSequences, querySequence, queryFilePath, fastaFilePath);

				System.out.println("Do you want to continue? (Y/N)");
				continueChoice = scanner.nextLine().toLowerCase();

			} while (continueChoice.equals("y"));

		} catch (IOException e) {
			System.err.println("An error occurred while reading the files.");
			e.printStackTrace();
		}
	}

	private static void runAlgorithm(int choice, Map<String, String> fastaSequences, String querySequence,
									 String queryFilePath, String fastaFilePath) throws IOException {
		// runs an algorithm based on the user's choice and outputs the results

		switch (choice) {
			case 1:
				int maxScore = Integer.MIN_VALUE;
				String index = "";

				Scanner scanner = new Scanner(System.in);
				boolean isValidInput = false;
				int matchScore = 0;
				int mismatchScore = 0;
				int gapScore = 0;
				while (!isValidInput) {
					try {
						System.out.print("Enter match score value: ");
						matchScore = scanner.nextInt();
						scanner.nextLine();
						System.out.print("Enter mismatch score  value: ");
						mismatchScore = scanner.nextInt();
						scanner.nextLine();
						System.out.print("Enter gap score value: ");
						gapScore = scanner.nextInt();
						scanner.nextLine();
						if (matchScore < 0 || mismatchScore > 0 || gapScore > 0) {
							System.out.println("Please enter a valid input");
						} else {
							isValidInput = true;
						}
					}
					catch (Exception e) {
						System.out.println("Please enter a valid input");
						scanner.nextLine();
					}
				}

				for (Map.Entry<String, String> entry : fastaSequences.entrySet()) {
					String sequenceId = entry.getKey();
					String sequence = entry.getValue();
					int score = nwAlgo(querySequence, sequence, matchScore, mismatchScore, gapScore);
					if (score > maxScore) {
						maxScore = score;
						index = sequenceId;
					}
				}
				System.out.println("Sequence with the highest score: " + index + "\nScore: " + maxScore);
				break;
			case 2:
				int maxLcsLength = 0;
				index = "";
				String maxLcs = "";

				for (Map.Entry<String, String> entry : fastaSequences.entrySet()) {
					String sequenceId = entry.getKey();
					String sequence = entry.getValue();
					String lcs = LCSubStr(querySequence.toCharArray(), sequence.toCharArray(), querySequence.length(),
							sequence.length());
					if (lcs.length() > maxLcsLength) {
						maxLcsLength = lcs.length();
						index = sequenceId;
						maxLcs = lcs;
					}
				}
				System.out.println("Sequence with the longest common substring: " + index + "\nLength: " + maxLcsLength);
				System.out.println("Longest Common Substring: " + maxLcs);

				break;
			case 3:
				int minDistance = Integer.MAX_VALUE;
				String minSequenceId = "";

				for (Map.Entry<String, String> entry : fastaSequences.entrySet()) {
					String sequenceId = entry.getKey();
					String sequence = entry.getValue();
					int distance = levenshteinDistance(querySequence, sequence);
					if (distance < minDistance) {
						minDistance = distance;
						minSequenceId = sequenceId;
					}
				}



				System.out.println(
						"Sequence with the smallest Levenshtein Distance: " + minSequenceId + "\nDistance: " + minDistance);

				break;
			case 4:
				maxLcsLength = 0;
				index = "";

				for (Map.Entry<String, String> entry : fastaSequences.entrySet()) {
					String sequenceId = entry.getKey();
					String sequence = entry.getValue();
					int m = sequence.length();
					int n = querySequence.length();
					int[][] dp = new int[m + 1][n + 1];
					for (int i = 0; i <= m; i++) {
						Arrays.fill(dp[i], -1);
					}
					int lcsLength = lcsDP(sequence, querySequence, m, n, dp);
					if (lcsLength > maxLcsLength) {
						maxLcsLength = lcsLength;
						index = sequenceId;
					}
				}
				System.out.println("Sequence with the longest common subsequence: " + index + "\nLength: " + maxLcsLength);
				break;
			default:
				System.out.println("Invalid choice");
				break;
		}
	}
}