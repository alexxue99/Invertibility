package invertibility;

import java.util.Collections;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;

import invertibility.TreeSimul.Edge;

public class Main {
	/**
	 * Basic example to show how to input the lengths of the sequences at the leaves
	 * of a tree, invert the length process, and print out the estimated parameters.
	 */
	public static void example1() {
		System.out.println("Starting example 1...\n");

		List<Integer> seqLeavesLengths = new LinkedList<Integer>();
		for (int i = 0; i < 50; i++)
			seqLeavesLengths.add(1000 + i);

		TreeLeaves exampleTree = new TreeLeaves(seqLeavesLengths);
		Invert exampleInverted = new Invert(exampleTree);

		// Print out estimated parameters.
		System.out.println("Estimated gamma is " + exampleInverted.getGamma() + ".");
		System.out.println("Estimated beta is " + exampleInverted.getBeta() + ".");
		System.out.println("Estimated M is " + exampleInverted.getM() + ".\n");
	}

	/**
	 * Example to show how to run trials on a variable number of leaves, using the
	 * length process simulator TreeSimul.
	 */
	public static void example2() {
		System.out.println("Starting example 2... HI");
		// set length process parameters
		double lambda = 1;
		double mu = 0.7;
		double nu = 3;
		double pi0 = 0.3;
		double time = 1;
		String root = "01000000";
		int rootVertex = 0;

		int[] Ns = {10000, 100000, 1000000};
		int numIter = Ns.length;

		// initialize arrays
		// N is the number of leaves for each trial (so it will contain 100, 200, 300,
		// ..., 5000)
		// g is the estimated gamma for each trial
		// b is the estimated beta for each trial
		// m is the estimated m for each trial
		int[] N = new int[numIter];

		double[] g = new double[numIter];
		double[] b = new double[numIter];
		double[] m = new double[numIter];
		double[][] Cs = new double[numIter][3];

		double[] gUpper = new double[numIter];
		double[] bUpper = new double[numIter];
		double[] mUpper = new double[numIter];
		double[][] CsUpper = new double[numIter][3];

		double[] gLower = new double[numIter];
		double[] bLower = new double[numIter];
		double[] mLower = new double[numIter];
		double[][] CsLower = new double[numIter][3];

		// set the actual values for M, gamma, beta, c1, c2, and c3
		int M = root.length();
		double gamma = lambda / mu;
		double beta = Math.exp((lambda - mu) * time);

		double[] cs = new double[3];
		cs[0] = beta * M;
		cs[1] = beta * (2 * gamma - beta * gamma - beta) / (1 - gamma) * M;
		cs[2] = 2 * beta * (beta * beta * gamma * gamma + beta * beta * gamma + beta * beta
				- 3 * beta * gamma * gamma - 3 * beta * gamma + 3 * gamma * gamma) / ((1 - gamma) * (1 - gamma)) * M;

		// run trials
		for (int trial = 0; trial < numIter; trial++) {
			int num = Ns[trial];
			System.out.println("Current trial: " + trial);
			// set the directed edges for the star tree
			LinkedList<Edge> edges = new LinkedList<Edge>();
			double[] t = new double[num];
			for (int i = 0; i < num; i++) {
				edges.add(new Edge(0, i + 1));
				t[i] = time;
			}

			N[trial] = (int) Math.round(Math.log10(num));

			TreeSimul example = new TreeSimul(lambda, mu, nu, pi0, rootVertex, root, edges, t);
			TreeLeaves exampleTree = example.toTreeLeaves();

			long[][] lengths = new long[3][num];
			List<Integer> seqLeavesLengths = exampleTree.getSeqLeavesLengths();
			Iterator<Integer> leafIter = seqLeavesLengths.iterator();

			for (int j = 0; j < num; j++) {
				lengths[0][j] = (long) leafIter.next();
				lengths[1][j] = lengths[0][j] * (lengths[0][j] - 1);
				lengths[2][j] = lengths[1][j] * (lengths[0][j] - 2);
			}

			double[][] avg_std = new double[3][2];
			for (int i = 0; i < 3; i++) {
				avg_std[i] = average_standardDeviation(lengths[i]);
			}

			double[] partials = new double[3];
			double[] partialsUpper = new double[3];
			double[] partialsLower = new double[3];

			for (int i = 0; i < 3; i++) {
				partials[i] = avg_std[i][0];
				partialsUpper[i] = avg_std[i][0] + avg_std[i][1];
				partialsLower[i] = avg_std[i][0] - avg_std[i][1];
			}
		
			Invert exampleInverted = new Invert(partials);

			g[trial] = exampleInverted.getGamma();
			b[trial] = exampleInverted.getBeta();
			m[trial] = exampleInverted.getM();
			Cs[trial] = exampleInverted.getCs();

			exampleInverted = new Invert(partialsUpper);

			gUpper[trial] = exampleInverted.getGamma();
			bUpper[trial] = exampleInverted.getBeta();
			mUpper[trial] = exampleInverted.getM();
			CsUpper[trial] = exampleInverted.getCs();

			exampleInverted = new Invert(partialsLower);

			gLower[trial] = exampleInverted.getGamma();
			bLower[trial] = exampleInverted.getBeta();
			mLower[trial] = exampleInverted.getM();
			CsLower[trial] = exampleInverted.getCs();
		}

		// for (int num = stepSize, trial = 0; num <= maxN; num+= stepSize, trial++) {
		// System.out.println("trial " + trial + ": " + g[trial]);
		// }

		// display charts
		Chart.ErrorChart("Gamma vs. log N", "gamma", g, gUpper, gLower, gamma, N);
		Chart.ErrorChart("Beta vs. log N", "beta", b, bUpper, bLower, beta, N);
		Chart.ErrorChart("M vs. log N", "M", m, mUpper, mLower, M, N);

		double[][] CsTranspose = transpose(Cs);
		double[][] CsUpperTranspose = transpose(CsUpper);
		double[][] CsLowerTranspose = transpose(CsLower);

		Chart.ErrorChart("C1 vs. log N", "C1", CsTranspose[0], CsUpperTranspose[0], CsLowerTranspose[0], cs[0], N);
		Chart.ErrorChart("C2 vs. log N", "C2", CsTranspose[1], CsUpperTranspose[1], CsLowerTranspose[1], cs[1], N);
		Chart.ErrorChart("C3 vs. log N", "C3", CsTranspose[2], CsUpperTranspose[2], CsLowerTranspose[2], cs[2], N);
	}

	public static double[][] transpose(double[][] data) {
		double[][] transpose = new double[data[0].length][data.length];

		for (int i = 0; i < data[0].length; i++)
			for (int j = 0; j < data.length; j++)
				transpose[i][j] = data[j][i];

		return transpose;
	}

	public static double[] average_standardDeviation(long[] data) {
		double average = 0;
		int divisor = data.length;

		for (long sample : data) {
			average += (double) sample / divisor;
		}

		double variance = 0;
		divisor--;

		for (long sample : data) {
			variance += (sample - average) * (sample - average) / divisor;
		}

		System.out.println("AVG: " + average);
		System.out.println("Variance: " + Math.sqrt(variance));

		return new double[] { average, Math.sqrt(variance) };
	}

	public static void main(String[] args) {
		// example1();
		example2();

	}
}
