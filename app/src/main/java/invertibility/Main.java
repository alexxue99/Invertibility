package invertibility;

import java.util.Arrays;
import java.util.Collections;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;
import java.util.Random;

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
		// set length process parameters
		double lambda = 1;
		double mu = 0.7;
		double nu = 3;
		double pi0 = 0.3;
		double time = 1;
		String root = "01010101";
		int rootVertex = 0;

		int[] Ns = { 1000, 10000, 100000, 1000000 };
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
		double[][] Cs = new double[3][numIter];

		double[] gUpper = new double[numIter];
		double[] bUpper = new double[numIter];
		double[] mUpper = new double[numIter];
		double[][] CsUpper = new double[3][numIter];

		double[] gLower = new double[numIter];
		double[] bLower = new double[numIter];
		double[] mLower = new double[numIter];
		double[][] CsLower = new double[3][numIter];

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
			final int NUM_SAMPLES = 50;

			double[] gSample = new double[NUM_SAMPLES];
			double[] bSample = new double[NUM_SAMPLES];
			double[] mSample = new double[NUM_SAMPLES];
			double[][] CsSample = new double[NUM_SAMPLES][3];

			for (int sample = 0; sample < NUM_SAMPLES; sample++) {
				TreeSimul example = new TreeSimul(lambda, mu, nu, pi0, rootVertex, root, edges, t);
				TreeLeaves exampleTree = example.toTreeLeaves();

				Invert exampleInverted = new Invert(exampleTree);

				gSample[sample] = exampleInverted.getGamma();
				bSample[sample] = exampleInverted.getBeta();
				mSample[sample] = exampleInverted.getM();

				CsSample[sample] = exampleInverted.getCs();
			}

			double[][] CsSampleTranspose = transpose(CsSample);

			double[] gTrimmed = trim(gSample);
			double[] bTrimmed = trim(bSample);
			double[] mTrimmed = trim(mSample);
			double[] CTrimmed1 = trim(CsSampleTranspose[0]);
			double[] CTrimmed2 = trim(CsSampleTranspose[1]);
			double[] CTrimmed3 = trim(CsSampleTranspose[2]);

			g[trial] = average(gTrimmed);
			gUpper[trial] = g[trial] + standardDeviation(g[trial], gTrimmed);
			gLower[trial] = g[trial] - standardDeviation(g[trial], gTrimmed);

			b[trial] = average(bTrimmed);
			bUpper[trial] = b[trial] + standardDeviation(b[trial], bTrimmed);
			bLower[trial] = b[trial] - standardDeviation(b[trial], bTrimmed);

			m[trial] = average(mTrimmed);
			mUpper[trial] = m[trial] + standardDeviation(m[trial], mTrimmed);
			mLower[trial] = m[trial] - standardDeviation(m[trial], mTrimmed);

			Cs[0][trial] = average(CTrimmed1);
			CsUpper[0][trial] = Cs[0][trial] + standardDeviation(Cs[0][trial], CTrimmed1);
			CsLower[0][trial] = Cs[0][trial] - standardDeviation(Cs[0][trial], CTrimmed1);

			Cs[1][trial] = average(CTrimmed2);
			CsUpper[1][trial] = Cs[1][trial] + standardDeviation(Cs[1][trial], CTrimmed2);
			CsLower[1][trial] = Cs[1][trial] - standardDeviation(Cs[1][trial], CTrimmed2);

			Cs[2][trial] = average(CTrimmed3);
			CsUpper[2][trial] = Cs[2][trial] + standardDeviation(Cs[2][trial], CTrimmed3);
			CsLower[2][trial] = Cs[2][trial] - standardDeviation(Cs[2][trial], CTrimmed3);

			System.out.println("Trial: " + trial);
			System.out.println("gamma: " + g[trial] + " " + gUpper[trial]);
			System.out.println("beta: " + b[trial] + " " + bUpper[trial]);
			System.out.println("M: " + m[trial] + " " + mUpper[trial]);
			System.out.println("C1: " + Cs[0][trial] + " " + CsUpper[0][trial]);
			System.out.println("C2: " + Cs[1][trial] + " " + CsUpper[1][trial]);
			System.out.println("C3: " + Cs[2][trial] + " " + CsUpper[2][trial]);
		}


		String xAxis = "log N";
		// display charts
		System.out.println("gamma: " + gamma + " beta: " + beta + " M: " + M);
		Chart.ErrorChart(xAxis, "gamma", g, gUpper, gLower, gamma, N);
		Chart.ErrorChart(xAxis, "beta", b, bUpper, bLower, beta, N);
		Chart.ErrorChart(xAxis, "M", m, mUpper, mLower, M, N);

		Chart.ErrorChart(xAxis, "C1", Cs[0], CsUpper[0], CsLower[0], cs[0], N);
		Chart.ErrorChart(xAxis, "C2", Cs[1], CsUpper[1], CsLower[1], cs[1], N);
		Chart.ErrorChart(xAxis, "C3", Cs[2], CsUpper[2], CsLower[2], cs[2], N);
	}

	private static double[][] transpose(double[][] data) {
		double[][] transpose = new double[data[0].length][data.length];

		for (int i = 0; i < data[0].length; i++)
			for (int j = 0; j < data.length; j++)
				transpose[i][j] = data[j][i];

		return transpose;
	}

	private static double[] trim(double[] data) {
		Arrays.sort(data);

		int index = data.length;
		while (index > 0 && Double.isNaN(data[--index]))
			;

		if (Double.isNaN(data[index]))
			return new double[] {};

		return Arrays.copyOfRange(data, 0, index + 1);
	}

	private static double average(double[] data) {
		double average = 0;
		int divisor = 0;

		for (double sample : data) {
			if (!Double.isNaN(sample))
				divisor++;
		}

		for (double sample : data) {
			if (!Double.isNaN(sample))
				average += (double) sample / divisor;
		}

		return average;
	}

	private static double standardDeviation(double average, double[] data) {
		double variance = 0;
		int divisor = data.length - 1;

		for (double sample : data) {
			variance += (sample - average) * (sample - average) / divisor;
		}

		return Math.sqrt(variance);
	}

	public static void main(String[] args) {
		// example1();
		example2();
	}
}
