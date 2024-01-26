package invertibility;

import java.util.ArrayList;
import java.util.Arrays;

public class Main {
	public static void starTreeLength() {
		// set length process parameters
		final double LAMBDA = 1;
		final double MU = 0.7;
		final double NU = 3;
		final double PI0 = 0.5;
		final String ROOT = "01010101";

		// set the actual values for M, gamma, beta
		final int M = ROOT.length();
		final double GAMMA = LAMBDA / MU;
		final double BETA = Math.exp(LAMBDA - MU);

		// int[] N = { 1000, 10000, 100000, 1000000 };
		int[] N = { 100, 1000 };
		int numIter = N.length;

		// initialize arrays
		// logN is the number of leaves for each trial
		// gamma is the estimated gamma for each trial
		// beta is the estimated beta for each trial
		// m is the estimated m for each trial
		int[] logN = new int[numIter];

		final int NUM_SAMPLES = 50;
		double[][] gamma = new double[numIter][NUM_SAMPLES];
		double[][] beta = new double[numIter][NUM_SAMPLES];
		double[][] m = new double[numIter][NUM_SAMPLES];

		// run trials
		for (int trial = 0; trial < numIter; trial++) {
			int num = N[trial];
			System.out.println("Current trial: " + trial);

			logN[trial] = (int) Math.round(Math.log10(num));

			double[] gammaSample = new double[NUM_SAMPLES];
			double[] betaSample = new double[NUM_SAMPLES];
			double[] mSample = new double[NUM_SAMPLES];

			for (int sample = 0; sample < NUM_SAMPLES; sample++) {
				TreeSimul sampleTree = new TreeSimul(LAMBDA, MU, NU, PI0, ROOT, N[trial]);
				TreeLeaves sampleLeaves = sampleTree.toTreeLeaves();

				InvertLength inverted = new InvertLength(sampleLeaves);

				gammaSample[sample] = inverted.getGamma();
				betaSample[sample] = inverted.getBeta();

				Integer sampleM = inverted.getM();
				if (sampleM == null)
					mSample[sample] = Double.NaN;
				else
					mSample[sample] = sampleM;
			}

			gamma[trial] = trim(gammaSample);
			beta[trial] = trim(betaSample);
			m[trial] = trim(mSample);
		}

		String xAxis = "log N";
		// display charts
		System.out.println("gamma: " + GAMMA + " beta: " + BETA + " M: " + M);
		Chart.BoxWhiskerChart(xAxis, "gamma", gamma, GAMMA, logN);
		Chart.BoxWhiskerChart(xAxis, "beta", beta, BETA, logN);
		Chart.BoxWhiskerChart(xAxis, "M", m, M, logN);
	}

	public static void starTree1Mer() {
		// set length process parameters
		final double LAMBDA = 1;
		final double MU = 0.7;
		final double NU = 0.2;
		final double PI0 = 0.3;
		final String ROOT = "101011";
		final int M = ROOT.length();

		int temp = 0;
		for (char c : ROOT.toCharArray()) {
			if (c == '1')
				temp++;
		}
		final int A = temp;

		int[] N = { 1000, 10000 };
		int numIter = N.length;

		// initialize arrays
		// logN is the number of leaves for each trial
		int[] logN = new int[numIter];

		double[] nu = new double[numIter];
		double[] a = new double[numIter];

		double[] nuDeviation = new double[numIter];
		double[] aDeviation = new double[numIter];

		// run trials
		for (int trial = 0; trial < numIter; trial++) {
			int num = N[trial];
			System.out.println("Current trial: " + trial);

			logN[trial] = (int) Math.round(Math.log10(num));
			final int NUM_SAMPLES = 50;

			double[] nuSample = new double[NUM_SAMPLES];
			double[] aSample = new double[NUM_SAMPLES];

			for (int sample = 0; sample < NUM_SAMPLES; sample++) {
				TreeSimul sampleTree = new TreeSimul(LAMBDA, MU, NU, PI0, ROOT, N[trial]);
				TreeLeaves sampleLeaves = sampleTree.toTreeLeaves();

				Invert1Mer inverted = new Invert1Mer(LAMBDA, MU, PI0, M, sampleLeaves);

				nuSample[sample] = inverted.getNu();
				aSample[sample] = inverted.getA();
			}

			System.out.println("TRIMMING");
			nuSample = trim(nuSample);
			aSample = trim(aSample);

			nu[trial] = average(nuSample);
			nuDeviation[trial] = standardDeviation(nu[trial], nuSample);

			a[trial] = average(aSample);
			aDeviation[trial] = standardDeviation(a[trial], aSample);

			System.out.println("Trial: " + trial);
			System.out.println("nu: " + nu[trial] + " " + nuDeviation[trial]);
			// System.out.println("pi0: " + pi0[trial] + " " + pi0Deviation[trial]);
			System.out.println("a: " + a[trial] + " " + aDeviation[trial]);
		}

		String xAxis = "log N";
		// display charts
		System.out.println("nu: " + NU + " pi0: " + PI0 + " a: " + A);
		Chart.ErrorChart(xAxis, "nu", nu, nuDeviation, NU, logN);
		Chart.ErrorChart(xAxis, "a", a, aDeviation, A, logN);
	}

	public static void starTreeState(int method) {
		// set length process parameters
		final double LAMBDA = 1;
		final double MU = 0.4;
		final double NU = 0.2;
		final double PI0 = 0.3;
		final String ROOT = "11010111";
		// final String ROOT = "1101";
		final int M = ROOT.length();

		int[] N = { 1000 };
		int numIter = N.length;

		// initialize arrays
		// logN is the number of leaves for each trial
		int[] logN = new int[numIter];

		String[] root = new String[numIter];

		double[] diff = new double[numIter];
		double[] diffDeviation = new double[numIter];

		// run trials
		for (int trial = 0; trial < numIter; trial++) {
			int num = N[trial];
			System.out.println("Current trial: " + trial);

			logN[trial] = (int) Math.round(Math.log10(num));
			final int NUM_SAMPLES = 50;

			double[] diffSample = new double[NUM_SAMPLES];

			int[] count = new int[M]; // counts the number of 1's across all samples at each index

			for (int sample = 0; sample < NUM_SAMPLES; sample++) {
				TreeSimul sampleTree = new TreeSimul(LAMBDA, MU, NU, PI0, ROOT, N[trial]);
				TreeLeaves sampleLeaves = sampleTree.toTreeLeaves();

				InvertState inverted = new InvertState(LAMBDA, MU, NU, M, PI0, sampleLeaves, method);

				String sampleRootState = inverted.getRootState();
				for (int i = 0; i < M; i++) {
					if (sampleRootState.charAt(i) == '1')
						count[i]++;
					if (sampleRootState.charAt(i) != ROOT.charAt(i))
						diffSample[sample]++;
				}
			}

			diff[trial] = average(diffSample);
			diffDeviation[trial] = standardDeviation(diff[trial], diffSample);

			System.out.println("Trial: " + trial);
			System.out.println("diff: " + diff[trial]);

			root[trial] = "";
			for (int i = 0; i < M; i++) {
				root[trial] += (count[i] > NUM_SAMPLES / 2) ? '1' : '0';
				// System.out.println(i + "\t" + count[i]);
			}

			System.out.println("Trial: " + trial);
			System.out.println("Root: " + root[trial]);
		}

		String xAxis = "log N";
		// display charts
		System.out.println("REAL ROOT: " + ROOT);
		Chart.ErrorChart(xAxis, "difference", diff, diffDeviation, -1, logN);
	}

	public static void starTreePairwiseDistance() {
		// set length process parameters
		final double LAMBDA = 1;
		final double MU = 0.7;
		final double NU = 3;
		final double PI0 = 0.5;
		final String ROOT = "01010101";

		final double tw = 1;
		final double t1 = 1;
		final double t2 = 2;

		// set the actual values for M
		final int M = ROOT.length();

		int[] N = { 1000, 10000 };
		// int[] N = { 100, 1000 };
		int numIter = N.length;

		// initialize arrays
		// logN is the number of leaves for each trial
		int[] logN = new int[numIter];

		final int NUM_SAMPLES = 30;
		double[][] pwd = new double[numIter][NUM_SAMPLES];
		double[][] wd = new double[numIter][NUM_SAMPLES];

		// run trials
		for (int trial = 0; trial < numIter; trial++) {
			int num = N[trial];
			System.out.println("Current trial: " + trial);

			logN[trial] = (int) Math.round(Math.log10(num));

			double[] pwdSample = new double[NUM_SAMPLES];
			double[] wdSample = new double[NUM_SAMPLES];

			TreeSimul sampleTree = new TreeSimul(LAMBDA, MU, NU, PI0, ROOT, 0);
			for (int sample = 0; sample < NUM_SAMPLES; sample++) {
				int[] products = new int[N[trial]];
				for (int i = 0; i < N[trial]; i++)
					products[i] = sampleTree.pairwiseDistance(tw, t1, t2);

				InvertPairwiseDistance inverted = new InvertPairwiseDistance(LAMBDA, MU, 2, 3, M,
						products);

				pwdSample[sample] = inverted.getPairwiseDistance();
				wdSample[sample] = inverted.getAncestorDistance();
			}

			pwd[trial] = trim(pwdSample);
			wd[trial] = trim(wdSample);
		}

		String xAxis = "log N";
		// display charts
		System.out.println("exact: " + (t1 + t2) * MU);
		Chart.BoxWhiskerChart(xAxis, "pwd", pwd, (t1 + t2) * MU, logN);
		Chart.BoxWhiskerChart(xAxis, "wd", wd, tw * MU, logN);
	}

	private static double[] trim(double[] data) {
		Arrays.sort(data);

		int index = data.length;
		while (index > 0 && Double.isNaN(data[--index]))
			;

		if (Double.isNaN(data[index]))
			return new double[] {};

		System.out.println(index);

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
		int method = 6;
		System.out.println("Method " + method);
		// starTreeState(method);
		starTreePairwiseDistance();
	}
}
