package invertibility;

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

		// int[] Ns = { 1000, 10000, 100000, 1000000 };
		int[] N = { 1000, 10000, 100000 };
		int numIter = N.length;

		// initialize arrays
		// logN is the number of leaves for each trial
		// gamma is the estimated gamma for each trial
		// beta is the estimated beta for each trial
		// m is the estimated m for each trial
		int[] logN = new int[numIter];

		double[] gamma = new double[numIter];
		double[] beta = new double[numIter];
		double[] m = new double[numIter];

		double[] gammaDeviation = new double[numIter];
		double[] betaDeviation = new double[numIter];
		double[] mDeviation = new double[numIter];

		// run trials
		for (int trial = 0; trial < numIter; trial++) {
			int num = N[trial];
			System.out.println("Current trial: " + trial);

			logN[trial] = (int) Math.round(Math.log10(num));
			final int NUM_SAMPLES = 50;

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

			gammaSample = trim(gammaSample);
			betaSample = trim(betaSample);
			mSample = trim(mSample);

			gamma[trial] = average(gammaSample);
			gammaDeviation[trial] = standardDeviation(gamma[trial], gammaSample);

			beta[trial] = average(betaSample);
			betaDeviation[trial] = standardDeviation(beta[trial], betaSample);

			m[trial] = average(mSample);
			mDeviation[trial] = standardDeviation(m[trial], mSample);

			System.out.println("Trial: " + trial);
			System.out.println("gamma: " + gamma[trial] + " " + gammaDeviation[trial]);
			System.out.println("beta: " + beta[trial] + " " + betaDeviation[trial]);
			System.out.println("M: " + m[trial] + " " + mDeviation[trial]);
		}

		String xAxis = "log N";
		// display charts
		System.out.println("gamma: " + GAMMA + " beta: " + BETA + " M: " + M);
		Chart.ErrorChart(xAxis, "gamma", gamma, gammaDeviation, GAMMA, logN);
		Chart.ErrorChart(xAxis, "beta", beta, betaDeviation, BETA, logN);
		Chart.ErrorChart(xAxis, "M", m, mDeviation, M, logN);
	}

	public static void starTree1Mer() {
		// set length process parameters
		final double LAMBDA = 1;
		final double MU = 0.7;
		final double NU = 0.2;
		final double PI0 = 0.3;
		final String ROOT = "01010101";
		final int M = ROOT.length();

		int temp = 0;
		for (char c : ROOT.toCharArray()) {
			if (c == '1')
				temp++;
		}
		final int A = temp;

		int[] N = { 1000, 10000, 100000 };
		int numIter = N.length;

		// initialize arrays
		// logN is the number of leaves for each trial
		int[] logN = new int[numIter];

		double[] nu = new double[numIter];
		double[] pi0 = new double[numIter];
		double[] a = new double[numIter];

		double[] nuDeviation = new double[numIter];
		double[] pi0Deviation = new double[numIter];
		double[] aDeviation = new double[numIter];

		// run trials
		for (int trial = 0; trial < numIter; trial++) {
			int num = N[trial];
			System.out.println("Current trial: " + trial);

			logN[trial] = (int) Math.round(Math.log10(num));
			final int NUM_SAMPLES = 50;

			double[] nuSample = new double[NUM_SAMPLES];
			double[] pi0Sample = new double[NUM_SAMPLES];
			double[] aSample = new double[NUM_SAMPLES];

			for (int sample = 0; sample < NUM_SAMPLES; sample++) {
				TreeSimul sampleTree = new TreeSimul(LAMBDA, MU, NU, PI0, ROOT, N[trial]);
				TreeLeaves sampleLeaves = sampleTree.toTreeLeaves();

				Invert1Mer inverted = new Invert1Mer(LAMBDA, MU, M, sampleLeaves);

				nuSample[sample] = inverted.getNu();
				pi0Sample[sample] = inverted.getPi0();
				aSample[sample] = inverted.getA();
			}

			System.out.println("TRIMMING");
			nuSample = trim(nuSample);
			pi0Sample = trim(pi0Sample);
			aSample = trim(aSample);

			nu[trial] = average(nuSample);
			nuDeviation[trial] = standardDeviation(nu[trial], nuSample);

			pi0[trial] = average(pi0Sample);
			pi0Deviation[trial] = standardDeviation(pi0[trial], pi0Sample);

			a[trial] = average(aSample);
			aDeviation[trial] = standardDeviation(a[trial], aSample);

			System.out.println("Trial: " + trial);
			System.out.println("nu: " + nu[trial] + " " + nuDeviation[trial]);
			System.out.println("pi0: " + pi0[trial] + " " + pi0Deviation[trial]);
			System.out.println("a: " + a[trial] + " " + aDeviation[trial]);
		}

		String xAxis = "log N";
		// display charts
		System.out.println("nu: " + NU + " pi0: " + PI0 + " a: " + A);
		Chart.ErrorChart(xAxis, "nu", nu, nuDeviation, NU, logN);
		Chart.ErrorChart(xAxis, "pi0", pi0, pi0Deviation, PI0, logN);
		Chart.ErrorChart(xAxis, "a", a, aDeviation, A, logN);
	}

	public static void starTreeState() {
		// set length process parameters
		final double LAMBDA = 1;
		final double MU = 0.7;
		final double NU = 0.2;
		final double PI0 = 0.3;
		final String ROOT = "10000000";
		final int M = ROOT.length();

		int[] N = { 1, 10, 100, 1000, 10000, 100000 };
		int numIter = N.length;

		// initialize arrays
		// logN is the number of leaves for each trial
		int[] logN = new int[numIter];

		String[] root = new String[numIter];

		// run trials
		for (int trial = 0; trial < numIter; trial++) {
			int num = N[trial];
			System.out.println("Current trial: " + trial);

			logN[trial] = (int) Math.round(Math.log10(num));
			final int NUM_SAMPLES = 50;

			int[] count = new int[M]; // counts the number of 1's across all samples at each index

			for (int sample = 0; sample < NUM_SAMPLES; sample++) {
				TreeSimul sampleTree = new TreeSimul(LAMBDA, MU, NU, PI0, ROOT, N[trial]);
				TreeLeaves sampleLeaves = sampleTree.toTreeLeaves();

				InvertState inverted = new InvertState(LAMBDA, MU, NU, M, PI0, sampleLeaves);

				String sampleRootState = inverted.getRootState();
				for (int i = 0; i < M; i++) {
					if (sampleRootState.charAt(i) == '1')
						count[i]++;
				}
			}

			root[trial] = "";
			for (int i = 0; i < M; i++) {
				root[trial] += (count[i] > NUM_SAMPLES / 2) ? '1' : '0';
				System.out.println(i + "\t" + count[i]);
			}

			System.out.println("Trial: " + trial);
			System.out.println("Root: " + root[trial]);
		}

		System.out.println("ROOT: " + ROOT);
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
		starTreeState();
	}
}
