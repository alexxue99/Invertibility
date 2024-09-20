package invertibility;

import java.util.Arrays;

public class Main {
	/**
	 * Tests the InvertLength class. Simulates the length process for multiple
	 * trials for each value of N and plots the estimated gamma, beta, and M values.
	 * 
	 * @param treeSimul   the TreeSimul object containing the TKF91 process
	 *                    parameters
	 * @param N           an array containing the values of N to simulate
	 * @param NUM_SAMPLES the number of trials to use for each N
	 */
	public static void starTreeLength(TreeSimul treeSimul, int[] N,
			int NUM_SAMPLES) {
		double GAMMA = treeSimul.getLambda() / treeSimul.getMu();
		double BETA = Math.exp(treeSimul.getLambda() - treeSimul.getMu());

		double[][] gamma = new double[N.length][NUM_SAMPLES];
		double[][] beta = new double[N.length][NUM_SAMPLES];
		double[][] m = new double[N.length][NUM_SAMPLES];

		// run trials
		for (int trial = 0; trial < N.length; trial++) {
			double[] gammaSample = new double[NUM_SAMPLES];
			double[] betaSample = new double[NUM_SAMPLES];
			double[] mSample = new double[NUM_SAMPLES];

			for (int sample = 0; sample < NUM_SAMPLES; sample++) {
				treeSimul.runTKF91Process(N[trial]);
				TreeLeaves sampleLeaves = treeSimul.toTreeLeaves();

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
		Chart.BoxWhiskerChart(xAxis, "gamma", gamma, GAMMA, N, false);
		Chart.BoxWhiskerChart(xAxis, "beta", beta, BETA, N, false);
		Chart.BoxWhiskerChart(xAxis, "M", m, treeSimul.getM(), N, false);
	}

	/**
	 * Tests the Invert1Mer class. Simulates the 1mer process for multiple
	 * trials for each value of N and plots the estimated nu and a values.
	 * 
	 * @param treeSimul   the TreeSimul object containing the TKF91 process
	 *                    parameters
	 * @param N           an array containing the values of N to simulate
	 * @param NUM_SAMPLES the number of trials to use for each N
	 */
	public static void starTree1Mer(TreeSimul treeSimul, int[] N,
			int NUM_SAMPLES) {
		int numOnes = 0;
		for (char c : treeSimul.getRoot().toCharArray()) {
			if (c == '1')
				numOnes++;
		}
		final int A = numOnes;

		int[] logN = new int[N.length];
		double[][] nu = new double[N.length][NUM_SAMPLES];
		double[][] a = new double[N.length][NUM_SAMPLES];

		// run trials
		for (int trial = 0; trial < N.length; trial++) {
			int num = N[trial];

			logN[trial] = (int) Math.round(Math.log10(num));

			double[] nuSample = new double[NUM_SAMPLES];
			double[] aSample = new double[NUM_SAMPLES];

			for (int sample = 0; sample < NUM_SAMPLES; sample++) {
				treeSimul.runTKF91Process(N[trial]);
				TreeLeaves sampleLeaves = treeSimul.toTreeLeaves();

				Invert1Mer inverted = new Invert1Mer(treeSimul.getLambda(), treeSimul.getMu(), treeSimul.getPi0(),
						treeSimul.getM(), sampleLeaves);

				nuSample[sample] = inverted.getNu();
				aSample[sample] = inverted.getA();
			}

			nu[trial] = trim(nuSample);
			a[trial] = trim(aSample);
		}

		String xAxis = "log N";
		// display charts
		Chart.BoxWhiskerChart(xAxis, "nu", nu, treeSimul.getNu(), N, false);
		Chart.BoxWhiskerChart(xAxis, "a", a, A, N, false);

	}

	/**
	 * Tests the InvertState class. Simulates the TKF91 process for multiple
	 * trials for each value of N and plots the difference between the estimated
	 * root state and the actual root state.
	 * 
	 * @param treeSimul   the TreeSimul object containing the TKF91 process
	 *                    parameters
	 * @param N           an array containing the values of N to simulate
	 * @param NUM_SAMPLES the number of trials to use for each N
	 */
	public static void starTreeState(TreeSimul treeSimul, int[] N, int NUM_SAMPLES) {
		String[] root = new String[N.length];

		double[][] diff = new double[N.length][NUM_SAMPLES];

		// run trials
		for (int trial = 0; trial < N.length; trial++) {
			double[] diffSample = new double[NUM_SAMPLES];

			int[] count = new int[treeSimul.getM()]; // counts the number of 1's across all samples at each index

			for (int sample = 0; sample < NUM_SAMPLES; sample++) {
				treeSimul.runTKF91Process(N[trial]);
				TreeLeaves sampleLeaves = treeSimul.toTreeLeaves();

				InvertState inverted = new InvertState(treeSimul.getLambda(), treeSimul.getMu(), treeSimul.getNu(),
						treeSimul.getM(), treeSimul.getPi0(), sampleLeaves);

				String sampleRootState = inverted.getRootState();
				for (int i = 0; i < treeSimul.getM(); i++) {
					if (sampleRootState.charAt(i) == '1')
						count[i]++;
					if (sampleRootState.charAt(i) != treeSimul.getRoot().charAt(i))
						diffSample[sample]++;
				}
			}

			diff[trial] = diffSample;
			root[trial] = "";
			for (int i = 0; i < treeSimul.getM(); i++) {
				root[trial] += (count[i] > NUM_SAMPLES / 2) ? '1' : '0';
			}
		}

		String xAxis = "log N";
		// display charts
		Chart.BoxWhiskerChart(xAxis, "difference", diff, -1, N, true);
	}

	/**
	 * Tests the InvertPairwiseDistance class. Simulates the TKF91 process for
	 * multiple
	 * trials for each value of N and plots the estimated pairwise distance between
	 * u and v and the distance between the root and w.
	 * 
	 * @param treeSimul   the TreeSimul object containing the TKF91 process
	 *                    parameters
	 * @param N           an array containing the values of N to simulate
	 * @param NUM_SAMPLES the number of trials to use for each N
	 */
	public static void starTreePairwiseDistance(TreeSimul treeSimul, int[] N, int NUM_SAMPLES) {
		final double tw = 1; // distance from root to w
		final double t1 = 2; // distance from w to u
		final double t2 = 3; // distance from w to v
		final double tu = tw + t1; // distance from root to u
		final double tv = tw + t2; // distance from root to v

		int[] logN = new int[N.length];

		double[][] pwd = new double[N.length][NUM_SAMPLES];
		double[][] wd = new double[N.length][NUM_SAMPLES];

		double mean1 = treeSimul.getM() * Math.exp(treeSimul.getLambda() * tu - treeSimul.getMu() * tu);
		double mean2 = treeSimul.getM() * Math.exp(treeSimul.getLambda() * tv - treeSimul.getMu() * tv);

		// run trials
		for (int trial = 0; trial < N.length; trial++) {
			int num = N[trial];

			logN[trial] = (int) Math.round(Math.log10(num));

			double[] pwdSample = new double[NUM_SAMPLES];
			double[] wdSample = new double[NUM_SAMPLES];

			for (int sample = 0; sample < NUM_SAMPLES; sample++) {
				double covariance = 0;
				for (int i = 0; i < N[trial];) {
					double product = treeSimul.covarianceComponent(tw, t1, t2, mean1, mean2);
					covariance += (product - covariance) / ++i;
				}

				InvertPairwiseDistance inverted = new InvertPairwiseDistance(treeSimul.getLambda(), treeSimul.getMu(),
						tu, tv, treeSimul.getM(),
						covariance);

				pwdSample[sample] = inverted.getPairwiseDistance();
				wdSample[sample] = inverted.getAncestorDistance();
			}

			pwd[trial] = trim(pwdSample);
			wd[trial] = trim(wdSample);
		}

		String xAxis = "log N";
		// display charts
		Chart.BoxWhiskerChart(xAxis, "pwd", pwd, (t1 + t2) * treeSimul.getMu(), N, false);
		Chart.BoxWhiskerChart(xAxis, "wd", wd, tw * treeSimul.getMu(), N, false);
	}

	/**
	 * Returns a copy of the given data array with any NaN elements deleted.
	 * 
	 * @param data the array to be trimmed
	 * @return a copy of data with NaN elements deleted.
	 */
	private static double[] trim(double[] data) {
		Arrays.sort(data);

		int index = data.length;
		while (index > 0 && Double.isNaN(data[--index]))
			;

		if (Double.isNaN(data[index]))
			return new double[] {};

		return Arrays.copyOfRange(data, 0, index + 1);
	}

	public static void main(String[] args) {
		// set TKF91 process parameters
		final double LAMBDA = 1; // insertion rate
		final double MU = 0.9; // deletion rate
		final double NU = 0.2; // substitution rate
		final double PI0 = 0.5; // probability a character is a 0 after substitution or insertion
		final String ROOT = "01010101"; // the root sequence

		// create TreeSimul object using TKF91 process parameters
		TreeSimul treeSimul = new TreeSimul(LAMBDA, MU, NU, PI0, ROOT);

		// values of N to use for the simulation
		// the 10^5 and 10^6 values will take some time to run
		// if you just want to quickly see that the code works, try only N = {10^2, 10^3, 10^4}
		int[] N = { (int) Math.pow(10, 3), (int) Math.pow(10, 4), (int) Math.pow(10,
				5), (int) Math.pow(10, 6) };

		int NUM_SAMPLES = 50; // number of trials for each N

		starTreeLength(treeSimul, N, NUM_SAMPLES);
	}
}
