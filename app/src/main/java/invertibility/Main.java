package invertibility;

import java.io.PrintWriter;
import java.util.Random;

public class Main {
	/**
	 * Tests the InvertLength class. Simulates the length process for multiple
	 * trials for each value of N and plots the estimated gamma, beta, and M values.
	 * 
	 * @param treeSimul   the TreeSimul object containing the TKF91 process
	 *                    parameters
	 * @param N           an array containing the values of N to simulate
	 * @param NUM_SAMPLES the number of trials to use for each N
	 * @param tag         the tag for the output files
	 */
	public static void testInvertLength(TreeSimul treeSimul, int[] N,
			int NUM_SAMPLES, String tag) {
		double GAMMA = treeSimul.getLambda() / treeSimul.getMu();
		double BETA = Math.exp(treeSimul.getLambda() - treeSimul.getMu());

		double[][] gamma = new double[N.length][NUM_SAMPLES];
		double[][] beta = new double[N.length][NUM_SAMPLES];
		double[][] m = new double[N.length][NUM_SAMPLES];

		// run trials
		for (int trial = 0; trial < N.length; trial++) {
			long t = System.nanoTime();
			double[] gammaSample = new double[NUM_SAMPLES];
			double[] betaSample = new double[NUM_SAMPLES];
			double[] mSample = new double[NUM_SAMPLES];

			for (int sample = 0; sample < NUM_SAMPLES; sample++) {
				treeSimul.runTKF91Process(N[trial]);
				LeafSamples samples = treeSimul.toLeafSamples();

				InvertLength inverted = new InvertLength(samples);

				gammaSample[sample] = inverted.getGamma();
				betaSample[sample] = inverted.getBeta();

				Integer sampleM = inverted.getM();
				if (sampleM == null)
					mSample[sample] = Double.NaN;
				else
					mSample[sample] = sampleM;
				System.out.println("sampleM: " + sampleM);
			}

			gamma[trial] = gammaSample;
			beta[trial] = betaSample;
			m[trial] = mSample;

			System.out.println((System.nanoTime() - t) / 1e9 + " seconds has passed, trial " + trial + ".");
			t = System.nanoTime();
		}

		saveData(gamma, GAMMA, "gamma_" + tag + ".csv");
		saveData(beta, BETA, "beta_" + tag + ".csv");
		saveData(m, treeSimul.getM(), "m_" + tag + ".csv");

		// String xaxis = "N";
		// display charts
		// Chart.BoxWhiskerChart(xAxis, "gamma", gamma, GAMMA, N, false);
		// Chart.BoxWhiskerChart(xAxis, "beta", beta, BETA, N, false);
		// Chart.BoxWhiskerChart(xAxis, "M", m, treeSimul.getM(), N, false);
	}

	/**
	 * Tests the Invert1Mer class. Simulates the 1mer process for multiple
	 * trials for each value of N and plots the estimated nu and a values.
	 * 
	 * @param treeSimul   the TreeSimul object containing the TKF91 process
	 *                    parameters
	 * @param N           an array containing the values of N to simulate
	 * @param NUM_SAMPLES the number of trials to use for each N
	 * @param tag         the tag for the output files
	 */
	public static void testInvert1Mer(TreeSimul treeSimul, int[] N,
			int NUM_SAMPLES, String tag) {
		int numOnes = 0;
		for (char c : treeSimul.getRoot().toCharArray()) {
			if (c == '1')
				numOnes++;
		}
		final int A = numOnes;

		double[][] nu = new double[N.length][NUM_SAMPLES];
		double[][] a = new double[N.length][NUM_SAMPLES];

		// run trials
		for (int trial = 0; trial < N.length; trial++) {
			long t = System.nanoTime();
			double[] nuSample = new double[NUM_SAMPLES];
			double[] aSample = new double[NUM_SAMPLES];

			for (int sample = 0; sample < NUM_SAMPLES; sample++) {
				treeSimul.runTKF91Process(N[trial]);
				LeafSamples samples = treeSimul.toLeafSamples();

				Invert1Mer inverted = new Invert1Mer(treeSimul.getLambda(), treeSimul.getMu(), treeSimul.getPi0(),
						treeSimul.getM(), samples);

				nuSample[sample] = inverted.getNu();
				aSample[sample] = inverted.getA();
			}

			nu[trial] = nuSample;
			a[trial] = aSample;
			System.out.println((System.nanoTime() - t) / 1e9 + " seconds has passed, trial " + trial + ".");
			t = System.nanoTime();
		}

		saveData(nu, treeSimul.getNu(), "nu_" + tag + ".csv");
		saveData(a, A, "a_" + tag + ".csv");

		// String xaxis = "N";
		// display charts
		// Chart.BoxWhiskerChart(xAxis, "nu", nu, treeSimul.getNu(), N, false);
		// Chart.BoxWhiskerChart(xAxis, "a", a, A, N, false);
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
	 * @param tag         the tag for the output files
	 */
	public static void testInvertState(TreeSimul treeSimul, int[] N, int NUM_SAMPLES, String tag) {
		String[] root = new String[N.length];

		double[][] diff = new double[N.length][NUM_SAMPLES];

		// run trials
		for (int trial = 0; trial < N.length; trial++) {
			long t = System.nanoTime();
			double[] diffSample = new double[NUM_SAMPLES];

			int[] count = new int[treeSimul.getM()]; // counts the number of 1's across all samples at each index

			for (int sample = 0; sample < NUM_SAMPLES; sample++) {
				treeSimul.runTKF91Process(N[trial]);
				LeafSamples samples = treeSimul.toLeafSamples();

				InvertState inverted = new InvertState(treeSimul.getLambda(), treeSimul.getMu(), treeSimul.getNu(),
						treeSimul.getM(), treeSimul.getPi0(), samples);

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
			// print out time taken and diff for this sample
			System.out.println((System.nanoTime() - t) / 1e9 + " seconds has passed, trial " + trial + ".");
			System.out.println("estimated root: " + root[trial]);
			System.out.println("actual root:    " + treeSimul.getRoot());
			t = System.nanoTime();
		}

		saveData(diff, -1, "diff_" + tag + ".csv");

		// String xaxis = "N";
		// display charts
		// Chart.BoxWhiskerChart(xAxis, "difference", diff, -1, N, true);
	}

	/**
	 * Tests the InvertLength, Invert1Mer, and InvertState classes in sequence,
	 * using invertibility results from one class into the next.
	 * 
	 * @param treeSimul   the TreeSimul object containing the TKF91 process
	 *                    parameters
	 * @param N           an array containing the values of N to simulate
	 * @param NUM_SAMPLES the number of trials to use for each N
	 * @param tag         the tag for the output files
	 */
	public static void testInvertLength1MerState(TreeSimul treeSimul, int[] N, int NUM_SAMPLES, String tag) {
		int numOnes = 0;
		for (char c : treeSimul.getRoot().toCharArray()) {
			if (c == '1')
				numOnes++;
		}
		final int A = numOnes;

		double[][] gamma = new double[N.length][NUM_SAMPLES];
		double[][] beta = new double[N.length][NUM_SAMPLES];
		double[][] m = new double[N.length][NUM_SAMPLES];

		double[][] nu = new double[N.length][NUM_SAMPLES];
		double[][] a = new double[N.length][NUM_SAMPLES];

		double[][] diff = new double[N.length][NUM_SAMPLES];

		for (int trial = 0; trial < N.length; trial++) {
			long t = System.nanoTime();
			double[] gammaSample = new double[NUM_SAMPLES];
			double[] betaSample = new double[NUM_SAMPLES];
			double[] mSample = new double[NUM_SAMPLES];

			double[] nuSample = new double[NUM_SAMPLES];
			double[] aSample = new double[NUM_SAMPLES];

			String[] rootSample = new String[NUM_SAMPLES];
			double[] diffSample = new double[NUM_SAMPLES];

			for (int sample = 0; sample < NUM_SAMPLES; sample++) {
				treeSimul.runTKF91Process(N[trial]);
				LeafSamples samples = treeSimul.toLeafSamples();

				InvertLength invertedLength = new InvertLength(samples);
				gammaSample[sample] = invertedLength.getGamma();
				betaSample[sample] = invertedLength.getBeta();
				double lambda = invertedLength.getLambda();
				double mu = invertedLength.getMu();
				Integer sampleM = invertedLength.getM();
				if (sampleM == null) {
					mSample[sample] = Double.NaN;
					nuSample[sample] = Double.NaN;
					aSample[sample] = Double.NaN;
					diffSample[sample] = Double.NaN;
				} else {
					mSample[sample] = sampleM;
					Invert1Mer inverted1Mer = new Invert1Mer(lambda, mu, treeSimul.getPi0(),
							sampleM, samples);
					nuSample[sample] = inverted1Mer.getNu();
					aSample[sample] = inverted1Mer.getA();

					InvertState invertedState = new InvertState(lambda, mu, nuSample[sample],
							Math.min(sampleM, treeSimul.getM()), treeSimul.getPi0(), samples);
					rootSample[sample] = invertedState.getRootState();
					for (int i = 0; i < sampleM && i < treeSimul.getM(); i++) {
						if (rootSample[sample].charAt(i) != treeSimul.getRoot().charAt(i))
							diffSample[sample]++;
					}
					// add difference in lengths to diffSample
					diffSample[sample] += Math.abs(sampleM - treeSimul.getM());
				}

			}

			gamma[trial] = gammaSample;
			beta[trial] = betaSample;
			m[trial] = mSample;
			nu[trial] = nuSample;
			a[trial] = aSample;
			diff[trial] = diffSample;
			System.out.println((System.nanoTime() - t) / 1e9 + " seconds has passed, trial " + trial + ".");
			t = System.nanoTime();
		}

		saveData(gamma, treeSimul.getLambda() / treeSimul.getMu(), "gamma_fullseq_" + tag + ".csv");
		saveData(beta, Math.exp(treeSimul.getLambda() - treeSimul.getMu()), "beta_fullseq_" + tag + ".csv");
		saveData(m, treeSimul.getM(), "m_fullseq_" + tag + ".csv");
		saveData(nu, treeSimul.getNu(), "nu_fullseq_" + tag + ".csv");
		saveData(a, A, "a_fullseq_" + tag + ".csv");
		saveData(diff, -1, "diff_fullseq_" + tag + ".csv");

		// String xaxis = "N";
		// display charts
		// Chart.BoxWhiskerChart(xAxis, "gamma", gamma, treeSimul.getLambda() /
		// treeSimul.getMu(), N, false);
		// Chart.BoxWhiskerChart(xAxis, "beta", beta, Math.exp(treeSimul.getLambda() -
		// treeSimul.getMu()), N, false);
		// Chart.BoxWhiskerChart(xAxis, "M", m, treeSimul.getM(), N, false);
		// Chart.BoxWhiskerChart(xAxis, "nu", nu, treeSimul.getNu(), N, false);
		// Chart.BoxWhiskerChart(xAxis, "a", a, A, N, false);
		// Chart.BoxWhiskerChart(xAxis, "difference", diff, -1, N, true);
	}

	/**
	 * Tests the InvertLength and Invert1Mer classes in sequence,
	 * using invertibility results from InvertLength into the next.
	 * 
	 * @param treeSimul   the TreeSimul object containing the TKF91 process
	 *                    parameters
	 * @param N           an array containing the values of N to simulate
	 * @param NUM_SAMPLES the number of trials to use for each N
	 * @param tag         the tag for the output files
	 */
	public static void testInvertLength1Mer(TreeSimul treeSimul, int[] N, int NUM_SAMPLES, String tag) {
		int numOnes = 0;
		for (char c : treeSimul.getRoot().toCharArray()) {
			if (c == '1')
				numOnes++;
		}
		final int A = numOnes;

		double[][] gamma = new double[N.length][NUM_SAMPLES];
		double[][] beta = new double[N.length][NUM_SAMPLES];
		double[][] m = new double[N.length][NUM_SAMPLES];

		double[][] nu = new double[N.length][NUM_SAMPLES];
		double[][] a = new double[N.length][NUM_SAMPLES];

		for (int trial = 0; trial < N.length; trial++) {
			long t = System.nanoTime();
			double[] gammaSample = new double[NUM_SAMPLES];
			double[] betaSample = new double[NUM_SAMPLES];
			double[] mSample = new double[NUM_SAMPLES];

			double[] nuSample = new double[NUM_SAMPLES];
			double[] aSample = new double[NUM_SAMPLES];

			for (int sample = 0; sample < NUM_SAMPLES; sample++) {
				treeSimul.runTKF91Process(N[trial]);
				LeafSamples samples = treeSimul.toLeafSamples();

				InvertLength invertedLength = new InvertLength(samples);
				gammaSample[sample] = invertedLength.getGamma();
				betaSample[sample] = invertedLength.getBeta();
				double lambda = invertedLength.getLambda();
				double mu = invertedLength.getMu();
				Integer sampleM = invertedLength.getM();
				if (sampleM == null) {
					mSample[sample] = Double.NaN;
					nuSample[sample] = Double.NaN;
					aSample[sample] = Double.NaN;
				} else {
					mSample[sample] = sampleM;
					Invert1Mer inverted1Mer = new Invert1Mer(lambda, mu, treeSimul.getPi0(),
							sampleM, samples);
					nuSample[sample] = inverted1Mer.getNu();
					aSample[sample] = inverted1Mer.getA();
				}
			}

			gamma[trial] = gammaSample;
			beta[trial] = betaSample;
			m[trial] = mSample;
			nu[trial] = nuSample;
			a[trial] = aSample;
			System.out.println((System.nanoTime() - t) / 1e9 + " seconds has passed, trial " + trial + ".");
			t = System.nanoTime();
		}

		saveData(gamma, treeSimul.getLambda() / treeSimul.getMu(), "gamma_shortseq_" + tag + ".csv");
		saveData(beta, Math.exp(treeSimul.getLambda() - treeSimul.getMu()), "beta_shortseq_" + tag + ".csv");
		saveData(m, treeSimul.getM(), "m_shortseq_" + tag + ".csv");
		saveData(nu, treeSimul.getNu(), "nu_shortseq_" + tag + ".csv");
		saveData(a, A, "a_shortseq_" + tag + ".csv");

		// String xaxis = "N";
		// display charts
		// Chart.BoxWhiskerChart(xAxis, "gamma", gamma, treeSimul.getLambda() /
		// treeSimul.getMu(), N, false);
		// Chart.BoxWhiskerChart(xAxis, "beta", beta, Math.exp(treeSimul.getLambda() -
		// treeSimul.getMu()), N, false);
		// Chart.BoxWhiskerChart(xAxis, "M", m, treeSimul.getM(), N, false);
		// Chart.BoxWhiskerChart(xAxis, "nu", nu, treeSimul.getNu(), N, false);
		// Chart.BoxWhiskerChart(xAxis, "a", a, A, N, false);
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
	 * @param tag         the tag for the output files
	 */
	public static void testInvertPairwiseDistance(TreeSimul treeSimul, int[] N, int NUM_SAMPLES, String tag) {
		final double tw = 1; // distance from root to w
		final double t1 = 2; // distance from w to u
		final double t2 = 3; // distance from w to v
		final double tu = tw + t1; // distance from root to u
		final double tv = tw + t2; // distance from root to v

		double[][] pwd = new double[N.length][NUM_SAMPLES];
		double[][] wd = new double[N.length][NUM_SAMPLES];

		double mean1 = treeSimul.getM() * Math.exp(treeSimul.getLambda() * tu - treeSimul.getMu() * tu);
		double mean2 = treeSimul.getM() * Math.exp(treeSimul.getLambda() * tv - treeSimul.getMu() * tv);

		// run trials
		for (int trial = 0; trial < N.length; trial++) {
			long t = System.nanoTime();

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

			pwd[trial] = pwdSample;
			wd[trial] = wdSample;
			System.out.println((System.nanoTime() - t) / 1e9 + " seconds has passed, trial " + trial + ".");
			t = System.nanoTime();
		}

		saveData(pwd, (t1 + t2) * treeSimul.getMu(), "pwd_" + tag + ".csv");
		saveData(wd, tw * treeSimul.getMu(), "wd_" + tag + ".csv");

		// String xaxis = "N";
		// display charts
		// Chart.BoxWhiskerChart(xAxis, "pwd", pwd, (t1 + t2) * treeSimul.getMu(), N,
		// false);
		// Chart.BoxWhiskerChart(xAxis, "wd", wd, tw * treeSimul.getMu(), N, false);
	}

	/**
	 * Saves the simulation data to a CSV file.
	 * 
	 * @param data     the simulation data to save
	 * @param exact    the exact value to include in the first column of every row
	 * @param filename the name of the output file to save in the data directory
	 */
	private static void saveData(double[][] data, double exact, String filename) {
		try (PrintWriter out = new PrintWriter("data/" + filename)) {
			for (double[] row : data) {
				out.print(exact);
				for (int j = 0; j < row.length; j++) {
					out.print(",");
					out.print(row[j]);
				}
				out.println();
			}

			out.flush();
			System.out.println("Data saved to data/" + filename);
		} catch (Exception e) {
			e.printStackTrace();
		}
	}

	public static void main(String[] args) {
		// set TKF91 process parameters
		boolean longM = false; // whether to use a long root sequence (250) or a short one (8)
		boolean swap = true; // whether to swap the insertion and deletion rates

		String tag = (longM) ? "longM" : "shortM"; // tag for output files
		if (swap) {
			tag += "_swap";
		}

		int M = (longM) ? 250 : 8; // desired length of root sequence

		double LAMBDA = .05; // insertion rate
		double MU = 0.075; // deletion rate
		if (swap) {
			double temp = LAMBDA;
			LAMBDA = MU;
			MU = temp;
		}
		double NU = 1; // substitution rate
		double PI0 = .5; // probability a character is a 0 after substitution or insertion
		
		Random rand = new Random(123); // fixed seed for reproducibility
		StringBuilder sb = new StringBuilder(M);

		for (int i = 0; i < M; i++) {
			sb.append(rand.nextInt(2)); // generates 0 or 1
		}

		String ROOT = sb.toString();
		// ROOT = "01010101"; // you can change the root sequence here if you want

		// create TreeSimul object using TKF91 process parameters
		TreeSimul treeSimul = new TreeSimul(LAMBDA, MU, NU, PI0, ROOT);

		// values of N to use for the simulation
		int[] N = (longM) ? new int[] { 25, 100, 400 } : new int[] { (int) 1e3, (int) 1e4, (int) 1e5, (int) 1e6};
		int NUM_SAMPLES = 50; // number of trials for each N

		// testInvertLength(treeSimul, N, NUM_SAMPLES, tag);
		// testInvert1Mer(treeSimul, N, NUM_SAMPLES, tag);
		// testInvertPairwiseDistance(treeSimul, N, NUM_SAMPLES, tag);
		// testInvertLength1Mer(treeSimul, N, NUM_SAMPLES, tag);

		testInvertState(treeSimul, N, NUM_SAMPLES, tag);
		testInvertLength1MerState(treeSimul, N, NUM_SAMPLES, tag);
	}
}
