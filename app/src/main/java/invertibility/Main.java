package invertibility;

import java.util.Arrays;

public class Main {
	public static void starTreeLength() {
		// set length process parameters
		double lambda = 1;
		double mu = 0.7;
		double nu = 3;
		double pi0 = 0.3;
		double time = 1;
		String root = "01010101";

		// int[] Ns = { 1000, 10000, 100000, 1000000 };
		int[] Ns = { 1000, 10000, 100000 };
		int numIter = Ns.length;

		// initialize arrays
		// logN is the number of leaves for each trial
		// g is the estimated gamma for each trial
		// b is the estimated beta for each trial
		// m is the estimated m for each trial
		int[] logN = new int[numIter];

		double[] g = new double[numIter];
		double[] b = new double[numIter];
		double[] m = new double[numIter];

		double[] gDeviation = new double[numIter];
		double[] bDeviation = new double[numIter];
		double[] mDeviation = new double[numIter];

		// set the actual values for M, gamma, beta
		int M = root.length();
		double gamma = lambda / mu;
		double beta = Math.exp((lambda - mu) * time);

		// run trials
		for (int trial = 0; trial < numIter; trial++) {
			int num = Ns[trial];
			System.out.println("Current trial: " + trial);

			logN[trial] = (int) Math.round(Math.log10(num));
			final int NUM_SAMPLES = 50;

			double[] gSample = new double[NUM_SAMPLES];
			double[] bSample = new double[NUM_SAMPLES];
			double[] mSample = new double[NUM_SAMPLES];

			for (int sample = 0; sample < NUM_SAMPLES; sample++) {
				TreeSimul example = new TreeSimul(lambda, mu, nu, pi0, root, Ns[trial]);
				TreeLeaves exampleTree = example.toTreeLeaves();

				Invert exampleInverted = new Invert(exampleTree);

				gSample[sample] = exampleInverted.getGamma();
				bSample[sample] = exampleInverted.getBeta();
				mSample[sample] = exampleInverted.getM();
			}

			gSample = trim(gSample);
			bSample = trim(bSample);
			mSample = trim(mSample);

			g[trial] = average(gSample);
			gDeviation[trial] = standardDeviation(g[trial], gSample);

			b[trial] = average(bSample);
			bDeviation[trial] = standardDeviation(b[trial], bSample);

			m[trial] = average(mSample);
			mDeviation[trial] = standardDeviation(m[trial], mSample);

			System.out.println("Trial: " + trial);
			System.out.println("gamma: " + g[trial] + " " + gDeviation[trial]);
			System.out.println("beta: " + b[trial] + " " + bDeviation[trial]);
			System.out.println("M: " + m[trial] + " " + mDeviation[trial]);
		}

		String xAxis = "log N";
		// display charts
		System.out.println("gamma: " + gamma + " beta: " + beta + " M: " + M);
		Chart.ErrorChart(xAxis, "gamma", g, gDeviation, gamma, logN);
		Chart.ErrorChart(xAxis, "beta", b, bDeviation, beta, logN);
		Chart.ErrorChart(xAxis, "M", m, mDeviation, M, logN);
	}

	public static void starTree1Mer() {
		// set length process parameters
		final double LAMBDA = 1;
		final double MU = 0.7;
		final double NU = 0.6;
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
				TreeSimul example = new TreeSimul(LAMBDA, MU, NU, PI0, ROOT, N[trial]);
				TreeLeaves exampleTree = example.toTreeLeaves();

				Invert exampleInverted = new Invert(MU, LAMBDA, M, exampleTree);

				nuSample[sample] = exampleInverted.getNu();
				pi0Sample[sample] = exampleInverted.getPi0();
				aSample[sample] = exampleInverted.getA();
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
		starTree1Mer();
	}
}
