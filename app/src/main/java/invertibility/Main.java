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
		double lambda = 1;
		double mu = 0.7;
		double nu = 0.4;
		double pi0 = 0.3;
		String root = "01010101";
		int M = root.length();

		int a = 0;
		for (char c: root.toCharArray()) {
			if (c == '1')
				a++;
		}

		int[] Ns = { 1000, 10000, 100000, 200000, 300000 };
		int numIter = Ns.length;

		// initialize arrays
		// logN is the number of leaves for each trial
		int[] logN = new int[numIter];

		double[] n = new double[numIter];
		double[] pi0Upper = new double[numIter];
		double[] pi0Lower = new double[numIter];
		double[] aUpper = new double[numIter];
		double[] aLower = new double[numIter];

		double[] nDeviation = new double[numIter];
		double[] pi0UpperDeviation = new double[numIter];
		double[] pi0LowerDeviation = new double[numIter];
		double[] aUpperDeviation = new double[numIter];
		double[] aLowerDeviation = new double[numIter];

		// run trials
		for (int trial = 0; trial < numIter; trial++) {
			int num = Ns[trial];
			System.out.println("Current trial: " + trial);

			logN[trial] = (int) Math.round(Math.log10(num));
			final int NUM_SAMPLES = 50;

			double[] nSample = new double[NUM_SAMPLES];
			double[] pi0UpperSample = new double[NUM_SAMPLES];
			double[] pi0LowerSample = new double[NUM_SAMPLES];
			double[] aUpperSample = new double[NUM_SAMPLES];
			double[] aLowerSample = new double[NUM_SAMPLES];

			for (int sample = 0; sample < NUM_SAMPLES; sample++) {
				TreeSimul example = new TreeSimul(lambda, mu, nu, pi0, root, Ns[trial]);
				TreeLeaves exampleTree = example.toTreeLeaves();

				Invert exampleInverted = new Invert(mu, lambda, M, exampleTree);

				nSample[sample] = exampleInverted.getNu();
				pi0UpperSample[sample] = exampleInverted.getPi0Upper();
				pi0LowerSample[sample] = exampleInverted.getPi0Lower();
				aUpperSample[sample] = exampleInverted.getAUpper();
				aLowerSample[sample] = exampleInverted.getALower();
			}

			System.out.println("TRIMMING");
			nSample = trim(nSample);
			pi0UpperSample = trim(pi0UpperSample);
			pi0LowerSample = trim(pi0LowerSample);
			aUpperSample = trim(aUpperSample);
			aLowerSample = trim(aLowerSample);

			n[trial] = average(nSample);
			nDeviation[trial] = standardDeviation(n[trial], nSample);
			
			pi0Upper[trial] = average(pi0UpperSample);
			pi0UpperDeviation[trial]  = standardDeviation(pi0Upper[trial], pi0UpperSample);

			pi0Lower[trial] = average(pi0LowerSample);
			pi0LowerDeviation[trial] = standardDeviation(pi0Lower[trial], pi0LowerSample);

			aUpper[trial] = average(aUpperSample);
			aUpperDeviation[trial]  = standardDeviation(aUpper[trial], aUpperSample);

			aLower[trial] = average(aLowerSample);
			aLowerDeviation[trial] = standardDeviation(aLower[trial], aLowerSample);

			System.out.println("Trial: " + trial);
			System.out.println("nu: " + n[trial] + " " + nDeviation[trial]);
			System.out.println("pi0Upper: " + pi0Upper[trial] + " " + pi0UpperDeviation[trial]);
			System.out.println("pi0Lower: " + pi0Lower[trial] + " " + pi0LowerDeviation[trial]);
			System.out.println("aUpper: " + aUpper[trial] + " " + aUpperDeviation[trial]);
			System.out.println("aLower: " + aLower[trial] + " " + aLowerDeviation[trial]);
		}

		double exp = Math.exp(-nu);
		String xAxis = "log N";
		// display charts
		System.out.println("exp: " + exp + " pi0: " + pi0 + " a: " + a);
		Chart.ErrorChart(xAxis, "exp", n, nDeviation, exp, logN);
		Chart.ErrorChart(xAxis, "pi0Upper", pi0Upper, pi0UpperDeviation, pi0, logN);
		Chart.ErrorChart(xAxis, "pi0Lower", pi0Lower, pi0LowerDeviation, pi0, logN);
		Chart.ErrorChart(xAxis, "aUpper", aUpper, aUpperDeviation, a, logN);
		Chart.ErrorChart(xAxis, "aLower", aLower, aLowerDeviation, a, logN);
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
		starTreeLength();
	}
}
