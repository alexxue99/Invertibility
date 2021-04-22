package invertibility;

import java.util.LinkedList;

import invertibility.TreeSimul.Edge;

public class Main {
	/**
	 * Basic example to show how to input a tree, invert it, and print out the
	 * estimated parameters.
	 */
	public static void example1() {
		System.out.println("Starting example 1...\n");
		// set length process parameters
		double lambda = 0.7;
		double mu = 1;
		double nu = 3;
		double pi0 = 0.3;
		double time = 1;
		String root = "0110";
		int rootVertex = 0;

		// set number of leaves for tree

		int num = 1000;
		// set the actual values for M, gamma, and beta
		int M = root.length();
		double gamma = lambda / mu;
		double beta = Math.exp((lambda - mu) * time);

		// set the directed edges for the star tree
		LinkedList<Edge> edges = new LinkedList<Edge>();
		double[] t = new double[num];
		for (int i = 0; i < num; i++) {
			edges.add(new Edge(0, i + 1));
			t[i] = time;
		}

		TreeSimul example = new TreeSimul(lambda, mu, nu, pi0, rootVertex, root, edges, t);
		Invert exampleInverted = new Invert(example);

		// Print out estimated and actual parameters.
		// Some of the estimated parameters may show to be NaN because the calculations
		// may take the square root of a negative number.
		System.out.println("Estimated gamma is " + exampleInverted.getGamma() + ".");
		System.out.println("Actual gamma is " + gamma + ".");
		System.out.println("Estimated beta is " + exampleInverted.getBeta() + ".");
		System.out.println("Actual beta is " + beta + ".");
		System.out.println("Estimated M is " + exampleInverted.getM() + ".");
		System.out.println("Actual M is " + M + ".\n");
	}

	/** Example to show how to run trials on a variable number of leaves. */
	public static void example2() {
		System.out.println("Starting example 2...");
		// set length process parameters
		double lambda = 0.7;
		double mu = 1;
		double nu = 3;
		double pi0 = 0.3;
		double time = 1;
		String root = "0110";
		int rootVertex = 0;

		// set max number of leaves and step size for number of leaves in each trial
		int maxN = 5000;
		int stepSize = 100;
		int numIter = maxN / stepSize;

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

		// set the actual values for M, gamma, and beta
		int M = root.length();
		double gamma = lambda / mu;
		double beta = Math.exp((lambda - mu) * time);

		// run trials
		for (int num = stepSize, trial = 0; num <= maxN; num += stepSize, trial++) {
			// set the directed edges for the star tree
			LinkedList<Edge> edges = new LinkedList<Edge>();
			double[] t = new double[num];
			for (int i = 0; i < num; i++) {
				edges.add(new Edge(0, i + 1));
				t[i] = time;
			}

			TreeSimul example = new TreeSimul(lambda, mu, nu, pi0, rootVertex, root, edges, t);
			Invert exampleInverted = new Invert(example);

			N[trial] = num;
			g[trial] = exampleInverted.getGamma();
			b[trial] = exampleInverted.getBeta();
			m[trial] = exampleInverted.getM();
		}

		// display charts for the data
		Chart.chart("Gamma vs. N", "gamma", g, gamma, N);
		Chart.chart("Beta vs. N", "beta", b, beta, N);
		Chart.chart("M vs. N", "M", m, M, N);
	}

	public static void main(String[] args) {
		example1();
		example2();
	}
}
