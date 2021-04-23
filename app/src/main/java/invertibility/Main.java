package invertibility;

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
			TreeLeaves exampleTree = example.toTreeLeaves();
			Invert exampleInverted = new Invert(exampleTree);

			N[trial] = num;
			g[trial] = exampleInverted.getGamma();
			b[trial] = exampleInverted.getBeta();
			m[trial] = exampleInverted.getM();
		}

		// display charts
		Chart.chart("Gamma vs. N", "gamma", g, gamma, N);
		Chart.chart("Beta vs. N", "beta", b, beta, N);
		Chart.chart("M vs. N", "M", m, M, N);
	}

	public static void main(String[] args) {
		example1();
		example2();
	}
}
