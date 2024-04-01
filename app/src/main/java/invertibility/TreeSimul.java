package invertibility;

import java.util.Random;

/**
 * Class used to simulate the length process on a star tree with N leaves.
 */
public class TreeSimul {
	private double lambda;
	private double mu;
	private double nu;
	private double pi0;
	private String root;
	private int M;
	private int N;

	private String[] seqLeaves; // contains information about the sequences at the leaves

	private Random rand;

	/**
	 * Constructor to initialize tree.
	 * 
	 * @param lambda
	 *               is insertion rate.
	 * @param mu
	 *               is deletion rate.
	 * @param nu
	 *               is substitution rate.
	 * @param pi0
	 *               is probability that a character is 0.
	 * @param root
	 *               is the string at the root.
	 */

	public TreeSimul(double lambda, double mu, double nu, double pi0, String root,
			int N) {
		this.lambda = lambda;
		this.mu = mu;
		this.nu = nu;
		this.pi0 = pi0;
		this.root = root;
		M = root.length();
		this.N = N;

		seqLeaves = new String[N];

		rand = new Random();

		runLengthProcess();
	}

	/** Returns a variable with exponential distribution with parameter var. */
	private double getExp(double var) {
		Random rand = new Random();
		return Math.log(1 - rand.nextDouble()) / (-var);
	}

	/** Returns a 0 with probability pi0 and a 1 with probability 1 - pi0. */
	private char getChar() {
		return (rand.nextDouble() <= pi0) ? '0' : '1';
	}

	/**
	 * Evolves string s along an edge, using the length process
	 * parameters.
	 */
	private String evolve(String s, double t) {
		while (t > 0) {
			double[] timings = new double[s.length() * 3];

			// get exponential variables under parameters lambda, mu, and nu, and then find
			// the smallest variable
			for (int i = 0; i < s.length(); i++) {
				timings[3 * i] = getExp(lambda);
				timings[3 * i + 1] = getExp(mu);
				timings[3 * i + 2] = getExp(nu);
			}

			double min = Double.MAX_VALUE;
			int minIndex = 0;

			for (int i = 0; i < timings.length; i++) {
				if (timings[i] < min) {
					minIndex = i;
					min = timings[i];
				}
			}

			t -= min;
			if (t < 0)
				break;

			int location = minIndex / 3;
			switch (minIndex % 3) {
				case 0: // insertion
					s = s.substring(0, location + 1) + getChar() + s.substring(location + 1);
					break;
				case 1: // deletion
					s = s.substring(0, location) + s.substring(location + 1);
					break;
				case 2: // substitution
					s = s.substring(0, location) + getChar() + s.substring(location + 1);
					break;
			}
		}

		return s;
	}

	public double approxLength(double t, int num) {
		int sum = 0;

		for (int i = 0; i < num; i++) {
			sum += evolve(root, t).length();
		}

		return (double) sum / num;
	}

	public double pairwiseDistance(double tw, double t1, double t2, double l1_approx, double l2_approx, double mean1, double mean2) {
		String s = evolve(root, tw);

		int l1 = evolve(s, t1).length();
		int l2 = evolve(s, t2).length();

		return (l1 - l1_approx - mean1) * (l2 - l2_approx - mean2);
	}

	/** Simulates the length process on the whole tree. */
	public void runLengthProcess() {
		for (int i = 0; i < N; i++) {
			seqLeaves[i] = evolve(root, 1);
		}
	}

	/**
	 * Transforms a TreeSimul object into a TreeLeaves object.
	 */
	public TreeLeaves toTreeLeaves() {
		return new TreeLeaves(seqLeaves);
	}

	public double getLambda() {
		return lambda;
	}

	public double getMu() {
		return mu;
	}

	public double getNu() {
		return nu;
	}

	public double getPi0() {
		return pi0;
	}

	public String getRoot() {
		return root;
	}

	public int getM() {
		return M;
	}

	public String[] getSeqLeaves() {
		return seqLeaves;
	}
}