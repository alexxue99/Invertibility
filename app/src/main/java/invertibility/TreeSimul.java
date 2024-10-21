package invertibility;

import java.util.Random;

/**
 * Class used to simulate the TKF91 process for a leaf of the phylogeny.
 */
public class TreeSimul {
	private double lambda;
	private double mu;
	private double nu;
	private double pi0;
	private String root;
	private int M;

	private String[] samples; // contains information about the sampled sequences

	private Random rand;

	/**
	 * Constructor to initialize a TreeSimul object.
	 * 
	 * @param lambda the insertion rate
	 * @param mu     the deletion rate
	 * @param nu     the substitution rate
	 * @param pi0    the probability that a character is 0
	 * @param root   the string at the root
	 */

	public TreeSimul(double lambda, double mu, double nu, double pi0, String root) {
		this.lambda = lambda;
		this.mu = mu;
		this.nu = nu;
		this.pi0 = pi0;
		this.root = root;
		M = root.length();
		rand = new Random();
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
	 * Evolves string s along an edge for time t, using the TKF91 process
	 * parameters.
	 * 
	 * @param s the string to evolve
	 * @param t the time
	 * @return the evolved string
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

	/**
	 * Approximates the average length of the root sequence after time t by taking
	 * an average over num trials.
	 * 
	 * @param t   the time
	 * @param num the number of trials
	 * @return the approximate length
	 */
	public double approxLength(double t, int num) {
		int sum = 0;

		for (int i = 0; i < num; i++) {
			sum += evolve(root, t).length();
		}

		return (double) sum / num;
	}

	/**
	 * Runs the tree process for a forked tree. The root node is connected a to a
	 * node w with length tw, and the node w is connected to two nodes u and v with
	 * lengths t1 and t2.
	 * Returns (Lu - mu)(Lv - mv).
	 * 
	 * @param tw    the distance from root to w
	 * @param t1    the distance from w to u
	 * @param t2    the distance from w to v
	 * @param mean1 the average length mu of the sequence at u
	 * @param mean2 the average length mv of the sequence at v
	 * @return the covariance component for this instance of a forked tree
	 */
	public double covarianceComponent(double tw, double t1, double t2, double mean1,
			double mean2) {
		String s = evolve(root, tw);

		int l1 = evolve(s, t1).length();
		int l2 = evolve(s, t2).length();

		return (l1 - mean1) * (l2 - mean2);
	}

	/**
	 * Resets the samples and simulates the TKF91 process on the whole
	 * tree.
	 */
	public void runTKF91Process(int n) {
		samples = new String[n];

		for (int i = 0; i < n; i++) {
			samples[i] = evolve(root, 1);
		}
	}

	/**
	 * Transforms a TreeSimul object into a LeafSamples object.
	 */
	public LeafSamples toLeafSamples() {
		return new LeafSamples(samples);
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

	public String[] getSamples() {
		return samples;
	}
}