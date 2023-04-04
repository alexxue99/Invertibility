package invertibility;

import java.util.Iterator;
import java.util.List;

/**
 * Class used to invert the length process parameters, given the lengths of the
 * sequences at the leaves of a tree.
 */
public class Invert {
	final static int DIVISION_PRECISION = 8; // number of decimal digits in divisions
	private TreeLeaves tree;

	private double mu;
	private double lambda;
	private double gamma;
	private double beta;
	private double M;
	private double nu;
	private double a;
	private double pi0;

	private double[] partial;
	private double[] C;

	private double C2prime;
	private double C3prime;
	private double[] D;

	/**
	 * Inverts the process based on given data on the lengths of the sequences at
	 * the leaves to estimate gamma, beta, and M.
	 */
	public Invert(TreeLeaves tree) {
		this.tree = tree;
		partial = new double[3];
		C = new double[3];

		calcPartials("length");
		updateCs();

		estimateGamma();
		estimateBeta();
		estimateM();
	}

	/**
	 * Inverts the process based on given data on the lengths of the sequences at
	 * the leaves, as well as on mu, lambda, and M, to estimate a, nu, and pi0.
	 */
	public Invert(double mu, double lambda, double M, TreeLeaves tree) {
		this.tree = tree;
		partial = new double[3];
		C = new double[3];
		D = new double[3];

		this.mu = mu;
		this.lambda = lambda;
		this.M = M;

		calcPartials("1-mer");
		updateCs();

		updateDs();
	}

	private void calcPartials(String type) {
		List<Integer> seqLeavesLengths = null;

		switch(type) {
			case "length": seqLeavesLengths = tree.getSeqLeavesLengths();
			break;
			case "1-mer": seqLeavesLengths = tree.get
		}
		 
		calcPartials(seqLeavesLengths);
	}

	private void calcPartials(List<Integer> seqLeavesLengths) {
		Iterator<Integer> leafIter = seqLeavesLengths.iterator();

		// estimate partials of G(z, t) with respect to z and evaluated at z = 1 and
		// time t,
		// by using the expected value of the kth factorial moment of the length process

		int num = 1;
		while (leafIter.hasNext()) {
			long length = (long) leafIter.next();
			long[] lengths = new long[] { length, length * (length - 1), length * (length - 1) * (length - 2) };

			for (int i = 0; i < 3; i++) {
				try {
					partial[i] += (lengths[i] - partial[i]) / num;
				} catch (ArithmeticException e) {
					System.out.println("PARTIALS");
				}
			}

			num++;
		}
	}

	/* Prereq: partials are calculated */
	private void updateCs() {
		C[0] = partial[0];
		C[1] = partial[1] - partial[0] * partial[0];
		C[2] = partial[2] + 2 * Math.pow(partial[0], 3) - 3 * partial[0] * partial[1];
		if (C[0] == 0) {
			C2prime = Double.NaN;
			C3prime = Double.NaN;
			System.out.println("here");
		} else {
			C2prime = C[1] / C[0];
			C3prime = C[2] / C[0];

		}
	}

	/* Prereq: Cs are calculated */
	private void updateDs() {
		D[0] = C[0] * Math.exp(mu);
	}

	/** Estimates gamma. */
	private void estimateGamma() {
		// System.out.println("gamma: " + -(C2prime + 1) * (C2prime + 1) * (3 * C2prime
		// * C2prime - 2 * C3prime));
		gamma = Math.sqrt(-(C2prime + 1) * (C2prime + 1) * (3 * C2prime * C2prime - 2 * C3prime));

		// if (Double.isNaN(gamma)) {
		// System.out.println(-(C2prime + 1) * (C2prime + 1) * (3 * C2prime * C2prime -
		// 2 * C3prime));
		// } else {
		// System.out.println(gamma);
		// }

		gamma += -C2prime * C2prime + C2prime + C3prime;

		try {
			gamma /= 2 * C2prime * C2prime + 2 * C2prime - C3prime + 2;
		} catch (ArithmeticException e) {
			System.out.println("here2");
			gamma = Double.NaN;
		}
	}

	/** Estimates beta. */
	private void estimateBeta() {
		beta = gamma * (2 + C2prime) - C2prime;

		try {
			beta /= 1 + gamma;
		} catch (ArithmeticException e) {
			System.out.println("here3");
			beta = Double.NaN;
		}
	}

	/** Estimates M. */
	private void estimateM() {
		if (Double.isNaN(beta) || beta == 0) {
			// System.out.println("here4");
			M = Double.NaN;
		} else {
			// System.out.println("M calculation: " + C[0] + " " + beta);
			M = Math.round(C[0] / beta);
		}
	}

	public double getGamma() {
		return gamma;
	}

	public double getBeta() {
		return beta;
	}

	public double getM() {
		return M;
	}

	public double[] getCs() {
		return C;
	}
}
