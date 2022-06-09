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
	private double[] partial;
	private double[] C;
	private double C2prime;
	private double C3prime;
	private double gamma;
	private double beta;
	private double M;

	public Invert(TreeLeaves tree) {
		this.tree = tree;
		partial = new double[3];
		C = new double[3];
		invert();
	}

	/**
	 * Inverts the process based on given data on partials
	 */
	public Invert(double[] partial) {
		this.partial = partial;
		C = new double[3];
		updateCs();

		estimateGamma();
		estimateBeta();
		estimateM();
	}
	
	/**
	 * Calculates C<SUB>1</SUB>, C<SUB>2</SUB>, C<SUB>3</SUB>, and stores them in
	 * array C. Also calculates C<SUB>2</SUB>' and C<SUB>3</SUB>', and stores them
	 * in C2prime and C3prime, respectively. Assumes that the time from the root to
	 * each leaf is some constant t (i.e. the tree is ultrametric).
	 */
	private void calcC() {
		calcPartials();
		updateCs();
	}

	private void calcPartials() {
		List<Integer> seqLeavesLengths = tree.getSeqLeavesLengths();
		Iterator<Integer> leafIter = seqLeavesLengths.iterator();

		// estimate partials of G(z, t) with respect to z and evaluated at z = 1 and
		// time t,
		// by using the expected value of the kth factorial moment of the length process

		long[] curSum = new long[3];

		while (leafIter.hasNext()) {
			long length = (long) leafIter.next();
			long[] lengths = new long[] { length, length * (length - 1), length * (length - 1) * (length - 2) };

			for (int i = 0; i < 3; i++) {
				try {
					curSum[i] = Math.addExact(curSum[i], lengths[i]);
				} catch (ArithmeticException e) {
					partial[i] += (double) curSum[i] / seqLeavesLengths.size();
					curSum[i] = lengths[i];
				}
			}
		}

		// System.out.println("partial: " + partial[0]);
		// System.out.println("size: " + seqLeavesLengths.size());
		// System.out.println("PARTIALS");
		for (int i = 0; i < 3; i++) {
			partial[i] += (double) curSum[i] / seqLeavesLengths.size();
			// System.out.println(partial[i]);
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
			// System.out.println("here2");
			gamma = Double.NaN;
		}
	}

	/** Estimates beta. */
	private void estimateBeta() {
		beta = gamma * (2 + C2prime) - C2prime;

		try {
			beta /= 1 + gamma;
		} catch (ArithmeticException e) {
			// System.out.println("here3");
			beta = Double.NaN;
		}
	}

	/** Estimates M. */
	private void estimateM() {
		if (Double.isNaN(beta) || beta == 0)
			M = Double.NaN;
		else {
			// System.out.println("M calculation: " + C[0] + " " + beta);
			M = Math.round(C[0] / beta);
		}
	}

	/**
	 * Inverts the parameters of the length process, estimating gamma, beta, and M.
	 */
	private void invert() {
		calcC();
		estimateGamma();
		estimateBeta();
		estimateM();
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
