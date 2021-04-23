package invertibility;

import java.math.BigDecimal;
import java.math.RoundingMode;
import java.util.Iterator;
import java.util.List;

/**
 * Class used to invert the length process parameters, given the lengths of the
 * sequences at the leaves of a tree.
 */
public class Invert {
	final static int DIVISION_PRECISION = 8; // number of decimal digits in divisions
	private TreeLeaves tree;
	private BigDecimal[] C;
	private double C2prime;
	private double C3prime;
	private double gamma;
	private double beta;
	private double M;

	public Invert(TreeLeaves tree) {
		this.tree = tree;
		C = new BigDecimal[3];
		invert();
	}

	/**
	 * Calculates C<SUB>1</SUB>, C<SUB>2</SUB>, C<SUB>3</SUB>, and stores them in
	 * array C. Also calculates C<SUB>2</SUB>' and C<SUB>3</SUB>', and stores them
	 * in C2prime and C3prime, respectively. Assumes that the time from the root to
	 * each leaf is some constant t (i.e. the tree is ultrametric).
	 */
	private void calcC() {
		List<Integer> seqLeavesLengths = tree.getSeqLeavesLengths();
		Iterator<Integer> leafIter = seqLeavesLengths.iterator();

		// estimate partials of G(z, t) with respect to z and evaluated at z = 1 and
		// time t,
		// by using the expected value of the kth factorial moment of the length process
		BigDecimal[] partial = new BigDecimal[3];
		for (int k = 0; k < 3; k++)
			partial[k] = new BigDecimal(0);

		while (leafIter.hasNext()) {
			int length = leafIter.next();
			partial[0] = partial[0].add(new BigDecimal(length));
			partial[1] = partial[1].add(new BigDecimal(length * (length - 1)));
			partial[2] = partial[2].add(new BigDecimal(length * (length - 1) * (length - 2)));
		}

		for (int k = 0; k < 3; k++)
			partial[k] = partial[k].divide(new BigDecimal(seqLeavesLengths.size()), DIVISION_PRECISION, RoundingMode.HALF_UP);

		C[0] = partial[0];
		C[1] = partial[1].subtract(partial[0].pow(2));
		C[2] = partial[2].add(new BigDecimal(2).multiply(partial[0].pow(3)))
				.subtract(new BigDecimal(3).multiply(partial[0].multiply(partial[1])));

		if (C[0].equals(new BigDecimal(0))) {
			C2prime = Double.NaN;
			C3prime = Double.NaN;
		} else {
			C2prime = C[1].divide(C[0], DIVISION_PRECISION, RoundingMode.HALF_UP).doubleValue();
			C3prime = C[2].divide(C[0], DIVISION_PRECISION, RoundingMode.HALF_UP).doubleValue();
		}
	}

	/** Estimates gamma. */
	private void estimateGamma() {
		gamma = Math.sqrt(-(C2prime + 1) * (C2prime + 1) * (3 * C2prime * C2prime - 2 * C3prime));
		gamma += -C2prime * C2prime + C2prime + C3prime;

		try {
			gamma /= 2 * C2prime * C2prime + 2 * C2prime - C3prime + 2;
		} catch (ArithmeticException e) {
			gamma = Double.NaN;
		}
	}

	/** Estimates beta. */
	private void estimateBeta() {
		beta = gamma * (2 + C2prime) - C2prime;

		try {
			beta /= 1 + gamma;
		} catch (ArithmeticException e) {
			beta = Double.NaN;
		}
	}

	/** Estimates M. */
	private void estimateM() {
		if (Double.isNaN(beta) || beta == 0)
			M = Double.NaN;
		else
			M = C[0].divide(new BigDecimal(beta), DIVISION_PRECISION, RoundingMode.HALF_UP).doubleValue();
	}

	/**
	 * Inverts the parameters of the length process, estimating gamma, beta, and M.
	 */
	public void invert() {
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
}
