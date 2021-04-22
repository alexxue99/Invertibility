package invertibility;

import java.util.Iterator;
import java.util.LinkedList;

import invertibility.TreeSimul.Leaf;

public class Invert {
	private TreeSimul tree;
	private double[] C;
	private double C2prime;
	private double C3prime;
	private double gamma;
	private double beta;
	private double nuT;
	private double pi0;
	private double M;

	public Invert(TreeSimul tree) {
		this.tree = tree;
		C = new double[3];
		invert();
	}

	/**
	 * Calculates C<SUB>1</SUB>, C<SUB>2</SUB>, C<SUB>3</SUB>, and stores them in
	 * array C. Also calculates C<SUB>2</SUB>' and C<SUB>3</SUB>', and stores them
	 * in C2prime and C3prime, respectively. Assumes that the time from the root to
	 * each leaf is some constant t (i.e. the tree is ultrametric).
	 */
	private void calcC() {
		LinkedList<Leaf> seqLeaves = tree.getSeqLeaves();
		Iterator<Leaf> leafIter = seqLeaves.iterator();

		// estimate partials of G(z, t) with respect to z and evaluated at z = 1 and
		// time t,
		// by using the expected value of the kth factorial moment of the length process
		double[] partials = new double[3];

		while (leafIter.hasNext()) {
			int length = leafIter.next().getSequence().length();
			partials[0] += length;
			partials[1] += length * (length - 1);
			partials[2] += length * (length - 1) * (length - 2);
		}

		for (int k = 0; k < 3; k++)
			partials[k] /= seqLeaves.size();

		C[0] = partials[0];
		C[1] = partials[1] - partials[0] * partials[0];
		C[2] = partials[2] + 2 * Math.pow(partials[0], 3) - 3 * partials[0] * partials[1];

		C2prime = C[1] / C[0];
		C3prime = C[2] / C[0];
	}

	/** Estimates gamma. */
	private void estimateGamma() {
		gamma = Math.sqrt(-(C2prime + 1) * (C2prime + 1) * (3 * C2prime * C2prime - 2 * C3prime));
		gamma += -C2prime * C2prime + C2prime + C3prime;
		gamma /= 2 * C2prime * C2prime + 2 * C2prime - C3prime + 2;
	}

	/** Estimates beta. */
	private void estimateBeta() {
		beta = gamma * (2 + C2prime) - C2prime;
		beta /= 1 + gamma;
	}

	/** Estimates M. */
	private void estimateM() {
		M = C[0] / beta;
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

	public double getNuT() {
		return nuT;
	}

	public double getPi0() {
		return pi0;
	}

	public double getM() {
		return M;
	}
}
