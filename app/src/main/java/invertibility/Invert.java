package invertibility;

/**
 * Class used to invert the length process parameters, given the lengths of the
 * sequences at the leaves of a tree.
 */
public abstract class Invert {
	protected TreeLeaves tree;
	protected int N;

	protected double mu;
	protected double lambda;
	protected double gamma;
	protected double beta;
	protected double M;

	protected double exp;
	protected double nu;
	protected double a;
	protected double pi0;

	protected double[] partials;
	protected double[] C;

	public Invert(TreeLeaves tree) {
		this.tree = tree;
		N = tree.getN();

		partials = new double[3];
		C = new double[3];
	}

	abstract protected void calcPartials();

	protected double[] factorialMoments(int[] array) {
		// estimate partials of G(z, t) with respect to z and evaluated at z = 1 and
		// time t,
		// by using the expected value of the kth factorial moment of the length process
		double[] partials = new double[3];
		for (int i = 0; i < N;) {
			long length = array[i];
			long[] lengths = new long[] { length, length * (length - 1), length * (length - 1) * (length - 2) };

			i++;
			for (int j = 0; j < 3; j++) {
				try {
					partials[j] += (lengths[j] - partials[j]) / i;
				} catch (ArithmeticException e) {
					System.out.println("PARTIALS");
				}
			}
		}

		return partials;
	}

	/* Prereq: partials are calculated */
	protected void updateCs() {
		// System.out.println("PARTIALS: " + partials[0] + " " + partials[1] + " " +
		// partials[2]);
		C[0] = partials[0];
		C[1] = partials[1] - partials[0] * partials[0];
		C[2] = partials[2] + 2 * Math.pow(partials[0], 3) - 3 * partials[0] * partials[1];
	}
}
