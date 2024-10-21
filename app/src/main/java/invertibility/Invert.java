package invertibility;

/**
 * Class used to invert the TKF91 process parameters, given the
 * sampled sequences.
 */
public abstract class Invert {
	protected LeafSamples tree;
	protected int N;

	protected double mu;
	protected double lambda;
	protected double gamma;
	protected double beta;
	protected Integer M;

	protected double nu;
	protected double a;
	protected double pi0;

	protected double[] partials;
	protected double[] C;

	public Invert() {
	}

	public Invert(LeafSamples tree) {
		this.tree = tree;
		N = tree.getN();

		partials = new double[3];
		C = new double[3];
	}

	abstract protected void calcPartials(); // implemented differently in Invert1Mer and InvertLength

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
					System.out.println("Error in partial computation.");
				}
			}
		}

		return partials;
	}
}
