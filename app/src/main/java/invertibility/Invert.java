package invertibility;

/**
 * Class used to invert the length process parameters, given the lengths of the
 * sequences at the leaves of a tree.
 */
public class Invert {
	final static int DIVISION_PRECISION = 8; // number of decimal digits in divisions
	private TreeLeaves tree;
	private int N;

	private double mu;
	private double lambda;
	private double gamma;
	private double beta;
	private double M;

	private double exp;
	private double nu;
	private double aUpper;
	private double aLower;
	private double pi0Upper;
	private double pi0Lower;

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
		N = tree.getN();
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
	 * the leaves, as well as on mu, lambda, and M, to estimate nu, pi, and a.
	 */
	public Invert(double mu, double lambda, double M, TreeLeaves tree) {
		this.tree = tree;
		N = tree.getN();
		partial = new double[3];
		C = new double[3];
		D = new double[3];

		this.mu = mu;
		this.lambda = lambda;
		this.M = M;

		calcPartials("ones");
		updateCs();
		
		for (int i = 0; i < 3; i++)
			D[i] = C[i];

		calcPartials("zeros");
		updateCs();

		// for (int i = 0; i < 3; i++) {
		//	C[i] = C[i] * (pi0 - 1) + D[i] * pi0;
		// }

		updateDs();
		estimateNu();
		estimatePi0();
		estimateA();
	}

	private void calcPartials(String type) {
		int[] array = null;

		switch (type) {
			case "length":
				array = tree.getSeqLeavesLengths();
				calcPartials(array);
				break;
			case "ones":
				array = tree.getSeqNumOnes();
				// d ln G / dz1 * pi0  - d ln G / dz2 * pi1
				break;
			case "zeros":
				array = tree.getSeqNumZeros();
				break;
		}

		calcPartials(array);
	}

	private void calcPartials(int[] seqLeavesLengths) {
		// estimate partials of G(z, t) with respect to z and evaluated at z = 1 and
		// time t,
		// by using the expected value of the kth factorial moment of the length process

		for (int i = 0 ; i < N;) {
			long length = seqLeavesLengths[i];
			long[] lengths = new long[] { length, length * (length - 1), length * (length - 1) * (length - 2) };

			i++;
			for (int j = 0; j < 3; j++) {
				try {
					partial[j] += (lengths[j] - partial[j]) / i;
				} catch (ArithmeticException e) {
					System.out.println("PARTIALS");
				}
			}
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
		D[1] = -C[1] * Math.exp(2 * mu);
		D[2] = C[2] * Math.exp(3 * mu) / 2;
	}

	/** Estimates gamma. */
	private void estimateGamma() {
		gamma = Math.sqrt(-(C2prime + 1) * (C2prime + 1) * (3 * C2prime * C2prime - 2 * C3prime));

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

	private void estimateNu() {
		exp = (D[0] * D[1] - M * D[2]) / (D[0] * D[0] - M * D[1]);
		nu = -Math.log(exp);
	}

	private void estimatePi0() {
		double firstNumerator = D[0] * D[1] - D[2] * M;

		double secondNumerator = Math.sqrt(-3 * D[0] * D[0] * D[1] * D[1] + 4 * Math.pow(D[0], 3) * D[2]
				+ 4 * Math.pow(D[1], 3) * M
				- 6 * D[0] * D[1] * D[2] * M + D[2] * D[2] * M * M);

		double denominator = 2 * D[0] * D[1] - 2 * D[2] * M;

		pi0Upper = (firstNumerator + secondNumerator) / denominator;
		pi0Lower = (firstNumerator - secondNumerator) / denominator;
	}

	private void estimateA() {
		aUpper = (D[0] - M * (1 - pi0Upper) * exp) / ((2 * pi0Upper - 1) * exp);
		aLower = (D[0] - M * (1 - pi0Lower) * exp) / ((2 * pi0Lower - 1) * exp);
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

	public double getNu() {
		// TODO
		return exp;
	}

	public double getPi0Upper() {
		return pi0Upper;
	}

	public double getPi0Lower() {
		return pi0Lower;
	}

	public double getAUpper() {
		return aUpper;
	}

	public double getALower() {
		return aLower;
	}
}
