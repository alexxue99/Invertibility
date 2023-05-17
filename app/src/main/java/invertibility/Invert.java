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
	private double discriminant;
	private boolean flip = false;
	private double a;
	private double pi0;

	private double[] partials;
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
		partials = new double[3];
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
		partials = new double[3];
		C = new double[3];
		D = new double[3];

		this.mu = mu;
		this.lambda = lambda;
		this.M = M;

		calcPartials("1-mer");
		updateCs();
		updateDs();

		estimateNu();
		estimatePi0();
		estimateA();
	}

	private void calcPartials(String type) {
		switch (type) {
			case "length":
				partials = calcPartials(tree.getSeqLeavesLengths());
				break;
			case "1-mer":
				int[] ones = tree.getSeqNumOnes();
				int[] zeros = tree.getSeqNumZeros();
				partials = calcPartials(ones, zeros, 0.1);
				break;
		}
	}

	private double[] calcPartials(int[] array) {
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

	private double[] calcPartials(int[] ones, int[] zeros, double PI0) {
		double[] partials_1 = calcPartials(ones);
		double[] partials_2 = calcPartials(zeros);
		
		double partial_12 = 0;
		double partial_112 = 0;
		double partial_122 = 0;

		for (int i = 0; i < N; i++) {
			partial_12 += (ones[i] * zeros[i] - partial_12) / (i + 1);
			partial_112 += (ones[i] * (ones[i] - 1) * zeros[i] - partial_112) / (i + 1);
			partial_122 += (ones[i] * zeros[i] * (zeros[i] - 1) - partial_122) / (i + 1);
		}

		double[] partials = new double[3];
		double PI1 = 1 - PI0;

		partials[0] = partials_1[0] * PI0 - partials_2[0] * PI1;
		partials[1] = partials_1[1] * PI0 * PI0 + partials_2[1] * PI1 * PI1 - 2 * partial_12 * PI0 * PI1;
		partials[2] = partials_1[2] * Math.pow(PI0, 3) - partials_2[2] * Math.pow(PI1, 3)
				- 3 * partial_112 * PI0 * PI0 * PI1 + 3 * partial_122 * PI0 * PI1 * PI1;

		return partials;
	}

	/* Prereq: partials are calculated */
	private void updateCs() {
		C[0] = partials[0];
		C[1] = partials[1] - partials[0] * partials[0];
		C[2] = partials[2] + 2 * Math.pow(partials[0], 3) - 3 * partials[0] * partials[1];
		if (C[0] == 0) {
			C2prime = Double.NaN;
			C3prime = Double.NaN;
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
		System.out.println(D[0] + "\t" + D[1] + "\t" + D[2]);
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
		discriminant = -3 * D[0] * D[0] * D[1] * D[1] + 4 * Math.pow(D[0], 3) * D[2] + 4 * M * D[1] * D[1]
				- 6 * M * D[0] * D[1] * D[2] + M * M * D[2] * D[2];
		exp = Math.sqrt(-discriminant);
		exp /= (D[0] * D[0] - M * D[1]);

		if (exp < 0) {
			exp *= -1;
			flip = true;
		}

		// System.out.println((D[0] * D[1] - M * D[2]) / (D[0] * D[0] - M * D[1]));
		System.out.println("EXP: " + exp + "\t" + discriminant);
		nu = -Math.log(exp);
	}

	private void estimatePi0() {
		pi0 = D[0] * D[1] - M * D[2];
		pi0 /= Math.sqrt(discriminant);

		if (flip)
			pi0 *= -1;

		pi0 += 0.5;
	}

	private void estimateA() {
		a = D[0] + M * (1 - pi0);
		a /= exp;
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
		return nu;
	}

	public double getPi0() {
		return pi0;
	}

	public double getA() {
		return a;
	}
}
