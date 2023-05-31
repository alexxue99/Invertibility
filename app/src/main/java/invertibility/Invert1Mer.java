package invertibility;

public class Invert1Mer extends Invert {
    private double discriminant;
	private boolean flip = false;
    
	private double[] D;
	/**
	 * Inverts the process based on given data on the lengths of the sequences at
	 * the leaves, as well as on mu, lambda, and M, to estimate nu, pi, and a.
	 */
	public Invert1Mer(double mu, double lambda, int M, TreeLeaves tree) {
		this.tree = tree;
		N = tree.getN();
		partials = new double[3];
		C = new double[3];
		D = new double[3];

		this.mu = mu;
		this.lambda = lambda;
		this.M = M;

		calcPartials();
		updateCs();
		updateDs();

		estimateNu();
		estimatePi0();
		estimateA();
	}

    protected void calcPartials() {
        int[] ones = tree.getSeqNumOnes();
        int[] zeros = tree.getSeqNumZeros();

		double[] partials_1 = factorialMoments(ones);
		double[] partials_2 = factorialMoments(zeros);
		
		double partial_12 = 0;
		double partial_112 = 0;
		double partial_122 = 0;

		for (int i = 0; i < N; i++) {
			partial_12 += (ones[i] * zeros[i] - partial_12) / (i + 1);
			partial_112 += (ones[i] * (ones[i] - 1) * zeros[i] - partial_112) / (i + 1);
			partial_122 += (ones[i] * zeros[i] * (zeros[i] - 1) - partial_122) / (i + 1);
		}

        double PI0 = 0.3;
		double PI1 = 1 - PI0;

		partials[0] = partials_1[0] * PI0 - partials_2[0] * PI1;
		partials[1] = partials_1[1] * PI0 * PI0 + partials_2[1] * PI1 * PI1 - 2 * partial_12 * PI0 * PI1;
		partials[2] = partials_1[2] * Math.pow(PI0, 3) - partials_2[2] * Math.pow(PI1, 3)
				- 3 * partial_112 * PI0 * PI0 * PI1 + 3 * partial_122 * PI0 * PI1 * PI1;

		System.out.println("PARTIALS: " + partials[0] + " " + partials[1] + " " + partials[2]);
	}
    
	/* Prereq: Cs are calculated */
	private void updateDs() {
		D[0] = C[0] * Math.exp(mu);
		D[1] = -C[1] * Math.exp(2 * mu);
		D[2] = C[2] * Math.exp(3 * mu) / 2;
		System.out.println(D[0] + "\t" + D[1] + "\t" + D[2]);
	}
 	
	private void estimateNu() {
		discriminant = -3 * D[0] * D[0] * D[1] * D[1] + 4 * Math.pow(D[0], 3) * D[2] + 4 * M * D[1] * D[1]
				- 6 * M * D[0] * D[1] * D[2] + M * M * D[2] * D[2];
		exp = Math.sqrt(discriminant);
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




