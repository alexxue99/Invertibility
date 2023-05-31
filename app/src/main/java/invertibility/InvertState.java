package invertibility;

import org.ejml.simple.*;

public class InvertState extends Invert {
	private String rootState;
	/**
	 * Inverts the process based on given data on the sequences at
	 * the leaves, as well as on mu, lambda, M, nu and pi0 to estimate the initial
	 * state.
	 */
	public InvertState(double lambda, double mu, double nu, int M, double pi0, TreeLeaves tree) {
		super(tree);

		this.lambda = lambda;
		this.mu = mu;
		this.nu = nu;
		this.M = M;
		this.pi0 = pi0;

		gamma = lambda / mu;
		beta = Math.exp(lambda - mu);

        rootState = "";
		estimateRootState();
	}

	protected void calcPartials() {
		return;
	}

	private double eta(double t) {
		return (1 - Math.pow(beta, t)) / (1 - gamma * Math.pow(beta, t));
	}

	private double psi(double t) {
		return Math.exp(-(mu + nu) * t);
	}

	private double phi(double t) {
		return (-psi(t) + 1 - eta(t)) / (1 - eta(t));
	}

	private double p(double t) {
		double sum = 0;

		double eta = eta(t);
		double psi = psi(t);
		double phi = phi(t);

		String[] seq = tree.getSeq();
		for (String s : seq) {
			sum += pi0 * phi * (1 - Math.pow(eta, s.length()));

			double sum2 = 0;
			double powerEta = 1;
			for (int i = 0; i < s.length(); i++) {
				sum2 += (s.charAt(i) == '0') ? powerEta : 0;
				powerEta *= eta;
			}

			sum += sum2 * psi;
		}

		return sum / seq.length;
	}

	private SimpleMatrix VInverse() {
		double [][] V = new double[M][M];
		for (int i = 0; i < M; i++) {
			double eta = eta(1.01 + i / 100.0);

			V[i][0] = 1;
			for (int j = 1; j < M; j++) {
				V[i][j] = V[i][j-1] * eta;
			}
		}

		return (new SimpleMatrix(V)).invert();
	}

	private SimpleMatrix PsiInverse() {
		double [] PsiInverse = new double[M];
		for (int i = 0; i < M; i++) {
			PsiInverse[i] = 1 / psi(1.01 + i / 100.0);
		}

		return SimpleMatrix.diag(PsiInverse);
	}

	private SimpleMatrix U() {
		double [][] U = new double[M][2];

		for (int i = 0; i < M; i++) {
			double p = p(1.01 + i / 100.0);
			double secondTerm = phi(1.01 + i / 100.0) * (1 - Math.pow(eta(1.01 + i / 100.0), M));
	
			U[i][0] = p - pi0 * secondTerm;
			U[i][1] = (1 - p) - (1 - pi0) * secondTerm;
		}

		return new SimpleMatrix(U);
	}

	private SimpleMatrix Y() {
		return VInverse().mult((PsiInverse().mult(U())));
	}

	private void estimateRootState() {
		SimpleMatrix Y = Y();
		for (int i = 0; i < M; i++) {
			rootState += (Y.get(i, 0) > Y.get(i, 1)) ? '0' : '1';
		}
	}

    public String getRootState() {
        return rootState;
    }
}
