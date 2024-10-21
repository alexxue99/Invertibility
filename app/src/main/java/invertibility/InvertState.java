package invertibility;

import java.util.Arrays;

import org.ejml.simple.SimpleMatrix;

/**
 * Class used to invert the TKF91 process based on given data on the sampled
 * sequences, as well as on mu, lambda, M, nu and pi0 to estimate the initial
 * state.
 */
public class InvertState extends Invert {
	private String rootState;

	/**
	 * Constructor to initialize an InvertState object.
	 * 
	 * @param lambda the insertion rate
	 * @param mu     the deletion rate
	 * @param nu     the substitution rate
	 * @param M      the length of the root sequence
	 * @param pi0    the probability a character is 0 after substitution or
	 *               insertion
	 * @param tree   a LeafSamples object containing the sampled sequences
	 */
	public InvertState(double lambda, double mu, double nu, int M, double pi0, LeafSamples tree) {
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

	private double p(double t, boolean zero) {
		double sum = 0;

		double eta = eta(t - 1);
		double psi = psi(t - 1);
		double phi = phi(t - 1);

		String[] seq = tree.getSeq();
		for (String s : seq) {
			sum += ((zero) ? pi0 : (1 - pi0)) * phi * (1 - Math.pow(eta, s.length()));

			double sum2 = 0;
			double powerEta = 1;
			for (int i = 0; i < s.length(); i++) {
				sum2 += ((s.charAt(i) == '1') ^ zero) ? powerEta : 0;
				powerEta *= eta;
			}

			sum += sum2 * psi;
		}

		return sum / seq.length;
	}

	private SimpleMatrix V() {
		double[][] V = new double[M][M];
		for (int i = 0; i < M; i++) {
			double eta = eta(1.01 + i / 1.0);

			V[i][0] = 1;
			for (int j = 1; j < M; j++) {
				V[i][j] = V[i][j - 1] * eta;
			}
		}
		return new SimpleMatrix(V);
	}

	private SimpleMatrix Psi() {
		double[] Psi = new double[M];
		for (int i = 0; i < M; i++) {
			Psi[i] = psi(1.01 + i / 1.0);
		}

		return SimpleMatrix.diag(Psi);
	}

	private double[] Uvec() {
		double[] U = new double[M];

		for (int i = 0; i < M; i++) {
			double p0 = p(1.01 + i / 1.0, true);
			double secondTerm = phi(1.01 + i / 1.0) * (1 - Math.pow(eta(1.01 + i / 1.0), M));

			U[i] = p0 - pi0 * secondTerm;
		}

		return U;
	}

	// calculates length of U0 - Psi V Y0
	private double length(boolean[] Y) {
		double dist = 0;
		for (int i = 0; i < M; i++) {
			double sum = 0;
			for (int j = 0; j < M; j++) {
				if (Y[j])
					sum += prod.get(i, j);
			}
			sum -= U[i];
			dist += sum * sum;
		}

		return dist;
	}

	private SimpleMatrix prod;
	private double[] U;
	private boolean[] root;
	private boolean[] Y;
	private double min;

	private void estimateRootState() {
		root = new boolean[M];
		Y = new boolean[M];

		prod = Psi().mult(V());
		U = Uvec();

		min = Double.MAX_VALUE;
		estimateRootState(0);

		for (int i = 0; i < M; i++) {
			rootState += (Y[i]) ? '0' : '1';
		}
	}

	private void estimateRootState(int index) {
		if (index == M) {
			if (length(root) < min) {
				min = length(root);
				Y = Arrays.copyOf(root, M);
			}
			return;
		}

		root[index] = false;
		estimateRootState(index + 1);

		root[index] = true;
		estimateRootState(index + 1);
	}

	public String getRootState() {
		return rootState;
	}
}
