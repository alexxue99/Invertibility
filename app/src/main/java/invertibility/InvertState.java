package invertibility;

import org.ejml.simple.*;

public class InvertState extends Invert {

	/**
	 * Inverts the process based on given data on the sequences at
	 * the leaves, as well as on mu, lambda, M, nu and pi0 to estimate the initial
	 * state.
	 */
	public InvertState(double mu, double lambda, int M, double nu, double pi0, TreeLeaves tree) {
		super(tree);

		this.mu = mu;
		this.lambda = lambda;
		this.M = M;
		this.nu = nu;
		this.pi0 = pi0;

		gamma = lambda / mu;
		beta = Math.exp(lambda - mu);
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

	private double estimateP(double t) {
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
}
