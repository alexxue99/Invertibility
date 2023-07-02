package invertibility;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;

import org.ejml.data.*;
import org.ejml.dense.row.decomposition.lu.*;
import org.ejml.interfaces.decomposition.*;
import org.ejml.simple.*;

import invertibility.Invert;
import invertibility.TreeLeaves;

public class InvertState extends Invert {
	private String rootState;

	/**
	 * Inverts the process based on given data on the sequences at
	 * the leaves, as well as on mu, lambda, M, nu and pi0 to estimate the initial
	 * state.
	 */
	public InvertState(double lambda, double mu, double nu, int M, double pi0, TreeLeaves tree, int method) {
		super(tree);

		this.lambda = lambda;
		this.mu = mu;
		this.nu = nu;
		this.M = M;
		this.pi0 = pi0;

		gamma = lambda / mu;
		beta = Math.exp(lambda - mu);

		rootState = "";

		switch (method) {
			case 1:
				estimateRootState(); // using inverse
				break;
			case 2:
				estimateRootState2(); // using LDU decomposition
				break;
			case 3:
				estimateRootState3(); // exhaustive search
				break;
			case 4:
				estimateRootState4(); // greedy algorithm (suboptimal results)
				break;
			case 5:
				estimateRootState5(); // greedy algorithm under certain conditions (mathematically correct)
				break;
			case 6:
				estimateRootState6(); // uses subsetsum approximation
				break;
		}
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
				sum2 += (s.charAt(i) == '1' ^ zero) ? powerEta : 0;
				powerEta *= eta;
			}

			sum += sum2 * psi;
		}

		return sum / seq.length;
	}

	private void check() {
		// double[][] Y = { { 0, 1 }, { 0, 1 }, { 1, 0 }, { 0, 1 }, { 1, 0 }, { 0, 1 },
		// { 0, 1 }, { 0, 1 } };
		double[][] Y = { { 0, 1 }, { 0, 1 }, { 1, 0 }, { 0, 1 } };

		double[][] V = new double[M][M];
		for (int i = 0; i < M; i++) {
			double eta = eta(1.01 + i / 1.0);

			V[i][0] = 1;
			for (int j = 1; j < M; j++) {
				V[i][j] = V[i][j - 1] * eta;
			}
		}

		double[] Psi = new double[M];
		for (int i = 0; i < M; i++) {
			Psi[i] = psi(1.01 + i / 1.0);
		}

		SimpleMatrix ans = (SimpleMatrix.diag(Psi)).mult(new SimpleMatrix(V)).mult(new SimpleMatrix(Y));
		System.out.println("--------");
		ans.print();
		// U().print();
	}

	private void check2() {
		System.out.println("p: " + p(1, true));
		System.out.println(pi0 * phi(1) * (1 - eta(1)) + psi(1));
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

	private SimpleMatrix VInverse() {
		double[][] V = new double[M][M];
		for (int i = 0; i < M; i++) {
			double eta = eta(1.01 + i / 1.0);

			V[i][0] = 1;
			for (int j = 1; j < M; j++) {
				V[i][j] = V[i][j - 1] * eta;
			}
		}
		// ((new SimpleMatrix(V)).invert().mult(new SimpleMatrix(V))).print();
		return (new SimpleMatrix(V)).invert();
	}

	private SimpleMatrix Psi() {
		double[] Psi = new double[M];
		for (int i = 0; i < M; i++) {
			Psi[i] = psi(1.01 + i / 1.0);
		}

		return SimpleMatrix.diag(Psi);
	}

	private SimpleMatrix PsiInverse() {
		double[] PsiInverse = new double[M];
		for (int i = 0; i < M; i++) {
			PsiInverse[i] = 1 / psi(1.01 + i / 1.0);
		}

		return SimpleMatrix.diag(PsiInverse);
	}

	private SimpleMatrix U() {
		double[][] U = new double[M][2];

		for (int i = 0; i < M; i++) {
			double p0 = p(1.01 + i / 1.0, true);
			double p1 = p(1.01 + i / 1.0, false);

			// System.out.println((p0 + p1) + "\t" + (1 - Math.pow(eta(1.01 + i / 1.0),
			// M)));
			double secondTerm = phi(1.01 + i / 1.0) * (1 - Math.pow(eta(1.01 + i / 1.0), M));

			U[i][0] = p0 - pi0 * secondTerm;
			U[i][1] = p1 - (1 - pi0) * secondTerm;
		}

		return new SimpleMatrix(U);
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

	private SimpleMatrix Y() {
		SimpleMatrix ans = VInverse().mult((PsiInverse().mult(U())));
		// ans.print();
		return ans;
	}

	private void estimateRootState2() {
		DMatrixRMaj matrix = Psi().mult(V()).getMatrix();
		double[] U = Uvec();

		LUDecompositionAlt_DDRM decomp = new LUDecompositionAlt_DDRM();
		decomp.decompose(matrix);

		DMatrixRMaj lower = decomp.getLower(null);
		DMatrixRMaj upper = decomp.getUpper(null);
		int[] pivot = decomp.getRowPivotV(null);

		SimpleMatrix l = SimpleMatrix.wrap(lower);
		SimpleMatrix u = SimpleMatrix.wrap(upper);

		// solve Ly = b

		double[] y = new double[M];
		for (int i = 0; i < M; i++) {
			double sum = 0;
			for (int j = 0; j < i; j++) {
				sum += y[j] * lower.get(i, j);
			}
			y[i] = (U[i] - sum) / lower.get(i, i);
		}

		// solve Ux = y
		double[] x = new double[M];
		for (int i = M - 1; i >= 0; i--) {
			double sum = 0;
			for (int j = i + 1; j < M; j++) {
				sum += x[j] * upper.get(i, j);
			}
			x[i] = (y[i] - sum) / upper.get(i, i);
		}

		// pivot
		for (int i = 0; i < pivot.length; i++) {
			System.out.println(pivot[i]);
		}

		System.out.println("x");
		for (int i = 0; i < M; i++) {
			System.out.println(x[pivot[i]]);
		}

		SimpleMatrix prod = l.mult(u);
		for (int i = 0; i < M; i++) {
			double sum = 0;
			for (int j = 0; j < M; j++) {
				sum += prod.get(i, j) * x[j];
			}
			System.out.println(sum);
		}

		System.out.println("U");
		for (int i = 0; i < M; i++) {
			System.out.println(U[i]);
		}

		for (int i = 0; i < M; i++) {
			rootState += (x[pivot[i]] - 0.5 > 0) ? '0' : '1';
		}
	}

	private double dist(boolean[] Y) {
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

	private void estimateRootState3() {
		root = new boolean[M];
		Y = new boolean[M];

		prod = Psi().mult(V());
		U = Uvec();

		min = Double.MAX_VALUE;
		estimateRootState3(0);

		System.out.println(min);
		for (int i = 0; i < M; i++) {
			rootState += (Y[i]) ? '0' : '1';
		}
	}

	// calculates min among all possiblities of root state, where everything from
	// 0...index-1 is predetermined
	private void estimateRootState3(int index) {
		if (index == M) {
			if (dist(root) < min) {
				min = dist(root);
				Y = Arrays.copyOf(root, M);
			}
			return;
		}

		root[index] = false;
		estimateRootState3(index + 1);

		root[index] = true;
		estimateRootState3(index + 1);
	}

	private void estimateRootState4() {
		root = new boolean[M];

		prod = Psi().mult(V());
		U = Uvec();

		min = Double.MAX_VALUE;
		for (int i = 0; i < M; i++) {
			root[i] = false;
			if (dist(root) < min) {
				min = dist(root);
			}

			root[i] = true;
			if (dist(root) < min) {
				min = dist(root);
			} else {
				root[i] = false;
			}
		}

		System.out.println(min);
		for (int i = 0; i < M; i++) {
			rootState += (root[i]) ? '0' : '1';
		}
	}

	// use if gamma > 1 and beta(2-g) < 1 (e.x. when gamma > 2)
	// OR if gamma < 1 and beta(2-g) > 1
	private void estimateRootState5() {
		root = new boolean[M];

		prod = Psi().mult(V());
		U = Uvec();
		prod.print();
		System.out.println("U0: " + U[0]);

		double cumSum = 0;
		for (int i = 0; i < M; i++) {
			if (cumSum + prod.get(0, i) < U[0]) {
				root[i] = true;
				cumSum += prod.get(0, i);
			}
		}

		for (int i = 0; i < M; i++) {
			rootState += (root[i]) ? '0' : '1';
		}
	}

	// uses subset-sum approximation with ratio epsilon = 1/10
	private void estimateRootState6() {
		double epsilon = 0.1;
		prod = Psi().mult(V());

		HashMap<Double, boolean[]> L = new HashMap<Double, boolean[]>();
		L.put(Double.valueOf(0), new boolean[M]);
		double T = Uvec()[0];
		System.out.println("T: " + T);
		for (int i = 0; i < M; i++) {
			HashMap<Double, boolean[]> U = new HashMap<Double, boolean[]>();
			for (double y : L.keySet()) {
				U.put(y, L.get(y));

				if (y + prod.get(0, i) <= T) {
					boolean[] copy = Arrays.copyOf(L.get(y), M);
					copy[i] = true;
					U.put(y + prod.get(0, i), copy);
				}
			}

			L = new HashMap<Double, boolean[]>();
			L.put(Double.valueOf(0), new boolean[M]);

			double y = 0;
			List<Double> keys = new ArrayList<Double>(U.keySet());
			Collections.sort(keys);

			for (double z : keys) {
				if (y + epsilon * T / M < z) {
					y = z;
					L.put(z, U.get(z));
				}
			}
		}

		for (double y : L.keySet()) {
			System.out.println(y);
			for (int j = 0; j < M; j++) {
				System.out.print(((L.get(y)[j]) ? '1' : '0'));
			}
			System.out.println("\n");
		}

		double max = 0;
		for (double y : L.keySet()) {
			if (y > max)
				max = y;
		}

		boolean[] root = L.get(max);

		for (int i = 0; i < M; i++) {
			rootState += (root[i]) ? '0' : '1';
		}
	}

	private void estimateRootState() {
		SimpleMatrix Y = Y();
		// System.out.println("Y: " );
		// Y.print();
		for (int i = 0; i < M; i++) {
			rootState += (Y.get(i, 0) > Y.get(i, 1)) ? '0' : '1';
		}
		// System.out.println(Y.get(0, 0) - Y.get(0, 1));
	}

	public String getRootState() {
		// check();
		// System.out.println("U: " );
		// U().print();
		return rootState;
	}
}
