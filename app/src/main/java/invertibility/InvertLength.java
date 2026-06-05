package invertibility;

/**
 * * Class used to invert the length process based on given data on the lengths
 * of * the sampled sequences to estimate gamma, beta, and M.
 */
public class InvertLength extends Invert {
        private double C2prime;
        private double C3prime;

        /**
         * Constructor to initialize an InvertLength object. 
         * @param tree a
         * LeafSamples object containing the sampled sequences
         */
        public InvertLength(LeafSamples tree) {
                super(tree);
                calcPartials();
                updateCs();
                updateCprimes();
                estimateGamma();
                estimateBeta();
                estimateM();
        }

        protected void calcPartials() {
                partials = factorialMoments(tree.getSamplesLengths());
        }

        private void updateCs() {
                C[0] = partials[0];
                C[1] = partials[1] - partials[0] * partials[0];
                C[2] = partials[2] + 2 * Math.pow(partials[0], 3) - 3 * partials[0] * partials[1];
        }

        private void updateCprimes() {
                if (C[0] == 0) {
                        C2prime = Double.NaN;
                        C3prime = Double.NaN;
                } else {
                        C2prime = C[1] / C[0];
                        C3prime = C[2] / C[0];
                }
        }

        /** Estimates gamma. */
        private void estimateGamma() {
                gamma = Math.sqrt(-(C2prime + 1) * (C2prime + 1) * (3 * C2prime * C2prime - 2 * C3prime));
                gamma += -C2prime * C2prime + C2prime + C3prime;
                try {
                        gamma /= 2 * C2prime * C2prime + 2 * C2prime - C3prime + 2;
                } catch (ArithmeticException e) {
                        gamma = Double.NaN;
                }
        }

        /** Estimates beta. */
        private void estimateBeta() {
                beta = gamma * (2 + C2prime) - C2prime;
                try {
                        beta /= 1 + gamma;
                } catch (ArithmeticException e) {
                        beta = Double.NaN;
                }
        }

        /** Estimates M. */
        private void estimateM() {
                if (Double.isNaN(beta) || beta == 0) {
                        M = null;
                } else {
                        Long val = Math.round(C[0] / beta);
                        M = val.intValue();
                }
        }

        public double getGamma() {
                return gamma;
        }

        public double getBeta() {
                return beta;
        }

        public Integer getM() {
                return M;
        }

        public double getLambda() {
                if (Double.isNaN(beta) || beta == 0 || Double.isNaN(gamma)) {
                        return Double.NaN;
                } else {
                        return Math.log(beta) + Math.log(beta) / (gamma - 1);
                }
        }

        public double getMu() {
                if (Double.isNaN(beta) || beta == 0 || Double.isNaN(gamma)) {
                        return Double.NaN;
                } else {
                        return Math.log(beta) / (gamma - 1);
                }
        }
}

// package invertibility;

// import org.apache.commons.math3.analysis.MultivariateFunction;
// import org.apache.commons.math3.optim.InitialGuess;
// import org.apache.commons.math3.optim.MaxEval;
// import org.apache.commons.math3.optim.PointValuePair;
// import org.apache.commons.math3.optim.SimpleBounds;
// import org.apache.commons.math3.optim.nonlinear.scalar.GoalType;
// import org.apache.commons.math3.optim.nonlinear.scalar.ObjectiveFunction;
// import
// org.apache.commons.math3.optim.nonlinear.scalar.noderiv.BOBYQAOptimizer;

// /**
// * Class used to invert the length process based on given data on the lengths
// of
// * the sampled sequences to estimate gamma, beta, and M.
// *
// * Uses nonlinear least squares instead of explicit symbolic inversion.
// */
// public class InvertLength extends Invert {

// /**
// * Constructor to initialize an InvertLength object.
// *
// * @param tree a LeafSamples object containing the sampled sequences
// */
// public InvertLength(LeafSamples tree) {
// super(tree);

// calcPartials();
// updateCs();

// estimateParameters();
// }

// protected void calcPartials() {
// partials = factorialMoments(tree.getSamplesLengths());
// }

// /**
// * Computes empirical C1, C2, C3 from sample moments.
// */
// private void updateCs() {
// C[0] = partials[0];

// // variance-like quantity
// C[1] = partials[1] - partials[0] * partials[0];

// // centered third moment-like quantity
// C[2] = partials[2]
// + 2 * Math.pow(partials[0], 3)
// - 3 * partials[0] * partials[1];
// }

// /**
// * Estimate parameters via nonlinear least squares.
// */
// private void estimateParameters() {
// MultivariateFunction loss = point -> {
// double betaGuess = point[0];
// double gammaGuess = point[1];

// // invalid region
// if (betaGuess <= 0 || gammaGuess <= 0
// || gammaGuess >= 1) {
// return Double.POSITIVE_INFINITY;
// }

// // theoretical moments
// double C2theory = C[0]
// * (2 * gammaGuess
// - betaGuess * gammaGuess
// - betaGuess)
// / (1 - gammaGuess);

// double numerator = 2 * C[0]
// * (betaGuess * betaGuess * gammaGuess * gammaGuess
// + betaGuess * betaGuess * gammaGuess
// + betaGuess * betaGuess
// - 3 * betaGuess * gammaGuess * gammaGuess
// - 3 * betaGuess * gammaGuess
// + 3 * gammaGuess * gammaGuess);

// double denominator = (1 - gammaGuess) * (1 - gammaGuess);

// double C3theory = numerator / denominator;

// /*
// * Weighted least squares loss.
// *
// * Third moments are noisier, so downweight them slightly.
// */

// double w1 = 1.0;
// double w2 = 0.1;

// double lossValue = w1 * Math.pow(C2theory - C[1], 2)
// + w2 * Math.pow(C3theory - C[2], 2);

// return lossValue;
// };

// /*
// * Initial guesses.
// *
// * These do not need to be perfect.
// */

// double[] initialGuess = new double[] {
// // M
// 0.9, // beta
// 0.5 // gamma
// };

// /*
// * Bounds:
// *
// * M > 0
// * beta > 0
// * 0 < gamma < 1
// */

// double[] lowerBounds = new double[] {
// //1e-6,
// 1e-6,
// 1e-6
// };

// double[] upperBounds = new double[] {
// // 1e9,
// 1e3,
// 0.999999
// };

// try {

// BOBYQAOptimizer optimizer = new BOBYQAOptimizer(5);

// PointValuePair result = optimizer.optimize(
// new MaxEval(10000),
// new ObjectiveFunction(loss),
// GoalType.MINIMIZE,
// new InitialGuess(initialGuess),
// new SimpleBounds(lowerBounds, upperBounds));

// double[] solution = result.getPoint();

// beta = solution[0];
// gamma = solution[1];
// M = (int) Math.round(C[0] / beta);
// System.out.println("Loss = " + result.getValue());
// System.out.println("M: " + M + ", beta: " + beta + ", gamma: " + gamma);

// double trueLambda = 0.01;
// double trueMu = 0.015;
// double trueGamma = trueLambda / trueMu;
// double trueBeta = Math.exp((trueLambda - trueMu)*5);
// System.out.println("Loss at true parameters: " + loss.value(new double[]
// {trueBeta, trueGamma}));
// } catch (Exception e) {

// M = null;
// beta = Double.NaN;
// gamma = Double.NaN;
// }
// }

// public double getGamma() {
// return gamma;
// }

// public double getBeta() {
// return beta;
// }

// public Integer getM() {
// return M;
// }
// }