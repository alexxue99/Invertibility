package invertibility;

/**
 * Class to invert the 1mer process based on given data on the sequences at the
 * leaves, as well as on mu, lambda, and M, to estimate nu and a.
 */
public class Invert1Mer extends Invert {
    private double pi0;
    private double pi1;

    private double[] partials_1;
    private double[] partials_2;

    private double partial_12;
    private double partial_112;
    private double partial_122;

    private double expmu;
    private double expnegnu;

    /**
     * Constructor to initialize an Invert1Mer object.
     * 
     * @param lambda the insertion rate
     * @param mu     the deletion rate
     * @param pi0    the probability that a new character is 0 after substitution or
     *               insertion
     * @param M      the length of the root sequence
     * @param tree   a TreeLeaves object containing the sequences at the leaves of
     *               the tree
     */
    public Invert1Mer(double lambda, double mu, double pi0, int M, TreeLeaves tree) {
        super(tree);

        this.lambda = lambda;
        this.mu = mu;
        expmu = Math.exp(mu);
        this.pi0 = pi0;
        pi1 = 1 - pi0;
        this.M = M;

        calcPartials();
        estimateNu();
        estimateA();
    }

    protected void calcPartials() {
        int[] ones = tree.getSeqNumOnes();
        int[] zeros = tree.getSeqNumZeros();

        partials_1 = factorialMoments(ones);
        partials_2 = factorialMoments(zeros);

        partial_12 = 0;
        partial_112 = 0;
        partial_122 = 0;

        for (int i = 0; i < N; i++) {
            partial_12 += (ones[i] * zeros[i] - partial_12) / (i + 1);
            partial_112 += (ones[i] * (ones[i] - 1) * zeros[i] - partial_112) / (i + 1);
            partial_122 += (ones[i] * zeros[i] * (zeros[i] - 1) - partial_122) / (i + 1);
        }

        partials[0] = partials_1[0] * pi0 - partials_2[0] * pi1;
        partials[1] = partials_1[1] * pi0 * pi0 + partials_2[1] * pi1 * pi1 - 2 * partial_12 * pi0 * pi1;
        partials[2] = partials_1[2] * Math.pow(pi0, 3) - partials_2[2] * Math.pow(pi1, 3)
                - 3 * partial_112 * pi0 * pi0 * pi1 + 3 * partial_122 * pi0 * pi1 * pi1;
    }

    private void estimateNu() {
        double expmu = Math.exp(mu);
        double term = pi0 * partials_1[0] - pi1 * partials_2[0];
        double A = -M * pi1 * pi1 + M * pi1 * (pi1 * pi1 - pi0 * pi0);
        double B = expmu * (pi1 * pi1 - pi0 * pi0) * term;
        double C = -expmu * expmu
                * (pi0 * pi0 * partials_1[1] - 2 * pi0 * pi1 * partial_12 + pi1 * pi1 * partials_2[1] - term * term);

        double discriminant = B * B - 4 * A * C;
        expnegnu = (-B - Math.sqrt(discriminant)) / (2 * A);
        nu = -Math.log(expnegnu);
    }

    private void estimateA() {
        a = partials[0] * expmu + M * pi1 * expnegnu;
        a /= expnegnu;
    }

    public double getNu() {
        return nu;
    }

    public double getpi0() {
        return pi0;
    }

    public double getA() {
        return a;
    }
}
