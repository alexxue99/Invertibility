package invertibility;

public class Invert1Mer extends Invert {
    private double[] D;

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
     * Inverts the process based on given data on the sequences at
     * the leaves, as well as on mu, lambda, and M, to estimate nu, pi, and a.
     */
    public Invert1Mer(double lambda, double mu, double pi0, int M, TreeLeaves tree) {
        super(tree);

        D = new double[3];

        this.lambda = lambda;
        this.mu = mu;
        expmu = Math.exp(mu);
        this.pi0 = pi0;
        pi1 = 1 - pi0;
        this.M = M;

        calcPartials();
        updateCs(); // maybe not needed, move updateCs in invert to invertlength?
        // updateDs(); only needed for other way of doing 1mer?

        estimateNu();
        estimatepi0();
        estimateA();
    }

    protected void calcPartials() {
        int[] ones = tree.getSeqNumOnes();
        int[] zeros = tree.getSeqNumZeros();

        partials_1 = factorialMoments(ones);
        partials_2 = factorialMoments(zeros);

        // System.out.println("Partials ones " + partials_1[0] + " " + partials_1[1] + "
        // " + partials_1[2]);
        // System.out.println("Partials zeros " + partials_2[0] + " " + partials_2[1] +
        // " " + partials_2[2]);

        partial_12 = 0;
        partial_112 = 0;
        partial_122 = 0;

        for (int i = 0; i < N; i++) {
            partial_12 += (ones[i] * zeros[i] - partial_12) / (i + 1);
            partial_112 += (ones[i] * (ones[i] - 1) * zeros[i] - partial_112) / (i + 1);
            partial_122 += (ones[i] * zeros[i] * (zeros[i] - 1) - partial_122) / (i + 1);
        }

        // System.out.println("Partial_12: " + partial_12);
        // System.out.println("Partial_112: " + partial_112);
        // System.out.println("Partial_122: " + partial_122);

        partials[0] = partials_1[0] * pi0 - partials_2[0] * pi1;
        partials[1] = partials_1[1] * pi0 * pi0 + partials_2[1] * pi1 * pi1 - 2 * partial_12 * pi0 * pi1;
        partials[2] = partials_1[2] * Math.pow(pi0, 3) - partials_2[2] * Math.pow(pi1, 3)
                - 3 * partial_112 * pi0 * pi0 * pi1 + 3 * partial_122 * pi0 * pi1 * pi1;

        // System.out.println("PARTIALS: " + partials[0] + " " + partials[1] + " " +
        // partials[2]);
    }

    /* Prereq: Cs are calculated */
    private void updateDs() {
        D[0] = C[0] * expmu;
        D[1] = -C[1] * expmu * expmu;
        D[2] = C[2] * Math.pow(expmu, 3) / 2;
        System.out.println("D: " + D[0] + "\t" + D[1] + "\t" + D[2]);
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
        // discriminant = -3 * D[0] * D[0] * D[1] * D[1] + 4 * Math.pow(D[0], 3) * D[2]
        // + 4 * M * D[1] * D[1]
        // - 6 * M * D[0] * D[1] * D[2] + M * M * D[2] * D[2];
        // exp = Math.sqrt(discriminant);
        // exp /= (D[0] * D[0] - M * D[1]);

        // double add = 2 * (D[0] * D[1] - M * D[2]);
        // System.out.println("ADD: " + add + "\tExp: " + exp);
        // if (add + exp < 0) {
        // exp = add - exp;
        // flip = true;
        // } else {
        // exp = add + exp;
        // }

        // System.out.println((D[0] * D[1] - M * D[2]) / (D[0] * D[0] - M * D[1]));
        // System.out.println("EXP: " + exp + "\t" + discriminant);
        nu = -Math.log(expnegnu);
        // System.out.println("NU: " + nu);
    }

    private void estimatepi0() {
        // pi0 = D[0] * D[1] - M * D[2];
        // pi0 /= Math.sqrt(discriminant);

        // if (flip)
        // pi0 *= -1;

        // pi0 += 0.5;
    }

    private void estimateA() {
        // a = D[0] + M * (1 - pi0);
        // a /= exp;

        a = partials[0] + M * pi1 * expnegnu / expmu;
        a *= expmu / expnegnu;
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
