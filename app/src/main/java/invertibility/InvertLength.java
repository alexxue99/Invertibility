package invertibility;

public class InvertLength extends Invert {
    private double C2prime;
    private double C3prime;

    /**
     * Inverts the process based on given data on the lengths of the sequences at
     * the leaves to estimate gamma, beta, and M.
     */
    public InvertLength(TreeLeaves tree) {
        super(tree);

        calcPartials();
        updateCs();
        updateCprimes();

        estimateGamma();
        estimateBeta();
        estimateM();
    }

    protected void calcPartials() {
        partials = factorialMoments(tree.getSeqLeavesLengths());
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

    public double getGamma() {
        return gamma;
    }

    public double getBeta() {
        return beta;
    }

    public double getM() {
        return M;
    }
}
