package invertibility;

public class InvertPairwiseDistance extends Invert {
    private TreeLeaves tree2; // star tree for second set of leaves
    private double lambda2; // lambda*t for second tree
    private double mu2; // mu*t for second tree

    private double pwd; // pairwise distance
    private double wd; // distance for the common ancestor w

    /**
     * Inverts the process based on given data on the lengths of the sequences at
     * the leaves to estimate the pairwise distances between leaves.
     */
    public InvertPairwiseDistance(double lambda, double mu, double lambda2, double mu2, int M, TreeLeaves tree,
            TreeLeaves tree2) {
        super(tree);

        this.tree2 = tree2;
        this.lambda = lambda;
        this.mu = mu;
        this.lambda2 = lambda2;
        this.mu2 = mu2;
        this.M = M;

        estimatePairwiseDistance();
        estimateAncestorDistance();
    }

    protected void calcPartials() {
    }

    private double covariance() {
        double cov = 0;
        int count = 0;

        for (int l1 : tree.getSeqLeavesLengths()) {
            for (int l2 : tree2.getSeqLeavesLengths()) {
                cov += (l1 * l2 - cov) / (++count);
            }
        }

        double mean1 = M * Math.exp(lambda - mu);
        double mean2 = M * Math.exp(lambda2 - mu2);

        cov -= mean1 * mean2;
        return cov;
    }

    // estimates (mu - lambda) t_uv
    private void estimatePairwiseDistance() {
        double cov = covariance();

        double exp = (mu - lambda) / M / (lambda + mu) * Math.exp((mu - lambda + mu2 - lambda2) / 2) * cov
                + Math.exp(-(mu - lambda + mu2 - lambda2) / 2);

        pwd = Math.log(exp) * -2;
        System.out.println(pwd);
    }

    // estimates (mu - lambda)t_w
    private void estimateAncestorDistance() {
        wd = (mu - lambda + mu2 - lambda2 - pwd) / 2;
    }

    public double getPairwiseDistance() {
        return pwd;
    }

    public double getAncestorDistance() {
        return wd;
    }
}
