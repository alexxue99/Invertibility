package invertibility;

import java.util.List;

public class InvertPairwiseDistance extends Invert {
    private double lambda2; // lambda*t for second tree
    private double mu2; // mu*t for second tree

    private int[] products; // products of lengths of leaves, computed through TreeSimul

    private double pwd; // pairwise distance
    private double wd; // distance for the common ancestor w

    /**
     * Inverts the process based on given data on the lengths of the sequences at
     * the leaves to estimate the pairwise distances between leaves.
     */
    public InvertPairwiseDistance(double lambda, double mu, double tu, double tv, int M, int[] products) {
        this.lambda = lambda * tu;
        this.mu = mu * tu;
        this.lambda2 = lambda * tv;
        this.mu2 = mu * tv;
        this.M = M;

        this.products = products;

        estimatePairwiseDistance();
        estimateAncestorDistance();
    }

    protected void calcPartials() {
    }

    private double covariance() {
        double cov = 0;
        int count = 0;

        for (int p : products) {
            cov += (p - cov) / (++count);
        }

        double mean1 = M * Math.exp(lambda - mu);
        double mean2 = M * Math.exp(lambda2 - mu2);

        cov -= mean1 * mean2;
        return cov;
    }

    // estimates mu * t_uv
    private void estimatePairwiseDistance() {
        double cov = covariance();

        double exp = (mu - lambda) / M / (lambda + mu) * Math.exp((mu - lambda + mu2 - lambda2) / 2) * cov
                + Math.exp(-(mu - lambda + mu2 - lambda2) / 2);

        pwd = Math.log(exp) * -2;
        pwd *= mu / (mu - lambda);
        System.out.println(pwd);
    }

    // estimates mu * t_w
    private void estimateAncestorDistance() {
        wd = (mu + mu2 - pwd) / 2;
    }

    public double getPairwiseDistance() {
        return pwd;
    }

    public double getAncestorDistance() {
        return wd;
    }
}
