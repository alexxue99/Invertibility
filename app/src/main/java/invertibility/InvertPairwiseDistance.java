package invertibility;

import java.util.List;

public class InvertPairwiseDistance extends Invert {
    private double lambda2; // lambda*t for second tree
    private double mu2; // mu*t for second tree

    private double covariance;

    private double pwd; // pairwise distance
    private double wd; // distance for the common ancestor w

    /**
     * Inverts the process based on given data on the lengths of the sequences at
     * the leaves to estimate the pairwise distances between leaves.
     */
    public InvertPairwiseDistance(double lambda, double mu, double tu, double tv, int M, double covariance) {
        this.lambda = lambda * tu;
        this.mu = mu * tu;
        this.lambda2 = lambda * tv;
        this.mu2 = mu * tv;
        this.M = M;

        this.covariance = covariance;

        estimatePairwiseDistance();
        estimateAncestorDistance();
    }

    protected void calcPartials() {
    }

    // estimates mu * t_uv
    private void estimatePairwiseDistance() {
        double exp = (mu - lambda) / M / (lambda + mu) * Math.exp((mu - lambda + mu2 - lambda2) / 2) * covariance
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
