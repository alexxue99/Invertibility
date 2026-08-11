# Invertibility

This application inverts the TKF91 process described in the paper available at [doi link](http://dx.doi.org/10.2139/ssrn.5020031). Throughout, we assume that the height of leaves are 1 for simplicity. Four classes are provided to invert various parameters of the process. InvertLength inverts the length process, estimating gamma, beta, and M. Invert1Mer inverts the 1mer process, estimating nu and a. InvertState partially inverts the TKF91 process, estimating the root sequence. Lastly, InvertPairwiseDistance partially inverts the TKF91 process based on the covariance of two nodes u and v, estimating the distance between u and v.

The class Main.java contains code that tests these four classes. There are four methods provided, each testing a class by creating trees for various values of N for some number of trials. We've provided (commented out) code that may be used to plot the results in charts. The plotting uses JFreeChart 2.0.0, a free library that can be found at https://www.jfree.org/index.html.

This project uses Gradle to handle the JFreeChart dependency. The following Gradle command may be used to compile Main.java and then run the created .jar file:<br/>
"gradlew.bat build runJar"<br/>

The data can be found in [app/data/](app/data/), and the provided MATLAB files [graph_results.m](app/src/graph_results.m) and [graph_all_results.m](app/src/graph_all_results.m) reproduce the boxplots found in the paper.
