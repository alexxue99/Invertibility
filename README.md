# Invertibility

This application partially inverts the length process described in the process for ultrametric trees. It estimates the length of the root sequence M, the ratio gamma = lambda / mu, and the exponential beta = exp(-(mu - lambda)t).

The class Main.java contains two examples on how to use the application. Example 1 is a basic example on how to input the lengths of the sequences at the leaves of a tree, invert the length process, and print out the estimated parameters. Example 2 shows how to run trials on a variable number of leaves by using the TreeSimul class to simulate the length process. It then plots the results in charts. The plotting uses JFreeChart 1.5.2, a free library that can be found at https://www.jfree.org/index.html.

This project uses Gradle to handle the JFreeChart dependency. With Gradle, the following two commands may be used to compile Main.java and run the created .jar file:<br/>
"gradlew.bat build"<br/>
"gradlew.bat runJar"
