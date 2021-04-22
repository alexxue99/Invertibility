# Invertibility

This application partially inverts the length process described in the process for ultrametric trees. It estimates the length of the root sequence M, the ratio gamma = lambda / mu, and the exponential beta = exp(-(mu - lambda)t).

The class Main.java contains two examples on how to use the application. Example 1 is a basic example on how to input a tree, invert it, and print out the estimated parameters. Example 2 shows how to run trials on a variable number of leaves, and then to plot the results on charts. The plotting uses JFreeChart 1.5.2, a free library that can be found at https://www.jfree.org/index.html.

This project uses Gradle to handle the JFreeChart dependency. If the user has no need to plot data on multiple trials, then the user does not need to import the whole project. The user will only need the files Chart.java, Invert.java, and TreeSimul.java, all found at \app\src\main\java\invertibility\. If the user wants to view the charts and already has JFreeChart installed, then the user can handle the dependency themselves. The user will only need the file Main.java, in addition to the three files mentioned above.

With Gradle, the following two commands may be used to compile Main.java and create a .jar file:<br/>
"gradle build"<br/>
"gradle jar"

The created .jar file is \app\build\libs\app.jar, which can be run using<br/>
"java -jar \app\build\libs\app.jar"<br/>
from the project directory.
