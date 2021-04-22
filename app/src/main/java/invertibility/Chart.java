package invertibility;

import java.awt.Shape;
import java.awt.geom.Ellipse2D;

import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartPanel;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.chart.plot.XYPlot;
import org.jfree.chart.renderer.xy.XYLineAndShapeRenderer;
import org.jfree.chart.ui.ApplicationFrame;
import org.jfree.data.xy.XYDataset;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;

public class Chart {

	@SuppressWarnings("serial")
	static class Chart_AWT extends ApplicationFrame {

		public Chart_AWT(String applicationTitle, String chartTitle, String var, double[] estimated, double actual,
				int[] N) {
			super(applicationTitle);
			JFreeChart chart = ChartFactory.createScatterPlot(chartTitle, "N", var, createDataset(estimated, actual, N),
					PlotOrientation.VERTICAL, true, true, false);

			XYPlot plot = (XYPlot) chart.getPlot();
			XYLineAndShapeRenderer renderer = new XYLineAndShapeRenderer();

			renderer.setSeriesLinesVisible(0, false); // remove lines for estimated data
			final double CIRCLE_RADIUS = 2.5;
			Shape circle = new Ellipse2D.Double(-CIRCLE_RADIUS, -CIRCLE_RADIUS, 2 * CIRCLE_RADIUS, 2 * CIRCLE_RADIUS);
			renderer.setSeriesShape(0, circle); // set estimated data shape to be circle (default is square)

			renderer.setSeriesShapesVisible(1, false); // remove actual data shapes, so only a line appears for the
														// actual data

			plot.setRenderer(renderer);

			ChartPanel chartPanel = new ChartPanel(chart);
			chartPanel.setPreferredSize(new java.awt.Dimension(560, 367));
			setContentPane(chartPanel);
		}

		private XYDataset createDataset(double[] estimated, double actual, int[] N) {
			XYSeriesCollection dataset = new XYSeriesCollection();
			XYSeries estimatedSeries = new XYSeries("estimated");
			XYSeries actualSeries = new XYSeries("actual");

			for (int i = 0; i < N.length; i++) {
				estimatedSeries.add(N[i], estimated[i]);
			}

			actualSeries.add(N[0], actual);
			actualSeries.add(N[N.length - 1], actual);

			dataset.addSeries(estimatedSeries);
			dataset.addSeries(actualSeries);
			return dataset;
		}
	}

	public static void chart(String chartTitle, String var, double[] estimated, double actual, int[] N) {
		Chart_AWT chart = new Chart_AWT(chartTitle, chartTitle, var, estimated, actual, N);
		chart.pack();
		chart.setVisible(true);
	}
}
