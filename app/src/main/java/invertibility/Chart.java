package invertibility;

import java.awt.Color;
import java.awt.Shape;
import java.awt.geom.Ellipse2D;
import java.util.ArrayList;
import java.util.List;

import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartPanel;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.chart.plot.ValueMarker;
import org.jfree.chart.plot.XYPlot;
import org.jfree.chart.renderer.xy.XYLineAndShapeRenderer;
import org.jfree.chart.ui.ApplicationFrame;
import org.jfree.data.statistics.BoxAndWhiskerCategoryDataset;
import org.jfree.data.statistics.DefaultBoxAndWhiskerCategoryDataset;
import org.jfree.data.xy.XYDataset;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;

/** Class used to chart estimated data vs. actual values. */
public class Chart {

	@SuppressWarnings("serial")
	static class Chart_AWT extends ApplicationFrame {

		public Chart_AWT(String applicationTitle, String chartTitle, String var, double[] estimated, double actual,
				int[] N) {
			super(applicationTitle);
			// JFreeChart chart = ChartFactory.createScatterPlot(chartTitle, "log N", var,
			// createDataset(estimated, N),
			// PlotOrientation.VERTICAL, true, true, false);
			JFreeChart chart = ChartFactory.createBoxAndWhiskerChart(chartTitle, "log N", var,
					createDataset(estimated, N),
					true);

			// XYPlot plot = (XYPlot) chart.getPlot();
			// XYLineAndShapeRenderer renderer = new XYLineAndShapeRenderer();

			// renderer.setSeriesLinesVisible(0, false); // remove lines for estimated data
			// final double CIRCLE_RADIUS = 2.5;
			// Shape circle = new Ellipse2D.Double(-CIRCLE_RADIUS, -CIRCLE_RADIUS, 2 * CIRCLE_RADIUS, 2 * CIRCLE_RADIUS);
			// renderer.setSeriesShape(0, circle); // set estimated data shape to be circle (default is square)

			// renderer.setSeriesShapesVisible(1, false); // set actual data shape to be a line
			// renderer.setSeriesPaint(1, Color.BLUE);

			// plot.setRenderer(renderer);

			// plot.setBackgroundPaint(Color.WHITE);

			// ValueMarker marker = new ValueMarker(actual);
			// marker.setPaint(Color.BLUE);

			// plot.addRangeMarker(marker);

			ChartPanel chartPanel = new ChartPanel(chart);
			chartPanel.setPreferredSize(new java.awt.Dimension(480, 367));

			setContentPane(chartPanel);
		}

		private XYDataset createDataset2(double[] estimated, int[] N) {
			XYSeriesCollection dataset = new XYSeriesCollection();
			XYSeries estimatedSeries = new XYSeries("estimated");

			for (int i = 0; i < N.length; i++) {
				estimatedSeries.add(N[i], estimated[i]);
			}

			dataset.addSeries(estimatedSeries);
			dataset.addSeries(new XYSeries("actual"));
			return dataset;
		}

		private BoxAndWhiskerCategoryDataset createDataset(double[] estimated, int[] N) {
			final int entityCount = 22;

			final DefaultBoxAndWhiskerCategoryDataset dataset = new DefaultBoxAndWhiskerCategoryDataset();

			for (int i = 0; i < N.length; i++) {
				final List<Integer> list = new ArrayList<Integer>();
				// add some values...
				for (int k = 0; k < entityCount; k++) {
					final int value1 = 10 + (int) (Math.random() * 3);
					list.add(value1);
					final int value2 = 15;
					list.add(value2);
				}

				dataset.add(list, "", " Type " + i);
				dataset.add(new ArrayList<Integer>(), "", "Actual");
			}

			return dataset;
		}

	}

	public static void chart(String chartTitle, String var, double[] estimated, double actual, int[] N) {
		Chart_AWT chart = new Chart_AWT(chartTitle, chartTitle, var, estimated, actual, N);
		chart.pack();
		chart.setVisible(true);
	}
}
