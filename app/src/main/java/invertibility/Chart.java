package invertibility;

import java.awt.Color;
import java.util.ArrayList;
import java.util.List;

import org.jfree.chart.JFreeChart;
import org.jfree.chart.axis.CategoryAxis;
import org.jfree.chart.axis.NumberAxis;
import org.jfree.chart.plot.CategoryPlot;
import org.jfree.chart.plot.ValueMarker;
import org.jfree.chart.swing.ApplicationFrame;
import org.jfree.chart.swing.ChartPanel;
import org.jfree.data.statistics.DefaultBoxAndWhiskerCategoryDataset;

/**
 * Class used to plot estimated data vs. exact values in a box and whisker plot.
 */
public class Chart {

	static class Chart_AWT extends ApplicationFrame {
		public Chart_AWT(String xAxis, String applicationTitle, String chartTitle, String var, double[][] estimated,
				double exact,
				int[] N, Boolean meanVisible) {
			super(applicationTitle);

			DefaultBoxAndWhiskerCategoryDataset<String, Integer> dataset = new DefaultBoxAndWhiskerCategoryDataset<String, Integer>();
			for (int i = 0; i < N.length; i++) {
				List<Double> list = new ArrayList<Double>();
				for (double d : estimated[i])
					list.add(d);
				dataset.add(list, "", (int) Math.round(Math.log10(N[i])));
			}

			CategoryAxis x = new CategoryAxis(xAxis);
			NumberAxis y = new NumberAxis(var);
			y.setAutoRangeIncludesZero(false);
			ExtendedBoxAndWhiskerRenderer renderer = new ExtendedBoxAndWhiskerRenderer();

			CategoryPlot<String, Integer> plot = new CategoryPlot<String, Integer>(dataset, x, y, renderer);

			ValueMarker marker = new ValueMarker(exact);
			marker.setPaint(Color.BLUE);
			plot.addRangeMarker(marker);

			plot.setBackgroundPaint(Color.WHITE);

			renderer.setSeriesPaint(0, Color.RED);
			renderer.setSeriesVisibleInLegend(0, false);

			renderer.setFillBox(false);
			renderer.setMeanVisible(meanVisible);

			renderer.setMaximumBarWidth(0.10);
			renderer.setWhiskerWidth(0.6);

			JFreeChart chart = new JFreeChart("", plot);
			chart.setBackgroundPaint(Color.WHITE);
			ChartPanel chartPanel = new ChartPanel(chart);
			chartPanel.setPreferredSize(new java.awt.Dimension(540, 360));

			setContentPane(chartPanel);
		}
	}

	public static void BoxWhiskerChart(String xAxis, String var, double[][] estimated,
			double exact, int[] N, Boolean meanVisible) {
		Chart_AWT chart = new Chart_AWT(xAxis, var + " vs. " + xAxis, "", var,
				estimated, exact,
				N, meanVisible);
		chart.pack();
		chart.setVisible(true);
	}
}
