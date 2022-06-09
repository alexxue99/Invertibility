package invertibility;

import java.awt.Color;
import java.util.stream.Collectors;
import java.util.stream.DoubleStream;

import org.jfree.chart.ChartFactory;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.chart.plot.ValueMarker;
import org.jfree.chart.plot.XYPlot;
import org.jfree.chart.renderer.xy.XYErrorRenderer;
import org.jfree.chart.swing.ApplicationFrame;
import org.jfree.chart.swing.ChartPanel;
import org.jfree.data.statistics.BoxAndWhiskerCalculator;
import org.jfree.data.statistics.DefaultBoxAndWhiskerCategoryDataset;
import org.jfree.data.xy.XYIntervalDataItem;
import org.jfree.data.xy.XYIntervalSeries;
import org.jfree.data.xy.XYIntervalSeriesCollection;

/** Class used to chart estimated data vs. actual values. */
public class Chart {

	static class Chart_AWT extends ApplicationFrame {

	//	// box and whisker chart
	// 	public Chart_AWT(String applicationTitle, String chartTitle, String var, double[][] estimated, double actual,
	// 			int[] N) {
	// 		super(applicationTitle);

	// 		DefaultBoxAndWhiskerCategoryDataset<String, String> data = new DefaultBoxAndWhiskerCategoryDataset<String, String>();

	// 		for (int i = 0; i < N.length; i++) {
	// 			data.add(
	// 					BoxAndWhiskerCalculator.calculateBoxAndWhiskerStatistics(
	// 							DoubleStream.of(estimated[i]).boxed().collect(Collectors.toList())),
	// 					String.valueOf(N[i]), var);
	// 		}

	// 		JFreeChart chart = ChartFactory.createBoxAndWhiskerChart(
	// 				chartTitle, "log N", var, data, false);

	// 		ChartPanel chartPanel = new ChartPanel(chart);
	// 		chartPanel.setPreferredSize(new java.awt.Dimension(480, 367));

	// 		setContentPane(chartPanel);
	// 	}

		// error bar chart
		public Chart_AWT(String applicationTitle, String chartTitle, String var, double[] estimated,
				double[] estimatedUpper, double[] estimatedLower, double actual,
				int[] N) {
			super(applicationTitle);

			XYIntervalSeriesCollection<String> valuesDataSet = new XYIntervalSeriesCollection<String>();
			XYIntervalSeries<String> values = new XYIntervalSeries<String>("values");
			System.out.println("------\n" + var);
			for (int i = 0; i < N.length; i++) {
				XYIntervalDataItem val;
				System.out.println(estimated[i] + " " + estimatedUpper[i] + " " + estimatedLower[i]);
				if (estimatedUpper[i] > estimatedLower[i])
					val = new XYIntervalDataItem(N[i], N[i], N[i], estimated[i], estimatedLower[i],
							estimatedUpper[i]);
				else
					val = new XYIntervalDataItem(N[i], N[i], N[i], estimated[i], estimatedUpper[i],
							estimatedLower[i]);
				values.add(val, true);
			}

			valuesDataSet.addSeries(values);

			XYIntervalSeries<String> actualValue = new XYIntervalSeries<String>("actual");
			valuesDataSet.addSeries(actualValue);

			JFreeChart chart = ChartFactory.createScatterPlot(chartTitle, "log N", var, valuesDataSet,
					PlotOrientation.VERTICAL,
					true, false, false);

			XYPlot<String> plot = (XYPlot<String>) chart.getPlot();
			XYErrorRenderer errorRenderer = new XYErrorRenderer();
			plot.setRenderer(errorRenderer);

			ValueMarker marker = new ValueMarker(actual);
			marker.setPaint(Color.BLUE);
			plot.addRangeMarker(marker);
			errorRenderer.setSeriesPaint(1, Color.BLUE);
			errorRenderer.setSeriesLinesVisible(1, true);
			errorRenderer.setSeriesShapesVisible(1, false);

			plot.setBackgroundPaint(Color.WHITE);

			ChartPanel chartPanel = new ChartPanel(chart);
			chartPanel.setPreferredSize(new java.awt.Dimension(480, 367));

			setContentPane(chartPanel);
		}
	}

	public static void ErrorChart(String chartTitle, String var, double[] estimated, double[] estimatedUpper,
			double[] estimatedLower,
			double actual, int[] N) {
		Chart_AWT chart = new Chart_AWT(chartTitle, chartTitle, var, estimated, estimatedUpper, estimatedLower, actual,
				N);
		chart.pack();
		chart.setVisible(true);
	}

	// public static void BoxAndWhiskerChart(String chartTitle, String var, double[][] estimated,
	// 		double actual, int[] N) {
	// 	Chart_AWT chart = new Chart_AWT(chartTitle, chartTitle, var, estimated, actual,
	// 			N);
	// 	chart.pack();
	// 	chart.setVisible(true);
	// }

}
