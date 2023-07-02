package invertibility;

import java.awt.Color;

import org.jfree.chart.ChartFactory;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.chart.plot.ValueMarker;
import org.jfree.chart.plot.XYPlot;
import org.jfree.chart.renderer.xy.XYErrorRenderer;
import org.jfree.chart.swing.ApplicationFrame;
import org.jfree.chart.swing.ChartPanel;
import org.jfree.data.category.DefaultCategoryDataset;
import org.jfree.data.xy.XYIntervalDataItem;
import org.jfree.data.xy.XYIntervalSeries;
import org.jfree.data.xy.XYIntervalSeriesCollection;

/** Class used to chart estimated data vs. exact values. */
public class Chart {

	static class Chart_AWT extends ApplicationFrame {
		// error bar chart
		public Chart_AWT(String xAxis, String applicationTitle, String chartTitle, String var, double[] estimated,
				double[] deviation, double exact,
				int[] N) {
			super(applicationTitle);

			XYIntervalSeriesCollection<String> valuesDataSet = new XYIntervalSeriesCollection<String>();
			XYIntervalSeries<String> values = new XYIntervalSeries<String>("values");
			System.out.println("------\n" + var);
			for (int i = 0; i < N.length; i++) {
				XYIntervalDataItem val;
				System.out.println(estimated[i] + " " + deviation[i]);
				val = new XYIntervalDataItem(N[i], N[i], N[i], estimated[i], estimated[i] + deviation[i],
						estimated[i] - deviation[i]);
				values.add(val, true);
			}

			valuesDataSet.addSeries(values);

			XYIntervalSeries<String> exactValue = new XYIntervalSeries<String>("");
			valuesDataSet.addSeries(exactValue);

			JFreeChart chart = ChartFactory.createScatterPlot(chartTitle, xAxis, var, valuesDataSet,
					PlotOrientation.VERTICAL,
					true, false, false);

			@SuppressWarnings("unchecked")
			XYPlot<String> plot = (XYPlot<String>) chart.getPlot();
			XYErrorRenderer errorRenderer = new XYErrorRenderer();
			plot.setRenderer(errorRenderer);

			ValueMarker marker = new ValueMarker(exact);
			marker.setPaint(Color.BLUE);
			plot.addRangeMarker(marker);
			errorRenderer.setSeriesPaint(0, Color.RED);
			errorRenderer.setSeriesVisibleInLegend(0, false);

			errorRenderer.setSeriesPaint(1, Color.BLUE);
			errorRenderer.setSeriesLinesVisible(1, true);
			errorRenderer.setSeriesShapesVisible(1, false);
			errorRenderer.setSeriesVisibleInLegend(1, false);

			plot.setBackgroundPaint(Color.WHITE);

			ChartPanel chartPanel = new ChartPanel(chart);
			chartPanel.setPreferredSize(new java.awt.Dimension(540, 360));

			setContentPane(chartPanel);
		}

		// for line chart
		public Chart_AWT(String xAxis, String applicationTitle, String chartTitle, String var, double[] estimated,
				int[] N) {
			super(applicationTitle);

			DefaultCategoryDataset<String, Integer> dataset = new DefaultCategoryDataset<String, Integer>();

			System.out.println("------\n" + var);
			for (int i = 0; i < N.length; i++) {
				dataset.addValue((Number) estimated[i], "difference", N[i]);
			}

			JFreeChart chart = ChartFactory.createLineChart(chartTitle, xAxis, var, dataset,
					PlotOrientation.VERTICAL,
					true, false, false);

			ChartPanel chartPanel = new ChartPanel(chart);
			chartPanel.setPreferredSize(new java.awt.Dimension(560, 367));
			setContentPane(chartPanel);
		}
	}

	public static void ErrorChart(String xAxis, String var, double[] estimated, double[] deviation,
			double exact, int[] N) {
		Chart_AWT chart = new Chart_AWT(xAxis, var + " vs. " + xAxis, "", var, estimated, deviation, exact,
				N);
		chart.pack();
		chart.setVisible(true);
	}

	
	public static void LineChart(String xAxis, String var, double[] estimated, int[] N) {
		Chart_AWT chart = new Chart_AWT(xAxis, var + " vs. " + xAxis, "", var, estimated,
				N);
		chart.pack();
		chart.setVisible(true);
	}
}
