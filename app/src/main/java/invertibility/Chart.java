package invertibility;

import java.awt.Color;
import java.awt.Font;
import java.awt.GraphicsEnvironment;
import java.io.File;
import java.util.ArrayList;
import java.util.List;

import org.jfree.chart.ChartFactory;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.axis.CategoryAxis;
import org.jfree.chart.axis.NumberAxis;
import org.jfree.chart.plot.CategoryPlot;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.chart.plot.ValueMarker;
import org.jfree.chart.plot.XYPlot;
import org.jfree.chart.renderer.category.BoxAndWhiskerRenderer;
import org.jfree.chart.renderer.xy.XYErrorRenderer;
import org.jfree.chart.swing.ApplicationFrame;
import org.jfree.chart.swing.ChartPanel;
import org.jfree.data.statistics.DefaultBoxAndWhiskerCategoryDataset;
import org.jfree.data.xy.XYIntervalDataItem;
import org.jfree.data.xy.XYIntervalSeries;
import org.jfree.data.xy.XYIntervalSeriesCollection;

/** Class used to chart estimated data vs. exact values. */
public class Chart {

	static class Chart_AWT extends ApplicationFrame {
		// error bar chart
		public Chart_AWT(String xAxis, String applicationTitle, String chartTitle, String var, double[] estimated,
				double[] deviation, double exact,
				int[] N, boolean plain) {
			super(applicationTitle);

			XYIntervalSeriesCollection<String> valuesDataSet = new XYIntervalSeriesCollection<String>();
			XYIntervalSeries<String> values = new XYIntervalSeries<String>("values");
			System.out.println("------\n" + var);
			for (int i = 0; i < N.length; i++) {
				XYIntervalDataItem val;
				System.out.println(estimated[i] + " " + deviation[i]);
				double lower = (exact == -1) ? estimated[i] : estimated[i] - deviation[i];
				val = new XYIntervalDataItem(N[i], N[i], N[i], estimated[i], estimated[i] + deviation[i],
						lower);
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

			String path = "C:\\Users\\alexx\\workspace\\Invertibility\\cm\\cmunrm.ttf";
			Font customFont = null;
			try {
				customFont = Font.createFont(Font.TRUETYPE_FONT, new File(path));
				GraphicsEnvironment ge = GraphicsEnvironment.getLocalGraphicsEnvironment();
				ge.registerFont(customFont);
				System.out.println("SUCCESS");
			} catch (Exception e) {
				System.out.println("Font error");
			}

			customFont = new Font("CMU Serif", plain ? Font.PLAIN : Font.ITALIC, 20);
			plot.getRangeAxis().setLabelFont(customFont);

			customFont = new Font("CMU Serif", Font.ITALIC, 20);
			plot.getDomainAxis().setLabelFont(customFont);

			ChartPanel chartPanel = new ChartPanel(chart);
			chartPanel.setPreferredSize(new java.awt.Dimension(540, 360));

			setContentPane(chartPanel);
		}

		// box and whisker plot
		public Chart_AWT(String xAxis, String applicationTitle, String chartTitle, String var, double[][] estimated,
				double exact,
				int[] N, boolean plain) {
			super(applicationTitle);

			DefaultBoxAndWhiskerCategoryDataset<String, Integer> dataset = new DefaultBoxAndWhiskerCategoryDataset<String, Integer>();
			for (int i = 0; i < N.length; i++) {
				List<Double> list = new ArrayList<Double>();
				for (double d : estimated[i])
					list.add(d);
				dataset.add(list, "", N[i]);
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
			renderer.setMeanVisible(false);
			//renderer.setMaxOutlierVisible(false);
			//renderer.setMinOutlierVisible(false);
			renderer.setMaximumBarWidth(0.10);

			String path = "C:\\Users\\alexx\\workspace\\Invertibility\\cm\\cmunrm.ttf";
			Font customFont = null;
			try {
				customFont = Font.createFont(Font.TRUETYPE_FONT, new File(path));
				GraphicsEnvironment ge = GraphicsEnvironment.getLocalGraphicsEnvironment();
				ge.registerFont(customFont);
				System.out.println("SUCCESS");
			} catch (Exception e) {
				System.out.println("Font error");
			}

			customFont = new Font("CMU Serif", plain ? Font.PLAIN : Font.ITALIC, 20);
			plot.getRangeAxis().setLabelFont(customFont);

			customFont = new Font("CMU Serif", Font.PLAIN, 20);
			plot.getDomainAxis().setLabelFont(customFont);

			JFreeChart chart = new JFreeChart("", customFont, plot, true);
			chart.setBackgroundPaint(Color.WHITE);
			ChartPanel chartPanel = new ChartPanel(chart);
			chartPanel.setPreferredSize(new java.awt.Dimension(540, 360));

			setContentPane(chartPanel);
		}
	}

	public static void ErrorChart(String xAxis, String var, double[] estimated, double[] deviation,
			double exact, int[] N) {
		boolean plain = false;
		switch (var) {
			case "gamma":
				var = "\u03B3";
				break;
			case "beta":
				var = "\u03B2";
				break;
			case "nu":
				var = "\u03BD";
				break;
			case "difference":
				plain = true;
		}

		Chart_AWT chart = new Chart_AWT(xAxis, var + " vs. " + xAxis, "", var,
				estimated, deviation, exact,
				N, plain);
		chart.pack();
		chart.setVisible(true);
	}

	public static void BoxWhiskerChart(String xAxis, String var, double[][] estimated,
			double exact, int[] N) {
		boolean plain = false;
		switch (var) {
			case "gamma":
				var = "\u03B3";
				break;
			case "beta":
				var = "\u03B2";
				break;
			case "nu":
				var = "\u03BD";
				break;
			case "difference":
				plain = true;
		}

		Chart_AWT chart = new Chart_AWT(xAxis, var + " vs. " + xAxis, "", var,
				estimated, exact,
				N, plain);
		chart.pack();
		chart.setVisible(true);
	}
}
