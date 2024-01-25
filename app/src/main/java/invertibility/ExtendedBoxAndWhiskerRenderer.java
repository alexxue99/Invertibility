// Source code is decompiled from a .class file using FernFlower decompiler.
package invertibility;

import java.awt.Graphics2D;
import java.awt.Paint;
import java.awt.Shape;
import java.awt.Stroke;
import java.awt.geom.Ellipse2D;
import java.awt.geom.Line2D;
import java.awt.geom.Rectangle2D;
import java.io.Serializable;

import org.jfree.chart.api.PublicCloneable;
import org.jfree.chart.api.RectangleEdge;
import org.jfree.chart.axis.CategoryAxis;
import org.jfree.chart.axis.ValueAxis;
import org.jfree.chart.entity.EntityCollection;
import org.jfree.chart.internal.PaintUtils;
import org.jfree.chart.plot.CategoryPlot;
import org.jfree.chart.renderer.category.BoxAndWhiskerRenderer;
import org.jfree.chart.renderer.category.CategoryItemRendererState;
import org.jfree.data.category.CategoryDataset;
import org.jfree.data.statistics.BoxAndWhiskerCategoryDataset;

public class ExtendedBoxAndWhiskerRenderer extends BoxAndWhiskerRenderer implements Cloneable, PublicCloneable, Serializable {
   public ExtendedBoxAndWhiskerRenderer() {
      super();
   }

   public void drawVerticalItem(Graphics2D g2, CategoryItemRendererState state, Rectangle2D dataArea, CategoryPlot plot, CategoryAxis domainAxis, ValueAxis rangeAxis, CategoryDataset dataset, int row, int column) {
      BoxAndWhiskerCategoryDataset bawDataset = (BoxAndWhiskerCategoryDataset)dataset;
      double categoryEnd = domainAxis.getCategoryEnd(column, this.getColumnCount(), dataArea, plot.getDomainAxisEdge());
      double categoryStart = domainAxis.getCategoryStart(column, this.getColumnCount(), dataArea, plot.getDomainAxisEdge());
      double categoryWidth = categoryEnd - categoryStart;
      int seriesCount = this.getRowCount();
      int categoryCount = this.getColumnCount();
      double xx;
      double yyAverage;
      double yyOutlier;
      if (seriesCount > 1) {
         yyAverage = dataArea.getWidth() * this.getItemMargin() / (double)(categoryCount * (seriesCount - 1));
         yyOutlier = state.getBarWidth() * (double)seriesCount + yyAverage * (double)(seriesCount - 1);
         double offset = (categoryWidth - yyOutlier) / 2.0;
         xx = categoryStart + offset + (double)row * (state.getBarWidth() + yyAverage);
      } else {
         yyAverage = (categoryWidth - state.getBarWidth()) / 2.0;
         xx = categoryStart + yyAverage;
      }

      Paint itemPaint = this.getItemPaint(row, column);
      g2.setPaint(itemPaint);
      Stroke s = this.getItemStroke(row, column);
      g2.setStroke(s);
      double aRadius = 0.0;
      RectangleEdge location = plot.getRangeAxisEdge();
      Number yQ1 = bawDataset.getQ1Value(row, column);
      Number yQ3 = bawDataset.getQ3Value(row, column);
      Number yMax = bawDataset.getMaxRegularValue(row, column);
      Number yMin = bawDataset.getMinRegularValue(row, column);
      Shape box = null;
      double maxAxisValue;
      double minAxisValue;
      double oRadius;
      if (yQ1 != null && yQ3 != null && yMax != null && yMin != null) {
         maxAxisValue = rangeAxis.valueToJava2D(yQ1.doubleValue(), dataArea, location);
         minAxisValue = rangeAxis.valueToJava2D(yQ3.doubleValue(), dataArea, location);
         oRadius = rangeAxis.valueToJava2D(yMax.doubleValue(), dataArea, location);
         double yyMin = rangeAxis.valueToJava2D(yMin.doubleValue(), dataArea, location);
         double xxmid = xx + state.getBarWidth() / 2.0;
         double halfW = state.getBarWidth() / 2.0 * this.getWhiskerWidth();
         box = new Rectangle2D.Double(xx, Math.min(maxAxisValue, minAxisValue), state.getBarWidth(), Math.abs(maxAxisValue - minAxisValue));
         if (this.getFillBox()) {
            g2.fill(box);
         }

         Paint outlinePaint = this.getItemOutlinePaint(row, column);
         if (this.getUseOutlinePaintForWhiskers()) {
            g2.setPaint(outlinePaint);
         }

         g2.draw(new Line2D.Double(xxmid, oRadius, xxmid, minAxisValue));
         g2.draw(new Line2D.Double(xxmid - halfW, oRadius, xxmid + halfW, oRadius));
         g2.draw(new Line2D.Double(xxmid, yyMin, xxmid, maxAxisValue));
         g2.draw(new Line2D.Double(xxmid - halfW, yyMin, xxmid + halfW, yyMin));
         g2.setStroke(this.getItemOutlineStroke(row, column));
         g2.setPaint(outlinePaint);
         g2.draw(box);
      }

      g2.setPaint(this.getArtifactPaint());
      Number yMedian;
      if (this.isMeanVisible()) {
         yMedian = bawDataset.getMeanValue(row, column);
         if (yMedian != null) {
            yyAverage = rangeAxis.valueToJava2D(yMedian.doubleValue(), dataArea, location);
            aRadius = state.getBarWidth() / 4.0;
            if (yyAverage > dataArea.getMinY() - aRadius && yyAverage < dataArea.getMaxY() + aRadius) {
               Ellipse2D.Double avgEllipse = new Ellipse2D.Double(xx + aRadius, yyAverage - aRadius, aRadius * 2.0, aRadius * 2.0);
               g2.fill(avgEllipse);
               g2.draw(avgEllipse);
            }
         }
      }

      if (this.isMedianVisible()) {
         yMedian = bawDataset.getMedianValue(row, column);
         if (yMedian != null) {
            double yyMedian = rangeAxis.valueToJava2D(yMedian.doubleValue(), dataArea, location);
            g2.draw(new Line2D.Double(xx, yyMedian, xx + state.getBarWidth(), yyMedian));
         }
      }

      maxAxisValue = rangeAxis.valueToJava2D(rangeAxis.getUpperBound(), dataArea, location) + aRadius;
      minAxisValue = rangeAxis.valueToJava2D(rangeAxis.getLowerBound(), dataArea, location) - aRadius;
      g2.setPaint(itemPaint);
      oRadius = state.getBarWidth() / 3.0;
      
      if (state.getInfo() != null && box != null) {
         EntityCollection entities = state.getEntityCollection();
         if (entities != null) {
            this.addItemEntity(entities, dataset, row, column, box);
         }
      }

   }
   
   public boolean equals(Object obj) {
      if (obj == this) {
         return true;
      } else if (!(obj instanceof ExtendedBoxAndWhiskerRenderer)) {
         return false;
      } else {
         ExtendedBoxAndWhiskerRenderer that = (ExtendedBoxAndWhiskerRenderer)obj;
         if (this.getFillBox() != that.getFillBox()) {
            return false;
         } else if (this.getItemMargin() != that.getItemMargin()) {
            return false;
         } else if (this.getMaximumBarWidth() != that.getMaximumBarWidth()) {
            return false;
         } else if (this.isMeanVisible()!= that.isMeanVisible()) {
            return false;
         } else if (this.isMedianVisible() != that.isMedianVisible()) {
            return false;
         } else if (this.isMinOutlierVisible() != that.isMinOutlierVisible()) {
            return false;
         } else if (this.isMaxOutlierVisible() != that.isMaxOutlierVisible()) {
            return false;
         } else if (this.getUseOutlinePaintForWhiskers() != that.getUseOutlinePaintForWhiskers()) {
            return false;
         } else if (this.getWhiskerWidth() != that.getWhiskerWidth()) {
            return false;
         } else {
            return !PaintUtils.equal(this.getArtifactPaint(), that.getArtifactPaint()) ? false : super.equals(obj);
         }
      }
   }
}
