package djr.util.gui;

import javax.swing.JFrame;
import cern.jet.stat.*;
import ptolemy.plot.*;
import corejava.*;

import djr.util.array.*;

/**
 * Class <code>MyPlot</code>
 *
 * @author <a href="mailto:reiss@uw.edu">David Reiss</a>
 * @version 1.9978 (Fri Nov 07 05:56:26 PST 2003)
 */
public class MyPlot {
   MyJFrame frame;
   String title, xaxis, yaxis, marks = "none";
   ObjVector legends;
   boolean xlog = false, ylog = false, connected = true, bars = false;
   Plot plot = null;
   int dataSet = 0;

   public static MyPlot Plot( String title, String xaxis, String yaxis, double data[] ) {
      MyPlot plot = new MyPlot( title, xaxis, yaxis );
      plot.SetConnected( true ).SetMarks( "none" ).SetLog( false, false );
      plot.PlotXY( data ).FillPlot().QuitOnClose();
      return plot;
   }

   public static MyPlot Plot( String title, String xaxis, String yaxis, double x[], 
			      double y[] ) {
      MyPlot plot = new MyPlot( title, xaxis, yaxis );
      plot.SetConnected( true ).SetMarks( "none" ).SetLog( false, false );
      plot.PlotXY( x, y ).FillPlot().QuitOnClose();
      return plot;
   }

   public static MyPlot PlotHistogram( String title, String xaxis, String yaxis, double data[] ) {
      MyPlot plot = new MyPlot( title, xaxis, yaxis );
      plot.SetConnected( true ).SetMarks( "none" ).SetLog( false, false );
      plot.PlotHistogram( data ).FillPlot().QuitOnClose();
      return plot;
   }

   public static MyPlot PlotHistogram( String title, String xaxis, String yaxis, 
				       double data[], int nbins ) {
      MyPlot plot = new MyPlot( title, xaxis, yaxis );
      plot.SetConnected( true ).SetMarks( "none" ).SetLog( false, false );
      plot.PlotHistogram( data, nbins, data.length ).FillPlot().QuitOnClose();
      return plot;
   }

   public MyPlot( String title, String xaxis, String yaxis ) {
      this.title = title; this.xaxis = xaxis; this.yaxis = yaxis; PutInFrame(); }

   public MyJFrame GetFrame() { return frame; }

   protected void Init() {
      plot = new Plot(); dataSet = 0;
      PutInFrame();
   }

   protected void CreateFrame() {
      frame = new MyJFrame( title ); 
      frame.setDefaultCloseOperation( javax.swing.JFrame.DISPOSE_ON_CLOSE );
      frame.createUI( 300, 300 );
   }

   protected void PutInFrame() {
      if ( frame == null ) CreateFrame();
      if ( plot == null ) Init();
      frame.setComponent( plot );
   }

   public MyPlot SetTitle( String t ) { title = t; 
   if ( plot != null ) plot.setTitle( t ); return this; }
   public MyPlot SetXAxisTitle( String xa ) { xaxis = xa; return this; }
   public MyPlot SetYAxisTitle( String ya ) { yaxis = ya; return this; }
   public MyPlot SetLog( boolean x, boolean y ) { xlog = x; ylog = y; return this; }
   public MyPlot SetConnected( boolean x ) { connected = x; return this; }
   public MyPlot SetBars( boolean x ) { bars = x; return this; }
   // marks can be "none", "points", "dots", "various", "pixels"
   public MyPlot SetMarks( String x ) { marks = x; return this; }
   public MyPlot FillPlot() { getPlot().fillPlot(); return this; }
   public MyPlot QuitOnClose() { if ( frame != null ) frame.QuitOnClose(); return this; }
   public MyPlot AddListener( java.awt.event.WindowListener l ) { 
      if ( frame != null ) frame.addWindowListener( l ); return this; }

   public MyPlot AddLegend( String legend ) { 
      if ( legends == null ) legends = new ObjVector(); 
      legends.addElement( legend );
      getPlot().addLegend( legends.size()-1, legend ); 
      return this; }

   public MyPlot ClearLegends() { this.legends = null; getPlot().clearLegends(); 
   return this; }
   public MyPlot Clear() { 
      dataSet = 0; if ( plot != null ) plot.clear( false ); return this; }
   public MyPlot ChangeColor() { dataSet ++; return this; }

   public Plot getPlot() { if ( plot == null ) Init(); return plot; }

   public MyPlot AddPoint( double x, double y, boolean connected ) {
      getPlot().addPoint( dataSet, x, y, connected ); return this; }

   public MyPlot PlotXY( double data[] ) {
      return PlotXY( DoubleUtils.Sequence( data.length ), data ); }

   public MyPlot PlotXY( double xdata[], double ydata[] ) {
      boolean changed = false;
      if ( plot == null ) Init();
      else changed = true;

      plot.setTitle( title );
      //plot.setXRange( DoubleUtils.Min( xdata ), DoubleUtils.Max( xdata ) );
      //plot.setYRange( DoubleUtils.Min( ydata ), DoubleUtils.Max( ydata ) );
      plot.setXLabel( xaxis );
      plot.setYLabel( yaxis );
      plot.setXLog( xlog );
      plot.setYLog( ylog );
      plot.setConnected( connected, dataSet );
      plot.setBars( bars );
      plot.setMarksStyle( marks, dataSet );
      plot.setButtons( true );
      for ( int i = 0; i < xdata.length; i ++ ) {
	 plot.addPoint( dataSet, xdata[ i ], ydata[ i ], connected && ! changed );
	 changed = false;
      }
      //if ( legend != null ) { plot.addLegend( dataSet, legend ); legend = null; }
      frame.setTitle( title );
      PutInFrame();
      if ( ! frame.isVisible() ) frame.show();
      dataSet ++;
      return this;
   }

   public MyPlot PlotHistogram( double data[] ) {
      return PlotHistogram( data, 20 ); }

   public MyPlot PlotHistogram( double data[], int nbins ) {
      return PlotHistogram( data, nbins, data.length ); }

   public MyPlot PlotHistogram( double data[], int nbins, int maxInd ) {
      double[][] hist = DoubleUtils.Histogram( data, nbins );
      return PlotXY( hist[ 0 ], hist[ 1 ] );
   }

   public static void main( String args[] ) {
      double y[] = new double[ 100 ];
      for ( int i = 0; i < 100; i ++ ) y[ i ] = Math.sqrt( (double) i );
      MyPlot plot = Plot( "Test plot", "X axis", "Y axis", y );
      final MyJFrame frame = (MyJFrame) plot.GetFrame();
      frame.printToPostScript( "test.ps" );
      //javax.swing.SwingUtilities.invokeLater( new Runnable() { public void run() { frame.close(); } } );
      frame.close();
   }
}
