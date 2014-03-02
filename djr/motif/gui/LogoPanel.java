package djr.motif.gui;
import java.awt.*;
import java.awt.font.*;
import java.awt.geom.*;
import java.awt.event.*;
import java.util.*;
import javax.swing.*;

import djr.motif.sampler.*;
import djr.motif.model.AlignmentMotifModel;
import djr.util.array.*;
import djr.util.gui.*;
import djr.util.bio.Sequence;

/**
 * Class <code>LogoPanel</code>.
 *
 * @author <a href="mailto:dreiss@systemsbiology.org">David Reiss</a>
 * @version 1.9978 (Fri Nov 07 05:56:26 PST 2003)
 */
public class LogoPanel extends JPanel {
   final static Font font = new Font( "Monospaced", Font.PLAIN, 24 );
   final static double CHARWIDTH = 30.0; // Width of column in logo
   final static String XAXIS = "Motif Position", YAXIS = "Bits";
   static int MARGIN = 40;
   static Font fontAxis = new Font( "Serif", Font.PLAIN, 12 );

   protected double logo[][], max, min, errBar[];
   protected int temp[][], width, height;
   protected boolean yesTitle;
   protected String title;

   public LogoPanel( Sampler samp, int mot, String title ) {
      this( samp.GetMotifModel().ComputeLogoScores( mot ), title );
   }

   public LogoPanel( double logo[][], String title ) { 
      this.logo = logo;
      this.title = title;
      this.width = (int) CHARWIDTH * logo.length + 2 * MARGIN;
      this.height = 400;

      errBar = DoubleUtils.New( logo.length );
      double JJ = (double) logo[ 0 ].length;
      DoubleUtils.Set( errBar, -1/JJ * DoubleUtils.Log2( 1/JJ ) );
      max = -Double.MAX_VALUE; min = Double.MAX_VALUE;
      for ( int i = 0; i < logo.length; i ++ ) {
	 max = DoubleUtils.Max( max, DoubleUtils.Sum( logo[ i ] ) );
	 for ( int j = 0; j < logo[ i ].length; j ++ )
	    min = DoubleUtils.Min( min, logo[ i ][ j ] );
      }
      temp = new int[ logo.length ][];
      for ( int ii = 0; ii < logo.length; ii ++ ) 
	 temp[ ii ] = DoubleUtils.Indexx( logo[ ii ], null );
   }

   public static JFrame showLogos( String title, AlignmentMotifModel mmodel, 
				   String titles[], JFrame frame ) {
      return showLogos( title, mmodel, titles, frame, 300 );
   }

   public static JFrame showLogos( String title, AlignmentMotifModel mmodel, 
				   String titles[], JFrame frame, int size ) {
      int nm = mmodel.GetMotifCount();
      JPanel panel = new JPanel();
      int grid1 = (int) Math.sqrt( (double) nm ), grid2 = nm / grid1;
      panel.setLayout( new GridLayout( grid2, grid1 ) );
      for ( int i = 0; i < nm; i ++ ) {
	 double pssm[][] = mmodel.GetMotifCounts( i );
	 LogoPanel lp = new LogoPanel( AlignmentMotifModel.ComputeLogoScores( pssm ), 
				       titles != null ? titles[ i ] : "Motif #" + (i+1) );
	 lp.setShowTitle( true );
	 panel.add( lp );
      }
      if ( frame == null ) {
	 frame = new MyJFrame( title );
	 int ww = size * grid1, hh = size * grid2;
	 if ( hh > 1.2 * ww ) hh = (int) ( 1.2 * ww );
	 frame.resize( ww, hh + 50 );
      }
      ( (MyJFrame) frame ).setTitle( title ); 
      ( (MyJFrame) frame ).setComponent( panel ); 
      /*if ( ! frame.isVisible() )*/ frame.show();
      return frame;
   }

   public static JFrame showLogos( String title, ObjVector states, 
				   String titles[], JFrame frame ) {
      return showLogos( title, states, titles, frame, 300 );
   }

   public static JFrame showLogos( String title, ObjVector states, 
				   String titles[], JFrame frame, int size ) {
      int nm = states.size();
      JPanel panel = new JPanel();
      int grid1 = (int) Math.sqrt( (double) nm ), grid2 = nm / grid1;
      panel.setLayout( new GridLayout( grid2, grid1 ) );
      for ( int i = 0; i < nm; i ++ ) {
	 Map state = (Map) states.get( i );
	 AlignmentMotifModel mmodel = (AlignmentMotifModel) state.get( "mm" );
	 double pssm[][] = mmodel.GetMotifCounts( 0 );
	 LogoPanel lp = new LogoPanel( AlignmentMotifModel.ComputeLogoScores( pssm ), 
				       titles != null ? titles[ i ] : "Motif #" + (i+1) );
	 lp.setShowTitle( true );
	 panel.add( lp );
      }
      if ( frame == null ) {
	 frame = new MyJFrame( title );
	 int ww = size * grid1, hh = size * grid2;
	 if ( hh > 1.2 * ww ) hh = (int) ( 1.2 * ww );
	 frame.resize( ww, hh + 50 );
      }
      ( (MyJFrame) frame ).setTitle( title ); 
      ( (MyJFrame) frame ).setComponent( panel ); 
      /*if ( ! frame.isVisible() )*/ frame.show();
      return frame;
   }

   public MyJFrame putInFrame() {
      return putInFrame( null, this.title ); }

   public MyJFrame putInFrame( JFrame frame, String titl ) {
      this.title = titl;
      if ( frame == null ) {
	 frame = new MyJFrame( title, this );
	 ( (MyJFrame) frame ).createUI( width, height + 50 );
      } else ( (MyJFrame) frame ).setComponent( this );
      frame.show();
      return (MyJFrame) frame;
   }

   public Dimension getPreferredSize() { return new Dimension( width, height ); }
   public Dimension getMinimumSize() { return new Dimension( width, height ); }

   public LogoPanel setShowTitle( boolean show ) { yesTitle = show; return this; }
 
   public void paint( Graphics g ) {
      int ww = size().width, hh = size().height;
      if ( ww < 100 || hh < 100 ) {
	 MARGIN = 10;	
	 fontAxis = new Font( "Sans", Font.PLAIN, 6 );      
	 this.width = (int) CHARWIDTH * logo.length + 2 * MARGIN;
      } else {
	 MARGIN = 40;	
	 fontAxis = new Font( "Sans", Font.PLAIN, 12 );      
	 this.width = (int) CHARWIDTH * logo.length + 2 * MARGIN;
      }

      Graphics2D g2 = (Graphics2D) g;
      g2.setFont( fontAxis );
      g2.setRenderingHint(RenderingHints.KEY_ANTIALIASING,
			  RenderingHints.VALUE_ANTIALIAS_ON);
      Rectangle r = g2.getClipRect();
      g2.setColor( Color.white );
      g2.fillRect( r.x, r.y, r.width, r.height );

      FontRenderContext frc = g2.getFontRenderContext();
      g2.translate( MARGIN, MARGIN );
      short J = (short) logo[ 0 ].length;
      String s = J == 4 ? Sequence.NUCLEOTIDES : Sequence.ACIDS;
      GlyphVector gv = font.createGlyphVector( frc, s );
      double totHeight = (double) ( hh - 2*MARGIN );
      double bitScale = totHeight / max; // Pixels per bit

      double stepSize = 0.1;
      FontMetrics fm = g.getFontMetrics( fontAxis );
      int yoff = fm.getAscent() / 4;
      while( bitScale * stepSize < fm.getHeight() ) stepSize *= 1.5;
      for ( double i = 0, pix = totHeight; i <= max; i += stepSize, pix -= bitScale * stepSize ) {
	 String str = Math.round( i * 10 ) / 10.0 + "";
	 g.drawString( str, -fm.stringWidth( str )-3, (int) pix + yoff );
	 g.setColor( Color.lightGray );
	 g.drawLine( 0, (int) pix, ww-2*MARGIN, (int) pix );
	 g.setColor( Color.black );
	 g.drawLine( 0, (int) pix, 3, (int) pix );
      }
      double charwidth = CHARWIDTH * ww / width;
      for ( int i = 0, pix = (int) charwidth / 2 - 5; i < logo.length; i ++, pix += (int) charwidth ) {
	 g.drawString( (i+1) + "", pix, (int) totHeight + fm.getHeight()+3 ); }
      g.drawString( XAXIS, ( ww - fm.stringWidth( XAXIS ) - MARGIN ) / 2, 
		    (int) totHeight + fm.getHeight()*2 + 3 );
      if ( yesTitle ) g.drawString( title, ( ww - fm.stringWidth( title ) - MARGIN ) / 2, 
		    fm.getHeight() + 3 - MARGIN );
      AffineTransform at = new AffineTransform();
      GlyphVector gl = fontAxis.createGlyphVector( frc, YAXIS );
      Shape shp = gl.getOutline();
      at.setToRotation( -Math.PI / 2 );
      Shape tshp = at.createTransformedShape( shp );
      Rectangle b = tshp.getBounds();
      at.setToTranslation( -fm.stringWidth( " 0.0" ), ( hh - b.height ) / 2 );
      tshp = at.createTransformedShape( tshp );
      g2.fill( tshp );

      Shape glyphs[] = new Shape[ J ];
      Color colors[] = new Color[ J ]; 
      for ( int i = 0; i < logo.length; i ++ ) {
	 double top = totHeight;
	 for ( int j = 0; j < logo[ 0 ].length; j ++ ) {
	    int ntide = temp[ i ][ j ]; 
	    //if ( logo[ i ][ ntide ] <= min ) continue;
	    if ( glyphs[ ntide ] == null ) 
	       glyphs[ ntide ] = gv.getGlyphOutline( ntide ); // Get the original glyph
	    Shape glyph = glyphs[ ntide ];
	    double charHeight = glyph.getBounds().height;
	    double scale = ( logo[ i ][ ntide ] /*- min*/ ) / max * totHeight / charHeight;
	    scale = Math.abs( scale );

	    at.setToScale( 2.0 * ww / width, scale ); // Scale it by the proper amount
	    Shape transformedGlyph = at.createTransformedShape( glyph );
	    double yyoff = transformedGlyph.getBounds().height + transformedGlyph.getBounds().y;
	    top -= yyoff;

	    at.setToTranslation( -transformedGlyph.getBounds().x + i * charwidth, top ); 
	    transformedGlyph = at.createTransformedShape( transformedGlyph ); // Translate it
	    
	    //Rectangle rr = transformedGlyph.getBounds();
	    //g2.drawRect(rr.x,rr.y,rr.width,rr.height);

	    String color = Sequence.getHTMLColorBG( (short) ntide, J );
	    if ( colors[ ntide ] == null ) 
	       colors[ ntide ] = Color.decode( "0x" + color.substring( 1 ) ).darker();
	    g2.setColor( colors[ ntide ] );
	    g2.fill( transformedGlyph );
	    top -= transformedGlyph.getBounds().height;
	    top += yyoff;
	 }
	 /*g.setColor( Color.black );
	 g.drawLine( (int) (i*charwidth+charwidth/2), (int) (top-errBar[ i ]/2*bitScale),
		     (int) (i*charwidth+charwidth/2), (int) (top+errBar[ i ]/2*bitScale) );
	 g.drawLine( (int) (i*charwidth+charwidth/4), (int) (top-errBar[ i ]/2*bitScale),
		     (int) (i*charwidth+charwidth*3/4), (int) (top-errBar[ i ]/2*bitScale) );
	 g.drawLine( (int) (i*charwidth+charwidth/4), (int) (top+errBar[ i ]/2*bitScale),
		     (int) (i*charwidth+charwidth*3/4), (int) (top+errBar[ i ]/2*bitScale) );
	 */
      }
      g.setColor( Color.black );
      g.drawLine( 0, 0, 0, (int) totHeight );
      g.drawLine( 0, (int) totHeight, ww-2*MARGIN, (int) totHeight );      
   }
}
