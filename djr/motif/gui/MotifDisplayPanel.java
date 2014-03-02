package djr.motif.gui;
import java.util.*;
import java.awt.*;
import java.awt.event.*;
import java.io.*;
import javax.swing.*;

import djr.motif.sampler.*;
import djr.motif.model.*;
import djr.util.array.*;
import djr.util.gui.*;
import djr.util.bio.Sequence;

/**
 * Class <code>MotifDisplayPanel</code>.
 *
 * @author <a href="mailto:dreiss@systemsbiology.org">David Reiss</a>
 * @version 1.9978 (Fri Nov 07 05:56:26 PST 2003)
 */
public class MotifDisplayPanel extends JPanel implements MouseListener, ActionListener {
   static double H_SCALE = 1; // Pixels per nucleotide/acid in sequence
   static double V_SCALE = 10; // Pixels per sequence
   static final int MAX_WIDTH = 800; // Max pixel width of window
   static final int MAX_HEIGHT = 600; // Max pixel height of window
   static final int MIN_WIDTH = 300; // Min pixel width of window
   static final int MARGIN = 40; // In pixels
   static final Color colors[] = { Color.red, Color.green, Color.blue, Color.yellow,
				   Color.orange, Color.cyan, Color.magenta, Color.pink,
				   Color.black, Color.darkGray, Color.lightGray };

   Sampler sampler;
   ObjVector states, listeners;
   int width, height, currMot = -1;
   ObjVector rects, mainRects, seqRects, alignments, pvalues;
   boolean selected[];
   IntVector rectInfos, mainRectInfos;
   ObjVector seqsByMotif;
   String title;
   JFrame frame;
   
   public MotifDisplayPanel( String title, ObjVector states, Sampler samp ) {
      this.states = states;
      this.sampler = samp;
      this.title = title;
      int maxLen = sampler.GetSites().GetMaxSequenceLength();
      this.width = (int) ( H_SCALE * maxLen + 2 * MARGIN );
      boolean worc = sampler.GetArgProcessor().getBooleanArg( "worc" );
      maxLen = worc ? maxLen / 2 : maxLen;
      this.height = (int) ( V_SCALE * ( sampler.GetSites().GetSequenceCount() + 1 ) + 2 * MARGIN );
      this.addMouseListener( this );

      alignments = new ObjVector();
      pvalues = new ObjVector();
      for ( int mot = 0, sz = states.size(); mot < sz; mot ++ ) {
	 Map mstate = (Map) states.elementAt( mot );
	 sampler.SetState( mstate );
	 AlignmentMotifModel mmodel = sampler.GetMotifModel();
	 int aa[][][] = mmodel.GetSites().GetAlignments();
	 double pvals[][] = mmodel.ComputePolyScores( aa );
	 alignments.addElement( aa );
	 pvalues.addElement( pvals );
      }
   }

   public MyJFrame putInFrame() {
      frame = new MyJFrame( title, this );
      int ww = (int) ( H_SCALE * sampler.GetSites().GetMaxSequenceLength() + 2 * MARGIN );
      if ( ww > MAX_WIDTH ) ww = MAX_WIDTH;
      if ( ww < MIN_WIDTH ) ww = MIN_WIDTH;
      int hh = (int) ( V_SCALE * ( sampler.GetSites().GetSequenceCount() + 1 ) + 2 * MARGIN );
      if ( hh > MAX_HEIGHT ) hh = MAX_HEIGHT;
      JScrollPane scrollPane = new JScrollPane( this );
      frame.getContentPane().add( BorderLayout.CENTER, scrollPane );
      ( ( MyJFrame) frame ).createUI( ww + 55, hh + 55 );
      ( ( MyJFrame) frame ).addButton( "Save", this );
      ( ( MyJFrame) frame ).addButton( "Info", this );
      frame.setVisible( true );
      return (MyJFrame) frame;
   }

   public MyJFrame getFrame() { return (MyJFrame) frame; }

   public void actionPerformed( ActionEvent evt ) {
      if ( "Save".equals( evt.getActionCommand() ) ) {
	 String toSave = ( (MyJFrame) frame ).promptForFileSave();
	 if ( toSave != null ) sampler.SaveStateToFile( toSave );
      } else if ( "Info".equals( evt.getActionCommand() ) ) {
	 Thread thr = new Thread() { public void run() {
	    sampler.SetUseHTMLOutput();
	    if ( currMot < 0 ) {
	       ( new Thread() { public void run() { 
		  sampler.PrintBestMotifs( states ); } } ).start();
	       LogoPanel.showLogos( "Motif Logo Window", states, null, null, 220 );
	    } else {
	       ( new Thread() { public void run() { 
		  sampler.PrintBestMotif( (Map) states.elementAt( currMot ), "#"+(currMot+1) );
	       } } ).start();
	       LogoPanel.showLogos( "Motif Logo Window", sampler.GetMotifModel(), 
				    new String[] { "Logo for Motif #"+(currMot+1) }, 
				    null, 220 );
	    }

	 } };
	 thr.setPriority( Thread.MIN_PRIORITY );
	 thr.setDaemon( true );
	 thr.start();
      } else {
	 transmitEvent( evt );
      }
   }

   public void addActionListener( ActionListener list ) {
      if ( listeners == null ) listeners = new ObjVector();
      if ( ! listeners.contains( list ) ) listeners.addElement( list );
   }

   public void removeActionListener( ActionListener list ) {
      if ( listeners == null ) return;
      if ( listeners.contains( list ) ) listeners.removeElement( list );
   }

   public void removeNotify() {
      super.removeNotify();
      listeners = states = rects = mainRects = seqRects = alignments = pvalues = null;
      selected = null;
      rectInfos = mainRectInfos = null;
      sampler = null;
      frame = null;
   }

   public Dimension getPreferredSize() { return new Dimension( width, height ); }
   public Dimension getMinimumSize() { return new Dimension( width, height ); }

   public void paint( Graphics g ) {
      Graphics2D g2 = (Graphics2D) g;
      g2.setRenderingHint(RenderingHints.KEY_ANTIALIASING,
			  RenderingHints.VALUE_ANTIALIAS_ON);
      Rectangle r = g2.getClipRect();
      g2.setColor( Color.white );
      g2.fillRect( r.x, r.y, r.width, r.height );

      rects = new ObjVector(); rectInfos = new djr.util.array.IntVector();
      mainRects = new ObjVector(); mainRectInfos = new djr.util.array.IntVector();
      boolean worc = sampler.GetArgProcessor().getBooleanArg( "worc" );
      int maxLen = worc ? sampler.GetSites().GetMaxSequenceLength() / 2 : sampler.GetSites().GetMaxSequenceLength();
      int NS = sampler.GetSites().GetSequenceCount(), len[] = (int[]) sampler.GetSites().GetSequenceLengths();
      Sequence S[] = (Sequence[]) sampler.GetSites().GetSequences();
      double xscale = (double) size().width / (double) width;
      double yscale = (double) size().height / (double) height;
      double vscale = V_SCALE * yscale;
      int margin = (int) ( MARGIN * Math.min( xscale, yscale ) );
      g.translate( margin, margin );
      double hscale = ( (double) ( size().width - 3 * margin ) ) / ( (double) maxLen );
      FontMetrics fm = g2.getFontMetrics( g2.getFont() );
      int fheight = fm.getAscent();
      boolean useAlpha = isDoubleBuffered(); // Not double buffered = printing to paper
      seqRects = new ObjVector();
      for ( int i = 0; i < NS; i ++ ) {
	 g2.setColor( Color.blue );
	 int llen = worc ? len[ i ] / 2 : len[ i ];
	 llen = (int) ( llen * hscale );
	 int ycur = (int) ( vscale / 2 + vscale * i );
	 g2.drawLine( 0, ycur, llen, ycur );
	 seqRects.addElement( new Rectangle( margin-3, margin+ycur-3, llen+7, 7 ) );
	 if ( useAlpha && selected != null && selected[ i ] ) {
	    g2.drawLine( 0, ycur - 1, llen, ycur - 1 );
	    g2.drawLine( 0, ycur + 1, llen, ycur + 1 );
	 }
	 g2.setColor( Color.black );
	 String ss = (i+1) + "";
	 int fwidth = fm.stringWidth( ss );
	 g2.drawString( ss, -fwidth, ycur + fheight/2 );
	 String head = S[ i ].GetHeader();
	 if ( head.length() > 9 ) head = head.substring( 0, 9 );
	 g2.drawString( head, llen + 3, ycur + fheight/2 );
      }
      if ( selected == null ) selected = new boolean[ seqRects.size() ];
      int stepsize = 50; 
      if ( hscale > 10 ) stepsize = 5;
      else if ( hscale > 5 ) stepsize = 10;
      else if ( hscale > 2 ) stepsize = 25;
      else if ( hscale > 1 ) stepsize = 50;
      else if ( hscale > 0.5 ) stepsize = 100;
      else stepsize = 200;
      for ( int i = 1; i <= maxLen; i += stepsize ) {
	 if ( i != 1 && ( maxLen - i ) * hscale > 10 ) {
	    g.setColor( Color.lightGray );
	    g.drawLine( (int) ( hscale * (i-1) ), 0, (int) ( hscale * (i-1) ), 
			(int) ( vscale * NS + 5 ) );
	 }
	 String ss = i + "";
	 int fwidth = fm.stringWidth( ss );
	 g.setColor( Color.black );
	 g2.drawString( ss, (int) ( hscale * (i-1) ) - fwidth/2, -fheight/2 );
      }
      if ( ( maxLen + 1 ) % stepsize > stepsize / 2 ) {
	 String ss = (maxLen+1) + "";
	 int fwidth = fm.stringWidth( ss );
	 g2.drawString( ss, (int) ( hscale * maxLen ) - fwidth/2, -fheight/2 );
      }
      for ( int mot = 0, sz = states.size(); mot < sz; mot ++ ) {
	 g2.setColor( Color.black );
	 g2.drawString( (mot+1)+"=", -10 + mot * (int) ( vscale * 2 + 
					    fm.stringWidth( (mot+1)+"=" ) + 3 ), 
		       (int) ( vscale * ( NS + 1 ) + 10 * vscale/V_SCALE ) );
	 g2.setColor( colors[ mot ] );
	 r = new Rectangle( -10 + mot * (int) ( vscale * 2 ) + 
			    (mot+1) * ( fm.stringWidth( (mot+1)+"=" ) + 3 ),  
			    (int) ( vscale * NS + 15 * vscale/V_SCALE ), 
			    (int) ( vscale * 2 - 2 ), (int) ( vscale / 2 ) );
	 g2.fillRect( r.x, r.y, r.width, r.height );
	 r.translate( margin, margin );
	 r.grow( 2, 2 );
	 mainRects.addElement( r );
	 mainRectInfos.addElement( mot );

	 int aa[][][] = (int[][][]) alignments.elementAt( mot );
	 double pvals[][] = (double[][]) pvalues.elementAt( mot );
	 for ( int m = 0; m < aa.length; m ++ ) {
	    if ( aa[ m ] == null ) continue;
	    for ( int i = 0; i < aa[ m ].length; i ++ ) {
	       int seq = aa[ m ][ i ][ 0 ];
	       if ( seq < 0 ) continue;
	       int locind = 1, loc = aa[ m ][ i ][ locind ];
	       while ( loc < 0 && locind ++ < aa[ m ][ i ].length - 1 ) 
		  loc = aa[ m ][ i ][ locind ];
	       if ( loc < 0 ) continue;
	       if ( worc && loc > len[ seq ] / 2 ) loc = len[ seq ] - loc;
	       int lastind = -1;
	       for ( int j = 1; j < aa[ m ][ i ].length; j ++ ) {
		  int jj = aa[ m ][ i ][ j ];
		  if ( jj > 0 ) lastind = jj;
	       }
	       if ( worc && lastind > len[ seq ] / 2 ) lastind = len[ seq ] - lastind;
	       if ( lastind < loc ) { int temp = loc; loc = lastind; lastind = temp; }
	       r = new Rectangle( (int) ( loc * hscale ),
				  (int) ( vscale/2 + seq * vscale - vscale/4 ) + 1,
				  (int) ( ( lastind - loc + 1 ) * hscale ), (int)( vscale / 2 ) );
	       if ( currMot == -1 || mot == currMot ) {
		  Color c = colors[ mot ];
		  g2.setColor( c );
		  if ( useAlpha ) {
		     double pp = pvals != null && pvals[ m ] != null ? pvals[ m ][ i ] : 
			sampler.GetMotifModel().ComputeMatchScore( aa, i, m );
		     pp = -DoubleUtils.Log10( pp == 0.0 ? 1.0 / Double.MAX_VALUE : pp );
		     if ( pp < 2.0 ) {
			int alpha = 255;
			for ( double q = 0; q < 2.0 - pp; q += 0.1 ) {	
			   alpha -= 12; if ( alpha < 0 ) { alpha = 0; break; }
			}
			g2.clearRect( r.x, r.y, r.width, r.height );
			g2.setColor( new Color( c.getRed(), c.getGreen(), c.getBlue(), alpha ) );
		     }
		  }
		  g2.fillRect( r.x, r.y, r.width, r.height );
	       }
	       r.translate( margin, margin );
	       r.grow( 2, 2 );
	       rects.addElement( r );
	       rectInfos.addElement( mot );
	    }
	 }
      }
   }   

   public void mouseClicked( MouseEvent e ) {
      if ( rects == null ) return;
      int x = e.getX(), y = e.getY();
      boolean gotClick = false;
      for ( int i = 0, size = mainRects.size(); i < size; i ++ ) {
	 Rectangle r = (Rectangle) mainRects.elementAt( i );
	 if ( r.contains( x, y ) ) {
	    currMot = mainRectInfos.elementAt( i );
	    getParent().repaint();
	    if ( listeners != null ) {
	       ActionEvent evt = new ActionEvent( new Integer( currMot ), 
				 ActionEvent.ACTION_PERFORMED, "select motif", currMot );
	       transmitEvent( evt );
	    }
	    gotClick = true;
	    break;
	 }
      }

      if ( ! gotClick ) {
	 for ( int i = 0, size = rects.size(); i < size; i ++ ) {
	    Rectangle r = (Rectangle) rects.elementAt( i );
	    if ( r.contains( x, y ) ) {
	       final int mot = rectInfos.elementAt( i );
	       if ( currMot != -1 && currMot != mot ) continue;
	       if ( listeners != null ) {
		  ActionEvent evt = new ActionEvent( new Integer( currMot ), 
						     ActionEvent.ACTION_PERFORMED, "display motif", 
						     currMot != -1 ? currMot : 0 );
		  transmitEvent( evt );
	       }
	       Thread thr = new Thread() { public void run() {
		  sampler.SetUseHTMLOutput();
		  ( new Thread() { public void run() { 
		     sampler.PrintBestMotif( (Map) states.elementAt( mot ), "#"+(mot+1) );
		  } } ).start();
		  //sampler.SetState( (Map) states.elementAt( mot ) );
		  LogoPanel.showLogos( "Motif Logo Window", sampler.GetMotifModel(), 
				       new String[] { "Logo for Motif #"+(mot+1) }, 
				       null, 220 );
	       } };
	       thr.setPriority( Thread.MIN_PRIORITY );
	       thr.setDaemon( true );
	       thr.start();
	       gotClick = true;
	       break;
	    }	
	 }
      }

      if ( ! gotClick ) {
	 for ( int i = 0, size = seqRects.size(); i < size; i ++ ) {
	    Rectangle r = (Rectangle) seqRects.elementAt( i );
	    if ( r.contains( x, y ) ) {
	       setSelected( i, ! selected[ i ], true );
	       gotClick = true;
	       break;
	    }
	 }
      }

      if ( ! gotClick ) currMot = -1;
      repaint();
   }

   public void setSelected( int seq, boolean sel, boolean repaint ) {
      if ( selected == null && seqRects != null ) selected = new boolean[ seqRects.size() ];
      if ( selected == null ) return;
      selected[ seq ] = sel;
      if ( listeners != null ) {
	 ActionEvent evt = new ActionEvent( sampler.GetSites().GetSequences()[ seq ], 
					    ActionEvent.ACTION_PERFORMED, 
					    sel ? "select sequence" : "deselect sequence", seq );
	 transmitEvent( evt );
      }
      if ( repaint ) repaint();
   }

   public void selectAll( boolean sel ) {
      if ( selected == null && seqRects == null ) return;
      for ( int i = 0; i < seqRects.size(); i ++ ) setSelected( i, sel, false );
      repaint();
   }

   public void selectAllWithCurrentMot() {
      for ( int mot = 0; mot < states.size(); mot ++ ) {
	 if ( currMot >= 0 && mot != currMot ) continue;
	 int aa[][][] = (int[][][]) alignments.elementAt( mot );
	 double pvals[][] = (double[][]) pvalues.elementAt( mot );
	 for ( int m = 0; m < aa.length; m ++ ) {
	    if ( aa[ m ] == null ) continue;
	    for ( int i = 0; i < aa[ m ].length; i ++ ) {
	       int seq = aa[ m ][ i ][ 0 ];
	       if ( selected != null && selected[ seq ] ) continue;
	       double pp = pvals != null && pvals[ m ] != null ? pvals[ m ][ i ] : 
		  sampler.GetMotifModel().ComputeMatchScore( aa, i, m );
	       pp = -DoubleUtils.Log10( pp == 0.0 ? 1.0 / Double.MAX_VALUE : pp );
	       if ( pp > 0.0 ) setSelected( seq, true, false );
	    }
	 }
      }
      repaint();
   }

   public ObjVector getStates() { return states; }

   protected void transmitEvent( ActionEvent evt ) {
      if ( listeners == null ) return;
      for ( int j = 0, s = listeners.size(); j < s; j ++ ) {
	 ActionListener al = (ActionListener) listeners.elementAt( j );
	 if ( al != this ) al.actionPerformed( evt );
      }
   }

   public void mouseEntered(MouseEvent e) { }
   public void mouseExited(MouseEvent e) { }
   public void mousePressed(MouseEvent e) { }
   public void mouseReleased(MouseEvent e) { }

   /*public static void main( String args[] ) {
      gnu.getopt.Getopt g = new gnu.getopt.Getopt( "MotifDisplayPanel", args, "f:" );
      String fname = "zzz"; int c;
      while ( ( c = g.getopt() ) != -1 ) {
	 switch( c ) {
	 case 'f': fname = g.getOptarg(); break;
	 case '?': break; // getopt() already printed an error
	 default: System.out.print( "getopt() returned " + c + "\n" );
	 }
      }
      try {
	 Sampler sampler = (Sampler) djr.util.MyUtils.ReadObject( fname ); 
	 ObjVector states = new ObjVector();
	 states.addElement( sampler.GetState() );
	 MotifDisplayPanel mp = new MotifDisplayPanel( fname, states, sampler );
	 mp.putInFrame();
      } catch( Exception e ) { e.printStackTrace(); }
   }
   */
}
