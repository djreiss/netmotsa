package djr.util.gui;

import java.io.*;
import java.awt.*;
import java.awt.event.*;

/**
 * Class <code>PrintStreamPanel</code>
 *
 * @author <a href="mailto:dreiss@systemsbiology.org">David Reiss</a>
 * @version 1.9978 (Fri Nov 07 05:56:26 PST 2003)
 */
public class PrintStreamPanel extends PrintStream implements ActionListener {
   protected MyFrame frame = null;
   protected TextArea area = null;
   protected Panel panel = null;
   protected ByteArrayOutputStream bos = null;

   public PrintStreamPanel() {
      super( System.out, true );
      area = new TextArea( "", 80, 24, TextArea.SCROLLBARS_VERTICAL_ONLY );

      panel = new Panel();
      panel.setLayout( new BorderLayout() );
      panel.add( BorderLayout.CENTER, area );

      bos = new ByteArrayOutputStream();
   }

   public void close() {
      if ( panel != null ) panel.removeAll(); 
      area = null;
      panel = null;
      if ( bos != null ) try { bos.close(); } catch( IOException e ) { e.printStackTrace(); } 
      bos = null;
      super.close(); 
   }

   public TextArea getArea() { return area; }
   public Panel getPanel() { return panel; }

   public void addButton( String name, ActionListener listener ) {
      if ( frame != null ) frame.addButton( name, listener );
   }

   public void write( int b ) { area.append( ( (char) b ) + "" ); }
   public void write( byte[] buf, int off, int len ) {
      bos.write( buf, off, len ); 
      if ( buf[ len-1 ] == (byte) '\n' || bos.size() > 128 ) flush(); }
   public void flush() {
      try { bos.flush(); 
      area.append( bos.toString() ); 
      bos.reset(); }
      catch( IOException e ) { };
   }	

   public void setForeground( Color c ) { area.setForeground( c ); }
   public void setBackground( Color c ) { area.setBackground( c ); }
   public void setFont( Font f ) { area.setFont( f ); }

   public void addStatus() { 
      if ( frame != null ) frame.addStatus(); }

   public void setStatus( String stat ) { 
      if ( frame != null ) frame.setStatus( stat ); }

   public void actionPerformed( ActionEvent evt ) {
      if ( "Save".equals( evt.getActionCommand() ) ) doSave();
   }

   protected void doSave() {
      PrintStream ps = frame.promptForFileSave();
      if ( ps != null ) {
	 String text = area.getText();
	 ps.println( text );
	 ps.close();
      }
   }

   public Frame putInFrame( String title ) {
      return putInFrame( title, 640, 480 );
   }

   public Frame putInFrame( String title, int width, int height ) {
      if ( frame != null ) return frame;
      frame = new MyFrame( title, area );
      frame.createUI( width, height + 50 );
      frame.add( BorderLayout.CENTER, this.getPanel() );
      frame.addButton( "Save", this );
      frame.show();
      return frame;
   }

   /*public static void main( String args[] ) {
      PrintStreamPanel psp = new PrintStreamPanel();
      psp.putInFrame( "PrintStreamPanel Test" );

      psp.println( "THIS IS A TEST OF THE PRINTSTREAM PANEL !!!" );
      for ( int i = 0; i < 100; i ++ ) psp.println( i );
   }   
   */
}
