package djr.util.gui;

import java.io.*;
import javax.swing.*;
import java.awt.*;
import java.awt.event.*;

/**
 * Class <code>JPrintStreamPanel</code>
 *
 * @author <a href="mailto:reiss@uw.edu">David Reiss</a>
 * @version 1.9978 (Fri Nov 07 05:56:26 PST 2003)
 */
public class JPrintStreamPanel extends PrintStream implements ActionListener {
   class MyEditorPane extends JEditorPane {
      public MyEditorPane() { super(); }
      public void scrollToRef( String ref ) { scrollToReference( ref ); }
   }

   protected JEditorPane area = null;
   protected JPanel panel = null;
   protected JScrollPane scrollPane = null;
   protected StringWriter sw = null;
   protected boolean isHTML = false;
   protected MyJFrame frame = null;

   public JPrintStreamPanel() {
      super( System.out, true );
      area = new MyEditorPane();
      scrollPane = new JScrollPane( area );

      panel = new JPanel();
      panel.setLayout( new BorderLayout() );
      panel.add( BorderLayout.CENTER, scrollPane );

      sw = new StringWriter();
   }

   public JEditorPane getArea() { return area; }
   public JPanel getPanel() { return panel; }

   public void addButton( String name, ActionListener listener ) {
      if ( frame != null ) frame.addButton( name, listener );
   }

   public void setContentHTML() { isHTML = true; area.setContentType( "text/html" ); }

   protected void appendText( String text ) {
      try {
	 sw.write( text ); sw.flush();
	 String s = sw.toString();
	 boolean hasRef = false, hasEndHtml = false;
	 if ( ! isHTML || ( hasRef = ( text.indexOf( "<a name=" ) >= 0 ) ) || 
	      text.indexOf( "<hr>" ) >= 0 ||
	      ( hasEndHtml = ( text.indexOf( "</html>" ) >= 0 ) ) ) { 
	    scrollPane.hide();
	    area.setText( s );
	    JScrollBar sb = scrollPane.getVerticalScrollBar();
	    sb.setValue( sb.getMaximum() + sb.getVisibleAmount() );
	    scrollPane.show();
	 }
      } catch( Exception e ) { e.printStackTrace(); }
   }

   public void write( int b ) { appendText( ( (char) b ) + "" ); }
   public void write( byte[] buf, int off, int len ) {
      appendText( new String( buf, off, len ) ); }

   public void setForeground( Color c ) { area.setForeground( c ); }
   public void setBackground( Color c ) { area.setBackground( c ); }
   public void setFont( Font f ) { area.setFont( f ); }

   public void addStatus() { 
      if ( frame != null ) frame.addStatus(); }

   public void setStatus( String stat ) { 
      if ( frame != null ) frame.setStatus( stat ); }

   public void close() {
      if ( panel != null ) panel.removeAll(); 
      area = null;
      panel = null;
      scrollPane = null;
      if ( sw != null ) try { sw.close(); } catch( IOException e ) { e.printStackTrace(); } 
      sw = null;
      super.close(); 
   }

   public void actionPerformed( ActionEvent evt ) {
      if ( "Save".equals( evt.getActionCommand() ) ) doSave();
   }

   protected void doSave() {
      PrintStream ps = frame.promptForFileSaveStream();
      if ( ps != null ) {
	 try { 
	    String text = sw.toString();
	    ps.println( text );
	    ps.flush(); ps.close();
	 } catch( Exception e ) { };
      }
   }

   public JFrame putInFrame( String title ) {
      return putInFrame( title, 640, 480 );
   }

   public JFrame putInFrame( String title, int width, int height ) {
      if ( frame != null ) return frame; 
      frame = new MyJFrame( title, this.getPanel() );
      frame.resize( width, height );
      frame.addButton( "Save", this );
      frame.center();
      frame.show();
      return frame;
   }

   /*public static void main( String args[] ) {
      JPrintStreamPanel psp = new JPrintStreamPanel();

      psp.putInFrame( "JPrintStreamPanel Test" );

      psp.println( "THIS IS A TEST OF THE PRINTSTREAM WINDOW !!!" );
      for ( int i = 0; i < 100; i ++ ) psp.println( i );
   }   
   */
}
