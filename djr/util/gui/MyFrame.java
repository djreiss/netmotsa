package djr.util.gui;
import java.awt.*;
import java.awt.event.*;
import java.awt.print.*;
import java.io.*;

/**
 * Class <code>MyFrame</code>.
 *
 * @author <a href="mailto:reiss@uw.edu">David Reiss</a>
 * @version 1.9978 (Fri Nov 07 05:56:26 PST 2003)
 */
public class MyFrame extends Frame implements WindowListener, ActionListener {
   int PANEL_HEIGHT = 30;

   protected String title;
   protected Panel buttonPanel;
   protected Label status;
   protected Component comp;
   
   public MyFrame( String title, Component comp ) {
      super( title );
      this.title = title;
      this.comp = comp;
   }

   public MyFrame( String title ) {
      this( title, null );
   }

   public void setComponent( Component c ) {
      if ( this.comp != null ) { remove( comp ); repaint(); }
      this.comp = c;
      add( BorderLayout.CENTER, comp );
   }

   public void actionPerformed( ActionEvent evt ) {
      if ( "Close".equals( evt.getActionCommand() ) ) close();
      else if ( "Print".equals( evt.getActionCommand() ) ) print();
   }

   public void print() {
      synchronized( this ) {
	 buttonPanel.hide();
	 PrinterJob pj = PrinterJob.getPrinterJob();
	 PageFormat pf = /*pj.pageDialog(*/ pj.defaultPage(); // );
	 ComponentPrintable cp = new ComponentPrintable( comp == null ? this : comp );
	 pj.setPrintable( cp, pf );
	 pj.setJobName( title );
	 if ( pj.printDialog() ) {
	    try { pj.print(); }
	    catch( PrinterException e ) { e.printStackTrace(); }
	 }
	 buttonPanel.show();
      }
   }

   /*public void printToImage( String file ) {
      synchronized( this ) {
	 String type = file.substring( file.lastIndexOf( '.' ) + 1 ).toLowerCase();
	 if ( "jpg".equals( type ) ) type = "jpeg";
	 try {
	    File temp = new File( ( new File( file ) ).getParent() );
	    temp.mkdirs();
	    if ( ! "gif".equals( type ) ) ;
	    //   com.sun.jimi.core.Jimi.putImage( "image/" + type, 
		//				createOffscreenImage(), file );
	    else {
	       FileOutputStream fos = new FileOutputStream( file );
	       ( new Acme.JPM.Encoders.GifEncoder( createOffscreenImage(), fos ) ).encode();
	    }
	 } catch( Exception e ) {
	    e.printStackTrace();
	 }
      }
   }
   */

   public Image createOffscreenImage() {
      synchronized( this ) {
	 Image offscreenImg = null;
	 try {
	    Component cp = comp == null ? this : (Component) comp;
	    offscreenImg = cp.createImage( cp.size().width, cp.size().height );
	    cp.paint( offscreenImg.getGraphics() );
	 } catch( Exception e ) {
	    offscreenImg = null;
	 }
	 return offscreenImg;
      }
   }

   public PrintStream promptForFileSave() {
      FileDialog fd = new FileDialog( this, "Choose file to save output to:", FileDialog.SAVE );
      fd.show();
      String file = fd.getFile();
      PrintStream ps = null;
      try { ps = new PrintStream( new FileOutputStream( file ) ); }
      catch( Exception e ) { ps = null; }
      return ps;
   }

   public void close() {
      hide();
      removeAll();
      buttonPanel = null;
      dispose();
   }

   public void addButton( String name, ActionListener listener ) {
      Button but = new Button( name );
      but.addActionListener( listener );
      buttonPanel.add( but );
      buttonPanel.repaint();
   }

   public Panel getButtonPanel() { return buttonPanel; }
   
   public void createUI( int width, int height ) {
      setSize( width, height + PANEL_HEIGHT );
      center();
      //if ( this instanceof JComponent ) ( (JComponent) this ).setDoubleBuffered( true ); 
      addWindowListener( this );
      setLayout( new BorderLayout() );
      buttonPanel = new Panel();
      addButton( "Print", this );
      addButton( "Close", this );
      this.add( BorderLayout.SOUTH, buttonPanel );
   }

   public void addStatus() { 
      status = new Label( "                   " );
      buttonPanel.add( status ); buttonPanel.doLayout(); }

   public void setStatus( String stat ) { 
      if ( status == null ) addStatus();
      int width = status.getGraphics().getFontMetrics().stringWidth( stat );
      if ( status.size().width < width ) status.setSize( width, status.size().height );
      status.setText( stat ); status.paint( status.getGraphics() ); }

   public void windowOpened( WindowEvent e ) { }
   public void windowClosing( WindowEvent e ) { this.close(); }
   public void windowClosed( WindowEvent e ) { this.close(); }
   public void windowIconified( WindowEvent e ) { }
   public void windowDeiconified( WindowEvent e ) { }
   public void windowActivated( WindowEvent e ) { }
   public void windowDeactivated( WindowEvent e ) { }
  
   public void center() {
      Dimension screenSize = Toolkit.getDefaultToolkit().getScreenSize();
      Dimension frameSize = getSize();
      int x = (screenSize.width - frameSize.width) / 2;
      int y = (screenSize.height - frameSize.height) / 2;
      setLocation(x, y);
   }
}
