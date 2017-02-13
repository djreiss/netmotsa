package djr.util.gui;
import java.awt.*;
import java.awt.event.*;
import java.io.*;
import javax.swing.*;
import java.awt.print.*;
import java.util.*;

import djr.util.MyUtils;

/**
 * Class <code>MyJFrame</code>.
 *
 * @author <a href="mailto:reiss@uw.edu">David Reiss</a>
 * @version 1.9978 (Fri Nov 07 05:56:26 PST 2003)
 */
public class MyJFrame extends JFrame implements WindowListener, ActionListener {
   final int PANEL_HEIGHT = 30;

   protected String title;
   protected JPanel buttonPanel;
   protected JLabel status = null;
   protected JComponent comp = null;
   protected boolean quitOnClose = false;
   protected Map buttons;
   
   public MyJFrame( String title, JComponent comp ) {
      super( title );
      this.title = title;
      getContentPane().setLayout( new BorderLayout() );
      if ( comp != null ) setComponent( comp );
   }

   public MyJFrame( String title ) {
      this( title, null );
   }

   public void setComponent( JComponent c ) {
      if ( buttonPanel == null ) createUI( getWidth(), getHeight() );
      if ( c != this.comp ) {
	 if ( this.comp != null ) { 
	    getContentPane().remove( comp ); getContentPane().repaint(); }
	 this.comp = c;
	 getContentPane().add( BorderLayout.CENTER, comp );
      }
      getContentPane().doLayout();
      getContentPane().repaint();
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
	 ComponentPrintable cp = new ComponentPrintable( comp == null ? 
							 this.getContentPane() : (Component) comp );
	 pj.setPrintable( cp, pf );
	 pj.setJobName( title );
	 if ( pj.printDialog() ) {
	    try { pj.print(); }
	    catch( PrinterException e ) { e.printStackTrace(); }
	 }
	 buttonPanel.show();
      }
   }

   public void printToPostScript( final String fname ) {
      synchronized( this ) {
	 try {
	    Thread.yield();
	    SwingUtilities.invokeLater( new Runnable() { public void run() {
	       buttonPanel.hide();
	       System.err.println( "Printing to postscript file " + fname + "..." ); 
	       ComponentPrintable cp = new ComponentPrintable( comp == null ? 
							       MyJFrame.this.getContentPane() : (Component) comp );
	       cp.printToPostScriptFile( fname );
	       System.err.println( "Done printing." ); 
	       buttonPanel.show();
	    } } );
	 } catch ( Exception e ) {
	    e.printStackTrace();
	 }
      }
   }

   //static String[] types = {"png","bmp","pict","jpeg","xbm","tga","psd","xpm","pcx"};

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
	    Component cp = comp == null ? this.getContentPane() : (Component) comp;
	    synchronized( cp ) {
	       offscreenImg = cp.createImage( cp.size().width, cp.size().height );
	       cp.paint( offscreenImg.getGraphics() );
	    }
	 } catch( Exception e ) {
	    offscreenImg = null;
	 }
	 return offscreenImg;
      }
   }

   public String promptForFileSave() {
      String currDir = System.getProperty( "user.dir" );
      JFileChooser fd = new JFileChooser( currDir );
      fd.setDialogTitle( "Choose file to save output to:" );
      fd.setDialogType( JFileChooser.SAVE_DIALOG );
      int returnVal = fd.showSaveDialog( this ); 
      if ( returnVal == JFileChooser.APPROVE_OPTION )
	 return fd.getSelectedFile().getName();
      return null;
   }

   public PrintStream promptForFileSaveStream() {
      String currDir = System.getProperty( "user.dir" );
      JFileChooser fd = new JFileChooser( currDir );
      fd.setDialogTitle( "Choose file to save output to:" );
      fd.setDialogType( JFileChooser.SAVE_DIALOG );
      int returnVal = fd.showSaveDialog( this ); 
      PrintStream ps = null;
      if ( returnVal == JFileChooser.APPROVE_OPTION ) {
	 String file = fd.getSelectedFile().getName();
	 try { ps = new PrintStream( new FileOutputStream( file ) ); }
	 catch( Exception e ) { ps = null; }
      }
      return ps;
   }

   public void close() {
      try {
	 Thread.yield();
	 SwingUtilities.invokeLater( new Runnable() { public void run() {
	    hide();
	    WindowListener ls[] = getWindowListeners();
	    //for ( int i = 0; i < ls.length; i ++ ) ls[ i ].windowClosing( null );
	    getContentPane().removeAll();
	    comp = null;
	    buttonPanel = null;
	    status = null;
	    //dispose();
	    if ( quitOnClose ) MyUtils.Exit( 0 );
	    //for ( int i = 0; i < ls.length; i ++ ) ls[ i ].windowClosed( null );
	 } } );
      } catch ( Exception e ) { e.printStackTrace(); }
   }

   public JButton addButton( String name, ActionListener listener ) {
      JButton but = new JButton( name );
      if ( listener == null && listener != comp && comp instanceof ActionListener ) 
	 but.addActionListener( (ActionListener) comp );
      but.addActionListener( listener );
      buttonPanel.add( but );
      buttonPanel.doLayout();
      buttonPanel.repaint();
      if ( buttons == null ) buttons = new HashMap();
      buttons.put( name, but );
      return but;
   }

   public JButton addButton( String name ) { return addButton( name, null ); }
   public JPanel getButtonPanel() { return buttonPanel; }
   public JButton getButton( String name ) { 
      return buttons == null ? null : (JButton) buttons.get( name ); }
   public JComponent getComponent() { return comp; }

   public JButton removeButton( String name ) { 
      if ( buttons == null ) return null;
      JButton but = (JButton) buttons.get( name );
      if ( but == null ) return but;
      buttonPanel.remove( but );
      buttonPanel.doLayout();
      buttonPanel.repaint();
      return but;
   }

   public void addStatus() { 
      status = new JLabel( "                   " );
      buttonPanel.add( status ); buttonPanel.doLayout(); 
   }
   
   public void setStatus( String stat ) { 
      if ( status == null ) addStatus(); 
      int width = status.getGraphics().getFontMetrics().stringWidth( stat );
      if ( status.size().width < width ) status.setSize( width, status.size().height );
      status.setText( stat ); 
   }   
   
   public void createUI( int width, int height ) {
      center( width, height + PANEL_HEIGHT );
      //if ( this instanceof JComponent ) ( (JComponent) this ).setDoubleBuffered( true ); 
      addWindowListener( this );
      //getContentPane().add( BorderLayout.CENTER, comp );
      buttonPanel = new JPanel();
      addButton( "Print", this );
      addButton( "Close", this );
      //JButton but = addButton( "Save PostScript", this );
      //but.hide();
      getContentPane().add( buttonPanel, BorderLayout.SOUTH );
      getContentPane().doLayout();
   }

   public void QuitOnClose() { 
      quitOnClose = true;
      setDefaultCloseOperation( EXIT_ON_CLOSE ); }

   public void windowOpened( WindowEvent e ) { }
   public void windowClosing( WindowEvent e ) { this.close(); }
   public void windowClosed( WindowEvent e ) { this.close(); }
   public void windowIconified( WindowEvent e ) { }
   public void windowDeiconified( WindowEvent e ) { }
   public void windowActivated( WindowEvent e ) { }
   public void windowDeactivated( WindowEvent e ) { }
  
   public void center( int w, int h ) {
      Dimension screenSize = Toolkit.getDefaultToolkit().getScreenSize();
      int x = (screenSize.width - w) / 2;
      int y = (screenSize.height - h) / 2;
      setBounds( x, y, w, h );
   }

   public void center() {
      Dimension screenSize = Toolkit.getDefaultToolkit().getScreenSize();
      int x = (screenSize.width - getWidth()) / 2;
      int y = (screenSize.height - getHeight()) / 2;
      setBounds( x, y, getWidth(), getHeight() );
   }
}
