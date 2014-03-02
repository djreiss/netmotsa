package djr.util.gui;
import java.awt.*;
import java.awt.print.*; 
import javax.print.*;
import javax.print.attribute.*;
import javax.print.attribute.standard.*;
import java.io.*;
 
import javax.swing.JComponent; 
 
/**
 * Class <code>ComponentPrintable</code>
 *
 * @author <a href="mailto:dreiss@systemsbiology.org">David Reiss</a>
 * @version 1.9978 (Fri Nov 07 05:56:26 PST 2003)
 */
public class ComponentPrintable implements Printable { 
   private Component mComponent; 
   
   public ComponentPrintable(Component c) { 
      mComponent = c; 
   } 
   
   public int print(Graphics g, PageFormat pageFormat, int pageIndex) { 
      if (pageIndex > 0) return NO_SUCH_PAGE; 
      Graphics2D g2 = (Graphics2D)g; 
      g2.translate( pageFormat.getImageableX(), pageFormat.getImageableY() ); 
      boolean wasBuffered = disableDoubleBuffering( mComponent ); 

      double xScale = pageFormat.getImageableWidth() / mComponent.getWidth();
      double yScale = pageFormat.getImageableHeight() / mComponent.getHeight();
      double scale = Math.min( xScale, yScale );
      g2.scale( scale, scale );

      mComponent.paint( g2 ); 
      restoreDoubleBuffering( mComponent, wasBuffered ); 
      return PAGE_EXISTS; 
   } 

   public void printToPostScriptFile( String fname ) {
      /* Use the pre-defined flavor for a Printable from an InputStream */
      DocFlavor flavor = DocFlavor.SERVICE_FORMATTED.PRINTABLE;

      /* Specify the type of the output stream */
      String psMimeType = DocFlavor.BYTE_ARRAY.POSTSCRIPT.getMimeType();

      /* Locate factory which can export a GIF image stream as Postscript */
      StreamPrintServiceFactory[] factories =
	 StreamPrintServiceFactory.lookupStreamPrintServiceFactories(flavor, psMimeType);
								    
      if (factories.length == 0) {
	 System.err.println("No suitable factories");
	 System.exit(0);
      }

      try {
	 /* Create a file for the exported postscript */
	 FileOutputStream fos = new FileOutputStream(fname);

	 /* Create a Stream printer for Postscript */
	 StreamPrintService sps = factories[0].getPrintService(fos);

	 /* Create and call a Print Job */
	 DocPrintJob pj = sps.createPrintJob();
	 PrintRequestAttributeSet aset = new HashPrintRequestAttributeSet();

	 Doc doc = new SimpleDoc(this, flavor, null);

	 pj.print(doc, aset);
	 fos.flush();
	 fos.close();

      } catch (PrintException pe) { 
	 System.err.println(pe);
      } catch (IOException ie) { 
	 System.err.println(ie);
      }
   }
 
   private boolean disableDoubleBuffering( Component c ) { 
      if ( c instanceof JComponent == false ) return false; 
      JComponent jc = (JComponent)c; 
      boolean wasBuffered = jc.isDoubleBuffered(); 
      jc.setDoubleBuffered( false ); 
      return wasBuffered; 
   } 
   
   private void restoreDoubleBuffering(Component c, boolean wasBuffered) { 
      if (c instanceof JComponent) 
	 ((JComponent)c).setDoubleBuffered(wasBuffered); 
   } 
}
