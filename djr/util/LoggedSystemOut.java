package djr.util;

import java.io.*;
import corejava.Format;

/**
 * Class <code>LoggedSystemOut</code>
 *
 * @author <a href="mailto:dreiss@systemsbiology.org">David Reiss</a>
 * @version 1.9978 (Fri Nov 07 05:56:26 PST 2003)
 */
public class LoggedSystemOut extends MyOut {
   protected int debug = 0;
   public PrintStream LOG;

   public LoggedSystemOut( String logFile ) {
      super( System.out );
      try {
	 LOG = new PrintStream( MyUtils.OpenOutputFile( logFile ), true );
      } catch( Exception e ) {
	 LOG = null;
	 e.printStackTrace();
      }
   }

   public void setDebugLevel( int lev ) {
      debug = lev;
   }

   public void write( byte b ) {
      super.write( b );
      if ( LOG != null ) LOG.write( b );
   }

   public void write( byte buf[], int off, int len ) {
      super.write( buf, off, len );
      if ( LOG != null ) LOG.write( buf, off, len );
   }

   public void close() {
      if ( LOG != null ) LOG.close();
      LOG = null;
      super.close();
   }

   protected void finalize() throws Throwable {
      close();
   }

   public void DEBUG( Object str ) { if ( debug > 0 ) System.err.println( str ); }
}
