package djr.motif.model.prior;
import java.util.*;
import djr.util.array.*;
import djr.motif.*;
import djr.motif.model.*;
import djr.util.MyUtils;

/**
 * Class <code>DirichletMixturePrior</code>
 *
 * @author <a href="mailto:reiss@uw.edu">David Reiss</a>
 * @version 1.9978 (Fri Nov 07 05:56:26 PST 2003)
 */
public class DirichletMixturePrior extends PseudocountPrior {
   protected static final double CACHE_EMPTY = Double.MIN_VALUE;
   protected static final int CACHE_SIZE = 1300 * 1000;
   protected static double gammaCache[] = DoubleUtils.New( CACHE_SIZE );
   static { DoubleUtils.Set( gammaCache, CACHE_EMPTY ); }

   protected int order;
   protected double alpha[][], q[], alphatot[], logBetaAlpha[], temp[];

   public DirichletMixturePrior( String fname, double pseudo, int alphSize ) {
      super( pseudo, alphSize );
      ReadDirichletData( fname );
   }   

   public DirichletMixturePrior( String fname, double pseudo, BackgroundModel bg ) {
      super( pseudo, bg );
      ReadDirichletData( fname );
   }
   
   public void GetProbs( double[] cts, double[] probs ) {
      double sum = 0.0, countstot = DoubleUtils.Sum( cts );
      for ( int i = 0; i < alphSize; i ++ ) {
	 double prob = 0.0, temp2;
	 double max = -Double.MAX_VALUE;
	 for ( int j = 0; j < order; j ++ ) {
	    if ( i == 0 ) temp[ j ] = LogBetaFunc2( alpha[ j ], alphatot[ j ], cts, countstot );
	    if ( temp[ j ] - logBetaAlpha[ j ] > max ) max = temp[ j ] - logBetaAlpha[ j ];
	 }
	 for ( int j = 0; j < order; j ++ ) {
	    double val = q[ j ] * ( alpha[ j ][ i ] + cts[ i ] ) / 
	       ( alphatot[ j ] + countstot );
	    double logb1 = temp[ j ], logb2 = logBetaAlpha[ j ];
	    prob += val * Math.exp( logb1 - logb2 - max );
	 }
	 sum += prob;
	 probs[ i ] = prob;
      }
      DoubleUtils.Divide( probs, sum );
   }

   /**
    * Compute the log-gamma function. Cache the value and
    * check if it has been cached before computing again.
    */
   protected static final double LogGamma( double xx ) {
      if ( xx < 0 ) xx = 0;
      //if ( xx < 0.001 ) return DoubleUtils.logGamma( xx );
      int ind = (int) ( xx * 1000 );
      if ( ind >= CACHE_SIZE ) return DoubleUtils.logGamma( xx );
      double out = gammaCache[ ind ];
      if ( out == CACHE_EMPTY ) {
	 out = DoubleUtils.logGamma( xx );
	 gammaCache[ ind ] = out;
      }
      return out;
   }

   protected static final double LogBetaFunc( double[] vec, double vectot ) {
      double sum = 0; 
      for ( int i = 0, size = vec.length; i < size; i ++ ) sum += LogGamma( vec[ i ] );
      sum -= LogGamma( vectot ); return sum;
   }

   protected static final double LogBetaFunc2( double[] vec1, double vectot1,
					       double[] vec2, double vectot2 ) {
      double sum = 0; 
      for ( int i = 0, size = vec1.length; i < size; i ++ ) 
	 sum += LogGamma( vec1[ i ] + vec2[ i ] );
      sum -= LogGamma( vectot1 + vectot2 ); return sum;
   }

   /**
    * Read in dirichlet mixture data file taken from Karplus web site 
    * http://www.soe.ucsc.edu/research/compbio/dirichlets/
    */
   protected void ReadDirichletData( String fname ) {
      try {
	 Enumeration lines = MyUtils.ReadFileLines( fname ).elements();
	 int line = 0;
	 while( lines.hasMoreElements() ) {
	    String str = (String) lines.nextElement();
	    if ( str == null || str.startsWith( "#" ) || str.startsWith( "//" ) || 
		 str.equals( "" ) ) continue;
	    if ( str.startsWith( "Alphabet" ) && ! str.endsWith( "ExtAA" ) ) {
	       throw new Exception( "Incorrect alphabet in dirichlet mixture model" );
	    } else if ( str.startsWith( "Order" ) && ! str.endsWith( "A C D E F G H I K L M N P Q R S T V W Y" ) ) {
	       throw new Exception ( "Incorrect alphabet order in dirichlet mixture model" );
	    } else if ( str.startsWith( "NumDistr" ) ) {
	       String tt[] = MyUtils.Tokenize( str, "= " );
	       order = Integer.parseInt( tt[ tt.length - 1 ] );
	       alpha = DoubleUtils.New( order, 20 );
	       q = DoubleUtils.New( 20 );
	       alphatot = DoubleUtils.New( 20 );
	       logBetaAlpha = DoubleUtils.New( 20 );
	       temp = DoubleUtils.New( order );
	    } else if ( str.startsWith( "Number" ) ) {
	       String tt[] = MyUtils.Tokenize( str, "= " );
	       int num = Integer.parseInt( tt[ tt.length - 1 ] );
	       str = (String) lines.nextElement();
	       if ( str.startsWith( "Mixture" ) ) {
		  tt = MyUtils.Tokenize( str, "= " );
		  q[ num ] = Double.parseDouble( tt[ tt.length - 1 ] );
	       }
	       str = (String) lines.nextElement();
	       if ( str.startsWith( "Alpha" ) ) {
		  double dd[] = 
		     DoubleUtils.FTokenize( str.substring( str.lastIndexOf( "=" ) + 1 ), " " );
		  alphatot[ num ] = dd[ 0 ];
		  System.arraycopy( dd, 1, alpha[ num ], 0, 20 );
		  logBetaAlpha[ num ] = LogBetaFunc( alpha[ num ], alphatot[ num ] );
	       }
	    }
	 }
      } catch( Exception e ) {
         System.err.println( "Could not load dirichlet mixtures file " + fname );
         e.printStackTrace();
	 MyUtils.Exit( 0 );
      }
   }
}

