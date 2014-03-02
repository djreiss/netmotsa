package djr.motif.sampler;

import djr.util.*;
import djr.util.array.*;

/**
 * Class <code>SiteSampler</code>.
 *
 * @author <a href="mailto:dreiss@systemsbiology.org">David Reiss</a>
 * @version 1.9978 (Fri Nov 07 05:56:26 PST 2003)
 */
public class SiteSampler extends Sampler {
   public SiteSampler() { }
   public SiteSampler( String argv[] ) { Initialize( argv ); }
   public SiteSampler( String args ) { this( MyUtils.Tokenize( args, " " ) ); }

   protected void IterateSampler( int niter ) {
      int itermax = niter == -1 ? NS : niter;
      for ( int i = 0; i < itermax; i ++ ) {
	 int ii = niter == -1 ? i : IntUtils.RandChoose( 0, NS-1 );
	 AddRemoveFromCounts( false, ii, 0 );
	 int max = ComputeProbsForSequence( ii, 0, qx, true );
	 if ( ! noSamp || ( max >= 0 && qx[ max ] <= 0 ) ) {
	    int mmax = DoubleUtils.Sample( qx, len[ ii ] - W + 1 );
	    if ( mmax != -1 ) max = mmax;
	 }
	 //if ( max < 0 ) max = IntUtils.Random( maxl );
	 //while( ! IsSiteValid( max, ii ) ) max = IntUtils.Random( maxl );
	 if ( max != -1 && IsSiteValid( max, ii ) ) {
	    GetSites().RemoveAllSites( ii );
	    GetSites().AddSite( ii, max, 0 );
	 }
	 AddRemoveFromCounts( true, ii, 0 );
      }
   }

   public static void main( String argv[] ) {
      SiteSampler samp = new SiteSampler( argv );
      if ( ! samp.IsInitialized() ) return;
      samp.Run();
   }
}
