package djr.motif.sh3;

import java.io.*;
import java.util.*;
import cern.jet.random.*;

import djr.util.*;
import djr.util.bio.*;
import djr.util.array.*;
import djr.motif.sampler.*;
import djr.motif.model.prior.*;
import djr.motif.sampler.sites.*;
import djr.motif.model.*;

/**
 * Class <code>NetMultiSampler</code>
 *
 * @author <a href="mailto:reiss@uw.edu">David Reiss</a>
 * @version 1.9978 (Fri Nov 07 05:56:26 PST 2003)
 */
public class NetMultiSampler extends NetSiteSampler {
   public NetMultiSampler() { };
   public NetMultiSampler( String argv[] ) { Initialize( argv ); }
   public NetMultiSampler( String args ) { this( MyUtils.Tokenize( args, " " ) ); }

   protected boolean Initialize( String argv[] ) {
      super.Initialize( argv );

      int maxPtrs = -1;
      for ( int i = 0; i < nInts; i ++ ) 
	 if ( sh3Ptrs[ i ] != null && sh3Ptrs[ i ].length > maxPtrs )
	    maxPtrs = sh3Ptrs[ i ].length;
      Sites sites = new Sites( S, W, nMotifs, 1, maxPtrs + 1 );

      double pseudo2 = argProc.getFloatArg( "P2" );
      double pseudo3 = argProc.getFloatArg( "P3" );
      mmodel = new MultiMixtureMotifModel( nMotifs, sites, W, intPtrs, sh3Ptrs, 
					   pseudo2, pseudo3 );
      mmodel.SetBackgroundAndForegroundModels( bgfile, bgorder, pseudo, fgType, 
					       dirichName, matrix );

      initialized = true;
      return initialized;
   }

   protected void IterateSampler( int niter ) {
      int itermax = niter == -1 ? NS : niter;
      MultiMixtureMotifModel mmod = (MultiMixtureMotifModel) mmodel;
      for ( int i = 0; i < itermax; i ++ ) {
	 int seqNum = niter == -1 ? i : IntUtils.RandChoose( 0, NS-1 );

	 int sh3p[] = sh3Ptrs[ seqNum ];
	 if ( sh3p == null ) continue;
	 for ( int j = 0, nsh3p = sh3p.length; j < nsh3p; j ++ ) {
	    int mot = sh3p[ j ];
	    mmod.AddRemoveFromCounts( false, seqNum, mot );
	    int max = ComputeProbsForSequence( seqNum, mot, qx, true );
	    if ( ! noSamp ) {
	       int mmax = DoubleUtils.Sample( qx, len[ seqNum ] - W + 1 );
	       if ( mmax != -1 ) max = mmax;
	    }
	    if ( max != -1 /*&& IsSiteValid( max, seqNum )*/ ) {
	       GetSites().RemoveAllSites( seqNum, mot );
	       GetSites().AddSite( seqNum, max, mot, true );
	    }
	    mmod.AddRemoveFromCounts( true, seqNum, mot );
	 }
      }
   }   

   public int ComputeProbsForSequence( int seqNum, int mot, double qx[], boolean normalize ) {
      MultiMixtureMotifModel mmod = (MultiMixtureMotifModel) mmodel;
      int max = mmod.ComputeProbsForSequence( seqNum, mot, qx, normalize );
      if ( power != 1.0 ) DoubleUtils.Pow( qx, power );
      if ( doDiscrimination() ) max = AdjustProbsUsingDiscrimination( qx, seqNum, normalize, max );
      return max;
   }

   protected void ShuffleSites() {
      for ( int ii = 0; ii < NS; ii ++ ) ShuffleSitesForSequence( ii );
   }

   protected boolean ShuffleSitesForSequence( int seq ) {
      int partners[] = sh3Ptrs[ seq ];
      if ( partners == null ) return false;
      for ( int mmot = 0; mmot < nMotifs; mmot ++ ) {
	 if ( ! IntUtils.Contains( partners, mmot ) ) continue;
	 GetSites().RemoveAllSites( seq, mmot );
	 GetSites().AddSite( seq, IntUtils.Random( len[ seq ] - W - 1 ), mmot, true );
      }
      return true;
   }

   public int[][][] GetAlignments() {
      int out[][][] = new int[ nMotifs ][][];
      for ( int mot = 0; mot < nMotifs; mot ++ ) {
         ObjVector tempout = new ObjVector();
         for ( int i = 0; i < NS; i ++ ) {
            int partners[] = sh3Ptrs[ i ];
            if ( partners == null || ! IntUtils.Contains( partners, mot ) ) continue;
            GetSites().GetAlignmentsForSequence( i, mot, tempout );
         }
         if ( tempout.size() > 0 ) out[ mot ] = new int[ tempout.size() ][];
         for ( int i = 0; i < tempout.size(); i ++ )
            out[ mot ][ i ] = (int[]) tempout.elementAt( i );
      }
      return out;
   }

   public void SetupArgs( ArgProcessor argProc ) {
      super.SetupArgs( argProc );
      argProc.AddArg( "NetMultiSampler Parameters" );
      argProc.AddArg( "P3", "<frange:0:100>", "0.5", "mixture model pseudocount fraction" );
   }

   public static void main( String argv[] ) {
      NetMultiSampler samp = new NetMultiSampler( argv );
      if ( ! samp.IsInitialized() ) return;
      samp.Run();
   }
}
