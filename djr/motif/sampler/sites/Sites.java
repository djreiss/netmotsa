package djr.motif.sampler.sites;
import java.io.PrintStream;
import java.util.*;

import djr.motif.model.prior.*;
import djr.motif.sampler.*;
import djr.util.array.*;
import djr.util.bio.Sequence;

/**
 * Class <code>Sites</code>
 *
 * @author <a href="mailto:reiss@uw.edu">David Reiss</a>
 * @version 1.9978 (Fri Nov 07 05:56:26 PST 2003)
 */
public class Sites implements java.io.Serializable {
   Sequence S[];
   short seqs[][];
   int NS = 0, a[][], mots[][], nMotifs, minPerSeq, maxPerSeq, nsites[], W;
   int maxLen, sumLen, siteCounts[][], motifCounts[], len[];
   
   public Sites( Sequence S[], int W, int nMotifs ) {
      this( S, W, nMotifs, 0, 1 );
   }

   public Sites( Sequence S[], int W, int nMotifs, int minPerSeq, int maxPerSeq ) {
      this.W = W;
      this.nMotifs = nMotifs;
      SetSequences( S );
      siteCounts = IntUtils.New( nMotifs, maxLen );
      motifCounts = IntUtils.New( nMotifs + 1 );
      motifCounts[ nMotifs ] = sumLen;
      SetPerSequenceLimits( minPerSeq, maxPerSeq );
   }

   public Sites( Sites s ) {
      this( s.S, s.W, s.nMotifs, s.minPerSeq, s.maxPerSeq );
      CopyFrom( s );
   }

   public void SetSequences( Sequence S[] ) {
      this.S = S;
      if ( this.NS != S.length ) {
	 this.NS = S.length;
	 a = new int[ NS ][];
	 mots = new int[ NS ][];
	 seqs = new short[ NS ][];
	 len = IntUtils.New( NS );
	 maxLen = -1; sumLen = 0;
	 for ( int i = 0; i < NS; i ++ ) {
	    len[ i ] = S[ i ].GetLength();
	    if ( len[ i ] > maxLen ) maxLen = len[ i ];
	    sumLen += len[ i ];
	    seqs[ i ] = S[ i ].GetResidues();
	 }
	 nsites = IntUtils.New( NS );
      }
   }

   public void ClearCounts() {
      IntUtils.Zero( siteCounts );
      IntUtils.Zero( motifCounts );
      motifCounts[ nMotifs ] = sumLen;
   }

   public void ClearCounts( int mot ) {
      IntUtils.Zero( siteCounts );
      motifCounts[ nMotifs ] += motifCounts[ mot ];
      motifCounts[ mot ] = 0;
   }

   public void CopyFrom( Sites s ) {
      if ( this.a == null || this.NS != s.NS ) {
	 this.a = IntUtils.New( s.a );
	 this.mots = IntUtils.New( s.mots );
	 if ( s.len != null ) this.len = IntUtils.New( s.len );
	 if ( s.seqs != null ) this.seqs = ShortUtils.New( s.seqs );
	 this.nsites = IntUtils.New( s.nsites );
	 this.siteCounts = IntUtils.New( s.siteCounts );
	 this.motifCounts = IntUtils.New( s.motifCounts );
      } else {
	 IntUtils.Copy( this.a, s.a );
	 IntUtils.Copy( this.mots, s.mots );
	 if ( s.len != null ) IntUtils.Copy( this.len, s.len );
	 if ( s.seqs != null ) ShortUtils.Copy( this.seqs, s.seqs );
	 IntUtils.Copy( this.nsites, s.nsites );
	 IntUtils.Copy( this.siteCounts, s.siteCounts );
	 IntUtils.Copy( this.motifCounts, s.motifCounts );
      }
      if ( s.S != null ) this.S = s.S;
      this.NS = s.NS;
      this.W = s.W;
      this.nMotifs = nMotifs;
      this.minPerSeq = s.minPerSeq;
      this.maxPerSeq = s.maxPerSeq;
      this.maxLen = s.maxLen;
      this.sumLen = s.sumLen;
   }

   public int GetSequenceLength( int seqNum ) { return len[ seqNum ]; }
   public int[] GetSequenceLengths() { return len; }
   public int GetMaxSequenceLength() { return maxLen; }
   public Sequence[] GetSequences() { return S; }
   public Sequence GetSequence( int seqNum ) { return S[ seqNum ]; }
   public short[] GetSequenceResidues( int seqNum ) { return seqs[ seqNum ]; }
   public short[][] GetSequenceResidues() { return seqs; }

   public void SetPerSequenceLimits( int min, int max ) {
      for ( int i = 0; i < NS; i ++ ) {
	 if ( a[ i ] == null ) {
	    a[ i ] = IntUtils.New( max );
	    IntUtils.Set( a[ i ], -1 );
	    mots[ i ] = IntUtils.New( max );
	    IntUtils.Set( mots[ i ], -1 );
	 } else if ( a.length < max ) {
	    a[ i ] = IntUtils.Resize( a[ i ], max );
	    IntUtils.Set( a[ i ], maxPerSeq, a[ i ].length - 1, -1 );
	    mots[ i ] = IntUtils.Resize( mots[ i ], max );
	    IntUtils.Set( mots[ i ], maxPerSeq, mots[ i ].length - 1, -1 );
	 }
      }
      minPerSeq = min;
      maxPerSeq = max;
   }

   public void Initialize() {
      IntUtils.Set( a, -1 );
      IntUtils.Set( mots, -1 );
      IntUtils.Zero( nsites );
   }

   public int GetSequenceCount() {
      return NS;
   }

   public int GetMotForSite( int seq, int site ) {
      int aa[] = a[ seq ], ns = nsites[ seq ];
      for ( int i = 0; i < ns; i ++ ) if ( aa[ i ] == site ) return mots[ seq ][ i ];
      return -1;
   }

   public int GetSiteForMot( int seq, int mot ) {
      int mm[] = mots[ seq ], ns = nsites[ seq ];
      for ( int i = 0; i < ns; i ++ ) if ( mm[ i ] == mot ) return a[ seq ][ i ];
      return -1;
   }

   public boolean IsSiteValid( int seq, int site ) {
      if ( S[ seq ].IsMaskedInRange( site - W, site + W ) ) return false;
      int aa[] = a[ seq ], mm[] = mots[ seq ], ns = nsites[ seq ];
      int mot = GetMotForSite( seq, site );
      for ( int i = 0; i < ns; i ++ ) {
	 int othersite = aa[ i ];
	 if ( othersite == site ) return false;
	 if ( site > othersite - W && site < othersite + W ) return false;
      }
      return true;
   }

   public void ShuffleSites() {
      ClearCounts();
      for ( int i = 0; i < NS; i ++ ) ShuffleSitesForSequence( i );
   }

   public boolean ShuffleSitesForSequence( int seq ) {
      IntUtils.Set( a[ seq ], -1 );
      IntUtils.Set( mots[ seq ], -1 );
      nsites[ seq ] = 0;
      int nmots = IntUtils.Random( minPerSeq, maxPerSeq ) + 1;
      short res[] = S[ seq ].GetResidues();
      for ( int i = 0; i < nmots; i ++ ) {
	 int site = IntUtils.Random( len[ seq ] - W );
	 while( res[ site ] < 0 )
	    site = IntUtils.Random( len[ seq ] - W );
	 AddSite( seq, site, IntUtils.Random( nMotifs ) );
      }
      return true;
   }

   public boolean HasSite( int seq, int site ) {
      int ns = nsites[ seq ], aa[] = a[ seq ];
      for ( int i = 0; i < ns; i ++ ) if ( aa[ i ] == site ) return true;
      return false;
   }

   public boolean HasMot( int seq, int mot ) {
      int ns = nsites[ seq ], mm[] = mots[ seq ];
      for ( int i = 0; i < ns; i ++ ) if ( mm[ i ] == mot ) return true;
      return false;
   }

   public boolean HasSite( int seq, int site, int mot ) {
      int ns = nsites[ seq ], aa[] = a[ seq ], mm[] = mots[ seq ];
      for ( int i = 0; i < ns; i ++ ) if ( aa[ i ] == site && mm[ i ] == mot ) return true;
      return false;
   }

   public boolean RemoveSiteAtIndex( int seq, int ind ) {
      int ns = nsites[ seq ];
      if ( ind > ns ) return false;
      int sites[] = a[ seq ], mm[] = mots[ seq ];
      int mot = mm[ ind ], site = sites[ ind ];
      sites[ ind ] = sites[ ns - 1 ];
      sites[ ns - 1 ] = -1;
      mm[ ind ] = mm[ ns - 1 ];
      mm[ ns - 1 ] = -1;
      nsites[ seq ] --;
      if ( siteCounts[ mot ][ site ] > 0 ) siteCounts[ mot ][ site ] --;
      if ( motifCounts[ mot ] > 0 ) { motifCounts[ mot ] --; motifCounts[ nMotifs ] ++; }
      return true;
   }

   public void RemoveAllSites( int seq ) {
      while ( nsites[ seq ] > 0 ) RemoveSiteAtIndex( seq, 0 );
   }

   public void RemoveAllSites( int seq, int mot ) {
      int mm[] = mots[ seq ], aa[] = a[ seq ], ns = nsites[ seq ], ind = 0;
      while( ind >= 0 ) {
	 ind = IntUtils.IndexOf( mm, mot );
	 if ( ind >= 0 ) RemoveSiteAtIndex( seq, ind );
      }
   }

   public boolean RemoveSite( int seq, int site ) {
      if ( ! HasSite( seq, site ) ) return false;
      int sites[] = a[ seq ];
      int ind = IntUtils.IndexOf( sites, site );
      return RemoveSiteAtIndex( seq, ind );
   }

   public boolean RemoveSite( int seq, int site, int mot ) {
      if ( ! HasSite( seq, site, mot ) ) return false;
      int sites[] = a[ seq ], mm[] = mots[ seq ], ns = nsites[ seq ];
      int ind = 0;
      while( ind >= 0 ) {
	 ind = IntUtils.IndexOf( sites, site );
	 if ( mm[ ind ] == mot ) RemoveSiteAtIndex( seq, ind );
      }
      return true;
   }

   public boolean AddSite( int seq, int site, int mot ) {
      return AddSite( seq, site, mot, false );
   }

   public boolean AddSite( int seq, int site, int mot, boolean force ) {
      int ns = nsites[ seq ];
      if ( ns >= maxPerSeq || ( ! force &&
				( ! IsSiteValid( seq, site ) || HasSite( seq, site ) ) ) ) return false;
      if ( HasSite( seq, site, mot ) ) return false;
      a[ seq ][ ns ] = site;
      mots[ seq ][ nsites[ seq ] ++ ] = mot;
      siteCounts[ mot ][ site ] ++;
      motifCounts[ mot ] ++; motifCounts[ nMotifs ] --;
      return true;
   }

   public int[] GetSites( int seqNum ) {
      return a[ seqNum ];
   }

   public int[][] GetSites() {
      return a;
   }

   public int[] GetMots( int seqNum ) {
      return mots[ seqNum ];
   }

   public int[][] GetMots() {
      return mots;
   }

   public int[] GetNSites() {
      return nsites;
   }

   public int GetNSites( int seqNum ) {
      return nsites[ seqNum ];
   }

   public int GetSiteCounts( int seqNum, int mot ) {
      return siteCounts[ mot ][ seqNum ];
   }

   public int[] GetSiteCounts( int mot ) {
      return siteCounts[ mot ];
   }

   public int[][] GetSiteCounts() {
      return siteCounts;
   }

   public int GetMotifCounts( int mot ) {
      return motifCounts[ mot ];
   }

   public int[] GetMotifCounts() {
      return motifCounts;
   }

   public void SetMotifCounts( int mc[] ) {
      this.motifCounts = mc;
   }

   public int GetNMotifs() { 
      return nMotifs;
   }

   public int GetSumLen() {
      return sumLen;
   }

   /** Format of output alignments aa[][][] is mot x num x (alignment length + 1) where
       num is the number of alignments attributed to motif number mot in all sequences.
       The alignment is in the last dimension. First index is the seq. number, and indices
       1-length give the indexes in that sequence that are part of the alignment. (Sorry!) **/
   public int[][][] GetAlignments() {
      int out[][][] = new int[ nMotifs ][][];
      for ( int mot = 0; mot < nMotifs; mot ++ ) out[ mot ] = GetAlignments( mot );
      return out;
   }   

   public int[][] GetAlignments( int mot ) {
      ObjVector tempout = new ObjVector( 10, 10 );
      for ( int i = 0; i < NS; i ++ ) GetAlignmentsForSequence( i, mot, tempout );
      int out[][] = new int[ tempout.size() ][];
      for ( int i = 0; i < tempout.size(); i ++ ) out[ i ] = (int[]) tempout.elementAt( i );
      return out;
   }

   public void GetAlignmentsForSequence( int seqNum, ObjVector out ) {
      int sites[] = a[ seqNum ], ns = nsites[ seqNum ];
      for ( int i = 0; i < ns; i ++ ) 
	 out.addElement( GetAlignmentForSequenceAndSite( seqNum, sites[ i ], null ) );
   }

   public void GetAlignmentsForSequence( int seqNum, int mot, ObjVector out ) {
      int sites[] = a[ seqNum ], ns = nsites[ seqNum ], mm[] = mots[ seqNum ];
      for ( int i = 0; i < ns; i ++ ) if ( mm[ i ] == mot )
	 out.addElement( GetAlignmentForSequenceAndSite( seqNum, mot, sites[ i ], null ) );
   }

   protected int[] GetAlignmentForSequenceAndSite( int seqNum, int site, int out[] ) {
      int llen = len[ seqNum ], mot = GetMotForSite( seqNum, site );
      if ( out == null || out.length < W + 1 ) out = IntUtils.New( W + 1 );
      out[ 0 ] = seqNum;
      for ( int j = 0; j < W; j ++ ) if ( site + j < llen ) out[ j+1 ] = site + j;
      return out;
   }

   protected int[] GetAlignmentForSequenceAndSite( int seqNum, int mot, int site, int out[] ) {
      int llen = len[ seqNum ];
      if ( out == null || out.length < W + 1 ) out = IntUtils.New( W + 1 );
      out[ 0 ] = seqNum;
      for ( int j = 0; j < W; j ++ ) if ( site + j < llen ) out[ j+1 ] = site + j;
      return out;
   }

   public String toString() {
      String out = "";
      for ( int seq = 0; seq < NS; seq ++ ) {
	 int ns = nsites[ seq ], aa[] = a[ seq ], mm[] = mots[ seq ];
	 out += "seq=" + seq + ": ";
	 for ( int i = 0; i < ns; i ++ ) out += mm[ i ] + " at " + aa[ i ] + ";   ";
	 out += "\n";
      }
      return out;
   }

   public double TryRewindowing( double G ) {
      return G;
   }

   public void ShiftSites( int amount ) {
      for ( int i = 0; i < NS; i ++ ) ShiftSitesForSequence( i, amount ); }

   protected void ShiftSitesForSequence( int seqNum, int amount ) {
      if ( nsites[ seqNum ] <= 0 ) return;
      int sites[] = a[ seqNum ];
      for ( int i = 0, ns = nsites[ seqNum ]; i < ns; i ++ ) {
	 sites[ i ] += amount;
	 if ( sites[ i ] < 0 ) sites[ i ] = 0;
	 else if ( sites[ i ] > len[ seqNum ] - W - 1 ) sites[ i ] = len[ seqNum ] - W - 1;
      }
   }

   public Map LocateAlignments( Sampler samp, int minPerSeq, int maxPerSeq, double alignThresh ) {
      double scores[][] = new double[ NS ][];
      int out[][][] = new int[ nMotifs ][][];
      double pvalues[][] = new double[ nMotifs ][];
      boolean flags[] = null;
      for ( int mot = 0; mot < nMotifs; mot ++ ) {
	 double distrib[] = samp.GetMotifModel().ComputePolyScoreDistribution( mot );
	 DoubleUtils.UnitNorm( distrib );
	 double minScore = 0, thresh = 1.0 - Math.pow( 10.0, -alignThresh ), ssum = 0;
	 for ( int i = 0; i < distrib.length; i ++ ) {
	    ssum += distrib[ i ];
	    if ( ssum >= thresh ) { minScore = (double) i; break; }
	 }

	 ObjVector tempout = new ObjVector();
	 DoubleVector temppval = new DoubleVector();
	 double cc[][] = samp.GetMotifModel().GetMotifCounts( mot );
	 double mean = DoubleUtils.Mean( cc ), bestScores[] = DoubleUtils.New( maxPerSeq ); 
	 int bestTemp[][] = new int[ maxPerSeq ][];
	 for ( int i = 0; i < NS; i ++ ) {
	    if ( scores[ i ] == null ) scores[ i ] = DoubleUtils.New( len[ i ] );
	    else DoubleUtils.Zero( scores[ i ] );
	    double sscores[] = scores[ i ];
	    if ( flags == null || flags.length < len[ i ] ) flags = BoolUtils.New( len[ i ] );
	    else Arrays.fill( flags, false );
	    samp.ComputeProbsForSequence( i, mot, sscores, false );
	    int indexx[] = DoubleUtils.Indexx( sscores, null );
	    IntUtils.Reverse( indexx );
	    short[] seq = seqs[ i ];
	    if ( maxPerSeq > 0 ) {
	       DoubleUtils.Set( bestScores, -Double.MAX_VALUE );
	       Arrays.fill( bestTemp, null );
	    }
	    int found = 0;
	    for ( int j = 0; j < len[ i ]; j ++ ) {
	       int temp[] = null;
	       temp = GetAlignmentForSequenceAndSite( i, mot, indexx[ j ], temp );
	       double score = 0.0;
	       for ( int ii = 1; ii < temp.length; ii ++ ) {
		  int site = temp[ ii ];
		  if ( site >= 0 ) score += cc[ ii-1 ][ seq[ site ] ];
		  else score -= mean; // Penalize gaps by subtracting the mean value in c[][].
	       }
	       if ( score >= minScore ) {
		  for ( int ii = 1; ii < temp.length; ii ++ ) {
		     int site = temp[ ii ];
		     if ( site >= 0 && flags[ site ] ) score -= cc[ ii-1 ][ seq[ site ] ];
		  }
	       }
	       if ( score >= minScore ) {
		  for ( int ii = 1; ii < temp.length; ii ++ ) {
		     int site = temp[ ii ];
		     if ( site >= 0 ) flags[ site ] = true;
		  }
		  tempout.addElement( temp );
		  double sum = DoubleUtils.Sum( distrib, 0, (int) score );
		  double pvalue = 1.0 - sum; if ( pvalue < 0.0 ) pvalue = 0.0;
		  temppval.addElement( pvalue );
		  found ++;
	       } else if ( maxPerSeq > 0 ) {
		  int min = DoubleUtils.WhereMin( bestScores );
		  if ( score > bestScores[ min ] ) {
		     bestTemp[ min ] = temp;
		     bestScores[ min ] = score;
		  }
	       }
	    }
	    if ( found >= minPerSeq ) {
	       if ( found > maxPerSeq ) found = maxPerSeq;
	       int iindexx[] = DoubleUtils.Indexx( bestScores, null );
	       IntUtils.Reverse( iindexx );
	       for ( int ss = 0; ss < /*minPerSeq -*/ found; ss ++ ) {
		  int iind = iindexx[ ss ];
		  tempout.addElement( bestTemp[ iind ] );
		  double sum = DoubleUtils.Sum( distrib, 0, (int) bestScores[ iind ] );
		  double pvalue = 1.0 - sum; if ( pvalue < 0.0 ) pvalue = 0.0;
		  temppval.addElement( pvalue );
	       }
	    }
	 }
	 if ( tempout.size() < 2 ) { out[ mot ] = null; pvalues[ mot ] = null; continue; }
	 out[ mot ] = new int[ tempout.size() ][];
	 pvalues[ mot ] = DoubleUtils.New( tempout.size() );
	 for ( int i = 0; i < tempout.size(); i ++ ) {
	    out[ mot ][ i ] = (int[]) tempout.elementAt( i );
	    pvalues[ mot ][ i ] = temppval.elementAt( i );
	 }
      }
      //System.gc();
      Map map = new java.util.HashMap();
      map.put( "a", out );
      map.put( "pvalues", pvalues );
      return map;
   }

   public void MaskAlignment( Sequence seqs[], int aa[][][] ) {
      if ( aa == null ) return;
      for ( int mot = 0; mot < nMotifs; mot ++ ) {
	 if ( aa[ mot ] == null ) continue;
	 for ( int i = 0; i < aa[ mot ].length; i ++ ) {
	    int seq = aa[ mot ][ i ][ 0 ];
	    if ( seq < 0 || seqs[ seq ] == null ) continue;
	    for ( int j = 1; j < aa[ mot ][ i ].length - 1; j ++ ) { 
	       int loc = aa[ mot ][ i ][ j ];
	       if ( loc >= 0 ) seqs[ seq ].SetMask( loc );
	    }
	 }
      }
   }
}
