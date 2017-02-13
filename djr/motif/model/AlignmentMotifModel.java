package djr.motif.model;
import java.util.Map;

import djr.motif.model.prior.*;
import djr.motif.sampler.sites.*;
import djr.motif.sampler.*;
import djr.util.*;
import djr.util.array.*;
import djr.util.bio.Sequence;

/**
 * Describe class <code>AlignmentMotifModel</code> here.
 *
 * @author <a href="mailto:reiss@uw.edu">David Reiss</a>
 * @version 1.9978 (Fri Nov 07 05:56:26 PST 2003)
 */
public class AlignmentMotifModel extends PSSMMotifModel {
   static Map bgModelCache = new java.util.HashMap();

   Sites sites;
   double c[][][];
   MotifModelPrior fgmodel;

   public AlignmentMotifModel() { };

   public AlignmentMotifModel( int nMotifs, Sites sites, int W ) {
      super( nMotifs, sites.GetSequence( 0 ).GetAlphabetSize(), W );
      this.sites = new Sites( sites );
      c = DoubleUtils.New( q );
   }

   public AlignmentMotifModel( AlignmentMotifModel m ) {
      this( m.nMotifs, m.sites, m.W );
      CopyFrom( m );
   }

   public PSSMMotifModel Duplicate() {
      return new AlignmentMotifModel( this );
   }

   public double[][] GetMotifCounts( int mot ) { return c[ mot ]; }
   public MotifModelPrior GetMotifModelPrior() { return fgmodel; }
   public Sites GetSites() { return sites; }

   public void SetMotifModels( double qq[][][] ) {
      super.SetMotifModels( qq );
      this.c = DoubleUtils.New( qq );
   }

   public void SetBackgroundAndForegroundModels( String bgfile, int bgorder, double pseudo,
						 String fgType, String dirichName, 
						 String matrix ) {
      if ( ! MyUtils.IsNullString( bgfile ) ) {
	 String key = bgfile + "_" + bgorder + "_" + pseudo;
	 if ( ! bgModelCache.containsKey( key ) ) {
	    bgmodel = new BackgroundModel( bgfile, (short) bgorder, pseudo );
	    bgModelCache.put( key, bgmodel );
	 } else bgmodel = (BackgroundModel) bgModelCache.get( key );
      }
      else bgmodel = new BackgroundModel( sites.GetSequences(), (short) bgorder, pseudo );
      fgmodel = MotifModelPrior.GetPrior( (short) J, fgType, dirichName, matrix, pseudo, 
					  bgmodel );
   }

   public void SetMotifModelPrior( MotifModelPrior fg ) {
      this.fgmodel = fg;
   }

   public void SetMotifCounts( double cc[][][] ) {
      this.c = cc;
      changed = new boolean[ cc.length ];
   }

   public void Initialize() {
      super.Initialize();
      DoubleUtils.Zero( c );
   }

   public void Copy( PSSMMotifModel to ) {
      super.Copy( to );
      AlignmentMotifModel ato = (AlignmentMotifModel) to;
      ato.sites = new Sites( sites );
      DoubleUtils.Copy( ato.c, c );
      ato.fgmodel = fgmodel;
   }

   public void SetSites( Sites s ) {
      this.sites = s;
   }

   public double ComputeMAPScore() {
      double F = 0;
      for ( int mot = 0; mot < nMotifs; mot ++ ) {
	 double qq[][] = GetMotifModel( mot ), cc[][] = GetMotifCounts( mot );
	 for ( int i = 0; i < W; i ++ ) {
	    double qi[] = qq[ i ], ci[] = cc[ i ], FF = 0.0;
	    for ( short j = 0; j < J; j ++ ) {
	       if ( qi[ j ] <= 0 || bgmodel.GetProb( j ) == 0 ) continue;
	       FF += ci[ j ] * DoubleUtils.Log2( qi[ j ] / bgmodel.GetProb( j ) );
	    }
	    if ( FF <= 0 || Double.isNaN( FF ) || Double.isInfinite( FF ) ) continue; 
	    F += DoubleUtils.Log2( FF * J );
	 }
      }
      if ( F == 0 || Double.isNaN( F ) || Double.isInfinite( F ) ) return -Double.MAX_VALUE;
      return F;
   }

   public void Print( Sampler samp, boolean ansi, boolean html, boolean printingQ,
                      String names[] ) {
      double save[][][] = q;
      if ( ! printingQ ) { q = c; FORMAT = " %4d"; }
      else FORMAT = " %5.3f";
      Print( samp.GetOutput(), ansi, html, names );
      //if ( samp.GetLog() != null ) Print( samp.GetLog(), ansi, html, names );
      q = save; save = null;
      if ( ! printingQ ) FORMAT = " %5.3f";
   }

   public void FillTheArrays() {
      int NS = sites.GetSequenceCount();
      DoubleUtils.Zero( c );
      for ( int mot = 0; mot < nMotifs; mot ++ )
	 for ( int i = 0; i < NS; i ++ ) AddRemoveFromCounts( true, i, mot );
   }

   public void AddRemoveFromCounts( boolean add, int seqNum, int mot ) {
      int locs[] = sites.GetSites( seqNum ), mots[] = sites.GetMots( seqNum );
      int ns = sites.GetNSites( seqNum );
      for ( int i = 0; i < ns; i ++ ) {
         if ( locs[ i ] < 0 || mots[ i ] < 0 ) continue;
	 if ( mot >= 0 && mots[ i ] != mot ) continue;
         AddRemoveFromCounts( add, seqNum, locs[ i ], mots[ i ] );
      }
   }

   public void AddRemoveFromCounts( boolean add, int seqNum, 
				    int loc, int mot ) {
      boolean hasIt = sites.HasSite( seqNum, loc, mot );
      if ( ! add && ! hasIt ) return;
      AddRemoveFromCounts( add, sites.GetSequenceResidues( seqNum ), loc,
			   c[ mot ], 1.0 );
      changed[ mot ] = true;
   }

   public void AddRemoveFromCounts( boolean add, short seq[], int loc,
				    double counts[][], double amt ) {
      double ad = add ? amt : -amt;
      int llen = seq.length, WW = counts.length;
      for ( int i = 0; i < WW; i ++ ) {
         if ( loc + i >= llen ) break;
         int j = (int) seq[ loc + i ];
	 if ( j < 0 ) break; // Masked
         double[] ci = counts[ i ];
         ci[ j ] += ad;
         if ( ci[ j ] < 0 ) ci[ j ] = 0;
      }
   }

   public int ComputeProbsForSequence( int seqNum, int mot, double[] qx, 
				       boolean normalize ) {
      return ComputeProbsForSequence( sites.GetSequenceResidues( seqNum ), GetMotifModel( mot ),
				      qx, normalize );
   }

   public double ComputeProbForSite( int seqNum, int site, int mot ) {
      return ComputeProbForSite( sites.GetSequenceResidues( seqNum ), site, GetMotifModel( mot ) );
   }
      
   public void FillQ( double cc[][], double qq[][] ) {
      for ( int i = 0, WWW = cc.length; i < WWW; i ++ ) 
	 fgmodel.GetProbs( cc[ i ], qq[ i ] );
   }

   public void FillQColumn( int mot, int col ) {
      fgmodel.GetProbs( c[ mot ][ col ], q[ mot ][ col ] );
   }

   protected void RefillFromAlignment( int aa[][][], double pvalues[][] ) {
      for ( int mot = 0; mot < nMotifs; mot ++ ) {
	 if ( aa[ mot ] == null ) continue;
	 int ns = aa[ mot ].length;
	 if ( ns <= 0 ) continue;
	 double cc[][] = c[ mot ];
	 DoubleUtils.Zero( cc );
	 for ( int i = 0; i < ns; i ++ ) {
	    int seqNum = aa[ mot ][ i ][ 0 ];
	    if ( seqNum < 0 || sites.GetSequence( seqNum ) == null ) continue;
	    short seq[] = sites.GetSequenceResidues( seqNum );
	    for ( int j = 1, size = aa[ mot ][ i ].length; j < size; j ++ ) {
	       int loc = aa[ mot ][ i ][ j ];
	       if ( loc >= 0 ) {
		  int res = (int) seq[ loc ];
		  cc[ j-1 ][ res ] ++;
	       } else {
		  for ( int jj = 0; jj < J; jj ++ ) 
		     cc[ j-1 ][ jj ] += bgmodel.GetProb( jj );
	       }
	    }
	 }
	 FillQ( mot );
      }
   }

   public double ComputeMatchScore( int aa[][][], int ii, int mot ) {
      int aaa[][] = aa[ mot ], seqNum = aaa[ ii ][ 0 ];
      double fg = 1.0, bg = 1.0, qmot[][] = GetMotifModel( mot );
      if ( qmot == null ) return 0.0;
      if ( seqNum < 0 ) return 0.0;
      short[] seq = sites.GetSequenceResidues( seqNum );
      int ww = aaa[ ii ].length - 1;
      for ( int i = 0; i < ww; i ++ ) {
	 int site = aaa[ ii ][ i + 1 ];
	 if ( site >= 0 && seq[ site ] >= 0 ) {
	    fg *= qmot[ i ][ seq[ site ] ];
	    bg *= bgmodel.GetProb( seq, site );
	 }
      }
      return DoubleUtils.Log2( fg / bg );
   }

   public double[] ComputeEValues( int aa[][][] ) {
      double out[] = DoubleUtils.New( nMotifs ); 
      for ( int mot = 0; mot < nMotifs; mot ++ ) out[ mot ] = ComputeEValue( aa, mot );
      return out;
   }

   public double ComputeEValue( int aa[][][], int mot ) {
      if ( aa[ mot ] == null ) return 0.0;
      double out = 0.0, qq[][] = GetMotifModel( mot );
      for ( int ii = 0; ii < aa[ mot ].length; ii ++ ) 
	 //out += ComputeMatchScore( samp, aa, ii, mot, q1 );
	 out += ComputeMatchScore( aa, ii, mot );
      return out;
   }

   public static double[] ComputePolyScoreDistribution( double cc[][], BackgroundModel bgmodel ) {
      int max = (int) Math.ceil( DoubleUtils.Max( cc ) ) + 1, WW = cc.length, J = cc[ 0 ].length;
      if ( max < 1 ) max = 1;
      double polya[] = DoubleUtils.New( max * WW + max ), 
	 polyb[] = DoubleUtils.New( max * WW + max ),
	 polyc[] = DoubleUtils.New( max * WW + max );
      if ( max == 1 ) return polyc;
      for ( short j = 0; j < J; j ++ ) polya[ (int) cc[ 0 ][ j ] ] += bgmodel.GetProb( j );
      int a = max, b = max, c = 0;
      for ( int i = 1; i < WW; i ++ ) {
	 DoubleUtils.Zero( polyb );
	 for ( short j = 0; j < J; j ++ ) polyb[ (int) cc[ i ][ j ] ] += bgmodel.GetProb( j );
	 c = a + b;
	 polyMul( polya, polyb, polyc, a, b, c );
	 a = c;
      }
      //( new MyPlot( "Scores", "score", "probability" ) ).SetConnected( false ).
      // SetLog( false, false ).SetMarks( "points" ).PlotXY( polyc );
      return polyc;
   }

   public double[] ComputePolyScoreDistribution( int mot ) {
      double cc[][] = c[ mot ];
      return ComputePolyScoreDistribution( cc, bgmodel );
   }

   protected static final void polyMul( double polya[], double polyb[], double polyc[],
					int a, int b, int c ) {
      DoubleUtils.Zero( polyc );
      for ( int i = 0; i <= a; i ++ )
	 for ( int j = 0; j <= b; j ++ )
	    polyc[ i + j ] += polya[ i ] * polyb[ j ];
      for ( int i = 0; i <= c; i ++ ) {
	 double z = polyc[ i ];
	 if ( z < 1e-10 ) z = 0.0;
	 polya[ i ] = z;
      }
   }

   public double[][] ComputePolyScores( int aa[][][] ) {
      double out[][] = new double[ nMotifs ][];
      for ( int mot = 0; mot < nMotifs; mot ++ ) {
	 int aaa[][] = aa[ mot ];
	 if ( aaa == null ) continue;
	 double distrib[] = ComputePolyScoreDistribution( mot );
	 DoubleUtils.UnitNorm( distrib );
	 out[ mot ] = DoubleUtils.New( aaa.length );
	 for ( int i = 0; i < aaa.length; i ++ ) {
	    double score = ComputePolyScore( aa, i, mot );
	    double sum = DoubleUtils.Sum( distrib, 0, (int) score );
	    out[ mot ][ i ] = 1.0 - sum; if ( out[ mot ][ i ] < 0.0 ) out[ mot ][ i ] = 0.0;
	 }
      }
      return out;
   }

   public static double ComputePolyScore( short seq[], int site, double cc[][] ) {
      int ww = cc.length;
      double score = 0.0;
      for ( int i = 0; i < ww; i ++ ) score += cc[ i ][ seq[ site + i ] ];
      return score;
   }

   public static double[] ComputePolyScores( short seq[], double cc[][], BackgroundModel bgmodel ) {
      double distrib[] = ComputePolyScoreDistribution( cc, bgmodel );
      DoubleUtils.UnitNorm( distrib );
      return ComputePolyScores( seq, cc, distrib, bgmodel, null );
   }

   public static double[] ComputePolyScores( short seq[], double cc[][], double distrib[], BackgroundModel bgmodel, double out[] ) {
      int WW = cc.length, max = seq.length - WW + 1;
      if ( out == null ) out = DoubleUtils.New( max );
      for ( int i = 0; i < max; i ++ ) {
	 double score = ComputePolyScore( seq, i, cc );
	 double sum = DoubleUtils.Sum( distrib, 0, (int) score );
	 out[ i ] = 1.0 - sum; if ( out[ i ] < 0.0 ) out[ i ] = 0.0;
      }
      return out;
   }

   public double ComputePolyScore( int aa[][][], int ii, int mot ) {
      int seqNum = aa[ mot ][ ii ][ 0 ];
      double c1[][] = c[ mot ];
      if ( seqNum < 0 ) return 0.0;
      short[] seq = sites.GetSequenceResidues( seqNum );
      int ww = aa[ mot ][ ii ].length - 1;
      double score = 0.0;
      for ( int i = 0; i < ww; i ++ ) {
	 int site = aa[ mot ][ ii ][ i + 1 ];
	 if ( site >= 0 && seq[ site ] >= 0 ) score += c1[ i ][ seq[ site ] ];
      }
      return score;
   }

   public double[][] ComputeLogoScores( int mot ) {
      return ComputeLogoScores( c[ mot ] );
   }

   public static double[][] ComputeLogoScores( double cc[][] ) {
      int W = cc.length, J = cc[ 0 ].length;
      double out[][] = DoubleUtils.New( cc ), defaultBg = 1.0 / (double) J;
      double sum[] = DoubleUtils.New( W ), q[][] = DoubleUtils.New( cc );
      double mmin = DoubleUtils.Min( q );
      if ( mmin < 0 ) DoubleUtils.Sub( q, mmin - 1e-3 );
      double maxSum = -Double.MAX_VALUE;
      for ( int i = 0; i < W; i ++ ) {
	 double ssum = DoubleUtils.Sum( q[ i ] );
	 if ( ssum > maxSum ) maxSum = ssum;
      }
      DoubleUtils.Divide( q, maxSum );
      for ( int i = 0; i < W; i ++ ) {
	 double[] q1i = q[ i ];
	 for ( short j = 0; j < J; j ++ ) {
	    if ( q1i[ j ] == 0 /*|| ( bgmodel.GetProb( j ) == 0 )*/ ) continue;
	    double qq = q1i[ j ]; // / bgmodel.GetProb( j ) / (double) J;
	    sum[ i ] += qq * DoubleUtils.Log2( qq );
	 }
	 sum[ i ] += DoubleUtils.Log2( (double) J );
	 for ( short j = 0; j < J; j ++ ) {
	    if ( q1i[ j ] == 0 /*|| ( bgmodel.GetProb( j ) == 0 )*/ ) continue;
	    double qq = q1i[ j ]; // / bgmodel.GetProb( j ) / (double) J;
	    out[ i ][ j ] = sum[ i ] * qq;
	 }
      }
      return out;
   }

   public double[] ComputeSPScores( int[][][] a ) {
      double out[] = DoubleUtils.New( nMotifs ); 
      ComputeSPScores( a, out );
      return out;
   }

   public void ComputeSPScores( int aa[][][], double SP[] ) {
      double probs[] = DoubleUtils.New( J );
      for ( int mot = 0; mot < nMotifs; mot ++ ) {
	 int aaa[][] = aa[ mot ];
	 if ( aaa == null ) break;
	 double cmot[][] = c[ mot ];
	 SP[ mot ] = 0.0;
	 for ( int i = 1; i < aaa.length; i ++ ) {
	    int seq1 = aaa[ i ][ 0 ];
	    if ( seq1 < 0 ) continue;
	    int aa1[] = aaa[ i ]; short ss1[] = sites.GetSequenceResidues( seq1 );
	    for ( int col = 1; col < aa1.length; col ++ ) {
	       double ci[] = cmot[ col-1 ];
	       double min = DoubleUtils.Min( ci );
	       if ( min < 0 ) DoubleUtils.Add( ci, -min );
	       fgmodel.GetProbs( ci, probs );
	       for ( int j = 0; j < i; j ++ ) {
		  int seq2 = aaa[ j ][ 0 ], aa2[] = aaa[ j ]; 
		  if ( seq2 < 0 ) continue;
		  short ss2[] = sites.GetSequenceResidues( seq2 );
		  if ( aa1[ col ] < 0 || aa2[ col ] < 0 ) continue;
		  if ( ss1[ aa1[ col ] ] < 0 || ss2[ aa2[ col ] ] < 0 ) continue;
		  if ( ss1[ aa1[ col ] ] == ss2[ aa2[ col ] ] ) SP[ mot ] += probs[ ss1[ aa1[ col ] ] ];
	       }
	    }
	 }
	 SP[ mot ] /= aa.length;
      }
   }

   public void ComputeLogOddsScores( int mot, double los[], 
				     int alignments[][][], double pvalues[], double alignThresh ) {
      for ( int ii = 0; ii < alignments[ mot ].length; ii ++ ) {
	 if ( pvalues != null ) {
	    double pp = pvalues != null ? pvalues[ ii ] : 
	       ComputeMatchScore( alignments, ii, mot );
	    if ( pp == 0 ) pp = 1.0 / Double.MAX_VALUE;
	    pp = DoubleUtils.Log10( pp );
	    if ( pp > -alignThresh ) { los[ ii ] = 0.0; continue; }
	 }
	 los[ ii ] = ComputeMatchScore( alignments, ii, mot );
      }
   }

   public double GetBestPValue( Sampler samp, short seq[], java.util.Map state ) {
      if ( state != null ) samp.SetState( state );
      int llen = seq.length;
      double out = 2.0;
      for ( int mot = 0; mot < nMotifs; mot ++ ) {
	 double distrib[] = ComputePolyScoreDistribution( mot );
	 DoubleUtils.UnitNorm( distrib );
	 double cc[][] = c[ mot ];
	 for ( int i = 0; i < llen - W; i ++ ) {
	    double score = 0.0;
	    for ( int j = 0; j < W; j ++ ) {
	       int site = i + j;
	       score += cc[ j ][ seq[ site ] ];
	    }
	    double sum = DoubleUtils.Sum( distrib, 0, (int) score );
	    double pvalue = 1.0 - sum; if ( pvalue < 0.0 ) pvalue = 0.0;
	    if ( pvalue < out ) out = pvalue;
	 }
      }
      return out;
   }
}
