package djr.motif.model;
import java.util.Map;
import java.io.PrintStream;

import cern.jet.random.*;

import djr.motif.model.prior.*;
import djr.util.*;
import djr.util.array.*;
import djr.util.bio.Sequence;

/**
 * Describe class <code>PSSMMotifModel</code> here.
 *
 * @author <a href="mailto:reiss@uw.edu">David Reiss</a>
 * @version 1.9978 (Fri Nov 07 05:56:26 PST 2003)
 */
public class PSSMMotifModel implements java.io.Serializable {
   public static String FORMAT = " %5.3f";

   protected int W, nMotifs, origW;
   protected int J;
   protected double Wdistrib[];
   transient protected double q[][][];
   transient protected boolean changed[];
   //boolean columns[][];
   protected BackgroundModel bgmodel;

   public PSSMMotifModel() { };

   public PSSMMotifModel( int nMotifs, int J, int W ) {
      this.nMotifs = nMotifs;
      this.W = W;
      this.J = J;
      this.origW = W;
      InitializeQ();
   }

   public PSSMMotifModel( PSSMMotifModel m ) {
      this( m.nMotifs, m.J, m.W );
      CopyFrom( m );
   }

   private void readObject( java.io.ObjectInputStream s ) throws java.io.IOException, 
								 ClassNotFoundException {
      s.defaultReadObject();
      InitializeQ();
      BoolUtils.Set( changed, true );
   }

   protected void InitializeQ() {
      q = new double[ nMotifs ][][];
      for ( int i = 0; i < nMotifs; i ++ ) q[ i ] = DoubleUtils.New( W, J );
      changed = new boolean[ nMotifs ];

      /*columns = BoolUtils.New( nMotifs, Wmax );
	for ( int i = 0; i < nMotifs; i ++ ) {
	int start = Wmax / 2 - W[ i ] / 2;
	for ( int j = start, jj = 0; jj < W[ i ]; j ++, jj ++ ) columns[ i ][ j ] = true;
	}*/
   }

   public PSSMMotifModel Duplicate() {
      return new PSSMMotifModel( this );
   }

   public BackgroundModel GetBackgroundModel() { return bgmodel; }
   public int GetMotifCount() { return nMotifs; }
   public int GetMotifWidth( int mot ) { return W; }
   public int GetAlphabetSize() { return J; }

   public double[][] GetMotifModel( int mot ) {
      if ( changed[ mot ] ) {
	 FillQ( mot );
	 changed[ mot ] = false;
      }
      return q[ mot ]; 
   }

   public void SetChanged( int mot, boolean ch ) { changed[ mot ] = ch; }

   public void SetBackgroundAndForegroundModels( String bgfile, int bgorder, double pseudo,
						 String fgType, String dirichName, 
						 String matrix ) {
      if ( ! MyUtils.IsNullString( bgfile ) )
	 bgmodel = new BackgroundModel( bgfile, (short) bgorder, pseudo );
      else bgmodel = null;
   }

   public void SetMotifModels( double qq[][][] ) {
      this.q = qq;
      changed = new boolean[ qq.length ];
      nMotifs = qq.length;
      W = qq[ 0 ].length;
   }

   public void Initialize() {
      if ( q != null ) DoubleUtils.Zero( q );
      if ( changed != null ) BoolUtils.Set( changed, true );
   }

   public void Copy( PSSMMotifModel to ) {
      if ( to.q != null ) DoubleUtils.Copy( to.q, q );
      else if ( q != null ) to.q = DoubleUtils.New( q );
      if ( to.changed != null ) BoolUtils.Copy( to.changed, changed );
      else if ( changed != null ) to.changed = BoolUtils.New( changed );
      to.W = W;
      to.bgmodel = bgmodel;
      to.nMotifs = nMotifs;
   }
   
   public void CopyFrom( PSSMMotifModel from ) {
      from.Copy( this );
   }

   public void Print( PrintStream out, boolean ansi, boolean html, String names[] ) {
      for ( int mot = 0; mot < nMotifs; mot ++ ) {
	 Print( out, ansi, html, names, mot );
	 if ( mot < nMotifs && html ) 
	    out.print( HTML.cellEnd() + HTML.NEW_ROW + HTML.cellStartA( "left" ) + 
		       HTML.SPACE + HTML.cellEnd() + HTML.NEW_ROW );
	 if ( ! html ) out.println();
      }
      if ( ! html ) out.println(); else out.print( HTML.reset() );
   }

   public void Print( PrintStream out, boolean ansi, boolean html, String names[], int mot ) {
      Print( out, ansi, html, names != null ? names[ mot ] : "" + mot, GetMotifModel( mot ) ); }

   public static void Print( PrintStream out, boolean ansi, boolean html, String name, 
			     double cc[][], String format ) {
      String save = FORMAT;
      FORMAT = format;
      Print( out, ansi, html, name, cc );
      FORMAT = save;
   }

   public static void Print( PrintStream out, boolean ansi, boolean html, String name, 
			     double cc[][] ) {
      boolean intFormat = " %4d".equals( FORMAT );
      int J = cc[ 0 ].length;
      short JJ = (short) J;
      int WW = cc.length;
      double[] hi = DoubleUtils.New( WW );
      if ( ansi || html ) for ( int i = 0; i < WW; i ++ ) { for ( int j = 0; j < J; j ++ )
	 if ( cc[ i ][ j ] > hi[ i ] ) hi[ i ] = cc[ i ][ j ]; }
      for ( int j = 0; j < J; j ++ ) {
	 char ntide = 0; 
	 if ( J == 4 ) ntide = (char) Sequence.SEQ_TABLE[ j ];
	 else if ( J != 4 ) ntide = (char) Sequence.PROT_TABLE[ j ];
	 if ( html ) out.print( HTML.cellStart() );
	 out.print( name + " " + ntide + " " );
	 if ( html ) out.print( HTML.SPACE + HTML.SPACE + HTML.cellEnd() );
	 for ( int i = 0; i < WW; i ++ ) {
	    if ( ansi ) {
	       if ( cc[ i ][ j ] / hi[ i ] >= 0.9 ) 
		  out.print( ANSI.setAttr( Sequence.getANSIColorBG( ntide, JJ ) + ";" + ANSI.BOLD ) 
			     );
	       else if ( cc[ i ][ j ] / hi[ i ] >= 0.75 ) 
		  out.print( ANSI.setAttr( Sequence.getANSIColorBG( ntide, JJ ) ) );
	       else if ( cc[ i ][ j ] / hi[ i ] >= 0.5 ) 
		  out.print( ANSI.setAttr( Sequence.getANSIColorFG( ntide, JJ ) + ";" + ANSI.BOLD ) 
			     );
	       else if ( cc[ i ][ j ] / hi[ i ] >= 0.25 ) 
		  out.print( ANSI.setAttr( Sequence.getANSIColorFG( ntide, JJ ) ) );
	    } else if ( html ) {
	       String colorfg = Sequence.getHTMLColorFG( (short) j, JJ );
	       String colorbg = Sequence.getHTMLColorBG( (short) j, JJ );
	       if ( cc[ i ][ j ] / hi[ i ] >= 0.9 ) 
		  out.print( HTML.cellStart( colorbg ) + HTML.BOLD );
	       else if ( cc[ i ][ j ] / hi[ i ] >= 0.75 ) 
		  out.print( HTML.cellStart( colorbg ) );
	       else if ( cc[ i ][ j ] / hi[ i ] >= 0.5 ) 
		  out.print( HTML.cellStart() + HTML.colorStart( colorfg ) + HTML.BOLD );
	       else if ( cc[ i ][ j ] / hi[ i ] >= 0.25 ) 
		  out.print( HTML.cellStart() + HTML.colorStart( colorfg ) );
	       else out.print( "<td>" );
	    }
	    if ( ! intFormat ) out.print( DoubleUtils.SPrintf( FORMAT, cc[ i ][ j ] ) );
	    else out.print( IntUtils.SPrintf( FORMAT, (int) ( cc[ i ][ j ] + 0.5 ) ) );
	    if ( ansi && cc[ i ][ j ] / hi[ i ] >= 0.25 )
	       out.print( ANSI.setAttr( ANSI.RESET ) );
	    else if ( html ) out.print( HTML.cellEnd() );
	 }
	 if ( html ) out.print( HTML.NEW_ROW ); else out.println(); 
      }
   }

   public double ComputeProbForSite( short seq[], int site, int mot ) {
      return ComputeProbForSite( seq, site, q[ mot ] );
   }

   public double ComputeProbForSite( short seq[], int site, double qq[][] ) {
      double qqq = 1.0;
      int WW = W; if ( site + WW > seq.length ) WW = seq.length - site;
      for ( int j = 0; j < WW; j ++, site ++ ) {
	 if ( seq[ site ] >= 0 ) // not masked
	    qqq *= qq[ j ][ seq[ site ] ] / bgmodel.GetProb( seq, site );
	 else { qqq = 0.0; break; } // masked => set prob to zero and break
      }
      return qqq;
   }

   public int ComputeProbsForSequence( short seq[], double qq[][], double[] qx, 
				       boolean normalize ) {
      int whereMax = -1, maxl = seq.length - qq.length + 1; 
      double max = -Double.MAX_VALUE;
      for ( int i = 0; i < maxl; i ++ ) {
	 double qqq = ComputeProbForSite( seq, i, qq );
	 if ( qqq > max ) { whereMax = i; max = qqq; }
	 qx[ i ] = qqq;
      }
      if ( normalize ) DoubleUtils.Divide( qx, max );
      return whereMax;
   }

   public int ComputeProbsForSequence( short seq[], int mot, double[] qx, 
				       boolean normalize ) {
      return ComputeProbsForSequence( seq, q[ mot ], qx, normalize );
   }
      
   public void FillQColumn( int mot, int col ) {
      DoubleUtils.UnitNorm( q[ mot ][ col ] );
   }

   public void FillQ( int mot ) {
      for ( int i = 0; i < W; i ++ ) FillQColumn( mot, i );
   }

   public void FillQ() {
      for ( int i = 0; i < nMotifs; i ++ ) FillQ( i );
   }

   public double ComputeSimilarityScore( int m1, int m2 ) {
      return ComputeSimilarityScore( GetMotifModel( m1 ), GetMotifModel( m2 ) );
   }

   public static double ComputeSimilarityScore( double q1[][], double q2[][] ) {
      return ComputeSimilarityScore( q1, q2, true );
   }

   /**
    * Compute similarity of two PSSM frequency matrices using the Pearson
    *	correlation coeff. (Pietrokovski (1996). Ranges from -1 (no match) to +1.
    */
   public static double ComputeSimilarityScore( double q1[][], double q2[][], 
						boolean shiftIt ) {
      int nsum = 0, WW = q1.length, J = q1[ 0 ].length;
      double scores[] = DoubleUtils.New( WW + 3 ), max = -Double.MAX_VALUE;
      for ( int shift = -WW/2-1; shift <= WW/2+1; shift ++ ) {
	 if ( ! shiftIt && shift != 0 ) continue;
	 int ncols = 0;
	 for ( int i = 0; i < WW; i ++ ) {
	    int ishift = i + shift;
	    if ( i >= q2.length || ishift < 0 || ishift >= q1.length ) continue;
	    ncols ++;
	    double sum1 = 0.0, sum2 = 0.0, sum3 = 0.0;
	    double mean1 = DoubleUtils.Mean( q1[ ishift ] );
	    double mean2 = DoubleUtils.Mean( q2[ i ] );
	    for ( int j = 0; j < J; j ++ ) {
	       double a = q1[ ishift ][ j ] - mean1;
	       double b = q2[ i ][ j ] - mean2;
	       sum1 += a * b;
	       sum2 += a * a;
	       sum3 += b * b;
	    }
	    scores[ nsum ] += sum1 / Math.sqrt( sum2 * sum3 );
	 }
	 scores[ nsum ] /= ncols;
	 if ( scores[ nsum ] > max ) max = scores[ nsum ];
	 nsum ++;
      }
      return max;
   }

   public boolean IsSimilarTo( PSSMMotifModel other, int mot, double thresh ) {
      return ComputeSimilarityScore( GetMotifModel( mot ), other.GetMotifModel( mot ) ) > thresh;
   }

   public boolean IsValid( int mot ) {
      return DoubleUtils.Sum( q[ mot ] ) > 0;
   }

   // Generate a sequence sampled from one of the PSSMs
   public short[] GenerateSequence( int mot ) {
      short res[] = ShortUtils.New( W );
      double temp[] = DoubleUtils.New( J );
      for ( int i = 0; i < W; i ++ ) {
	 DoubleUtils.Copy( temp, q[ mot ][ i ] );
	 DoubleUtils.MaxNorm( temp );
	 res[ i ] = (short) DoubleUtils.Sample( temp );
      }
      return res;
   }

   // Generate a sequence sampled from a mixture of two of the PSSMs
   public short[] GenerateSequence( int mot1, int mot2 ) {
      short res[] = ShortUtils.New( W );
      double temp[] = DoubleUtils.New( J );
      for ( int i = 0; i < W; i ++ ) {
	 DoubleUtils.Copy( temp, q[ mot1 ][ i ] );
	 DoubleUtils.Add( temp, q[ mot2 ][ i ] );
	 DoubleUtils.MaxNorm( temp );
	 res[ i ] = (short) DoubleUtils.Sample( temp );
      }
      return res;
   }

   // Generate a sequence sampled from a mixture of any number of the PSSMs
   public short[] GenerateSequence( int mots[] ) {
      if ( mots.length <= 1 ) return GenerateSequence( mots[ 0 ] );
      short res[] = ShortUtils.New( W );
      double temp[] = DoubleUtils.New( J );
      for ( int i = 0; i < W; i ++ ) {
	 DoubleUtils.Copy( temp, q[ mots[ 0 ] ][ i ] );
	 for ( int j = 1; j < mots.length; j ++ ) 
	    DoubleUtils.Add( temp, q[ mots[ j ] ][ i ] );
	 DoubleUtils.MaxNorm( temp );
	 res[ i ] = (short) DoubleUtils.Sample( temp );
      }
      return res;
   }

   public static PSSMMotifModel GenerateFromConsensi( String consensi ) {
      return GenerateFromConsensi( consensi, 1000 );
   }

   // Generate a PSSM model from a set of consensi (separated by commas or spaces)
   public static PSSMMotifModel GenerateFromConsensi( String consensi, double cnt ) {
      String toks[] = MyUtils.Tokenize( consensi, ", " );
      boolean prot = false;
      int len[] = IntUtils.New( toks.length );
      for ( int i = 0; i < toks.length; i ++ ) {
	 len[ i ] = toks[ i ].length();
	 int sum = 0;
	 if ( toks[ i ].indexOf( '[' ) >= 0 ) {
	    int ind = toks[ i ].indexOf( '[' ), ind2 = toks[ i ].indexOf( ']', ind );
	    sum += ( ind2 - ind );
	    while ( toks[ i ].indexOf( '[', ind2 ) >= 0 ) {
	       ind = toks[ i ].indexOf( '[', ind2 );
	       ind2 = toks[ i ].indexOf( ']', ind );
	       sum += ( ind2 - ind );
	    }
	 }
	 len[ i ] -= sum;
	 short res[] = ShortUtils.New( len[ i ] );
	 Sequence s = new Sequence( res );
	 for ( int x = 0; x < 100; x ++ ) {
	    s.InsertMotifs( toks[ i ] );
	    if ( s.GetType() == Sequence.PROTEIN ) prot = true;
	    //System.out.println(s);
	 }
      }
      int J = prot ? 20 : 4;
      PSSMMotifModel out = new PSSMMotifModel( toks.length, J, IntUtils.Max( len ) );

      double c[][][] = DoubleUtils.New( out.q );
      for ( int i = 0; i < toks.length; i ++ ) {
	 DoubleUtils.Zero( c[ i ] );
	 short res[] = ShortUtils.New( len[ i ] );
	 Sequence s = new Sequence( res );
	 for ( int x = 0; x < 1000; x ++ ) {
	    s.InsertMotifs( toks[ i ] );
	    res = s.GetResidues();
	    for ( int j = 0; j < res.length; j ++ ) c[ i ][ j ][ res[ j ] ] ++;
	 }
      }

      BackgroundModel bgmodel = new BackgroundModel( "sequences/yeast.fst.gz", 
						     (short) 0, 0.0001 );
      MotifModelPrior fgmodel = null;
      if ( J == 20 ) fgmodel = MotifModelPrior.GetPrior( (short) J, "dirich", 
				       "dirichlet/hydro-cons.3comp", "", 0.0001, bgmodel );
      else fgmodel = MotifModelPrior.GetPrior( (short) J, "pseudo", "", "", 0.0001, bgmodel );

      for ( int i = 0; i < toks.length; i ++ )
	 for ( int j = 0; j < c[ i ].length; j ++ ) {
	    DoubleUtils.Mult( c[ i ][ j ], cnt / 1000.0 );
	    fgmodel.GetProbs( c[ i ][ j ], out.q[ i ][ j ] );
	 }

      return out;
   }

   // Generate a set of background PSSM models from a set of sequences
   public static PSSMMotifModel GenerateFromSequences( Sequence S[], int cnt, int W ) {
      PSSMMotifModel out = new PSSMMotifModel( cnt, S[ 0 ].GetAlphabetSize(), W );
      double c[][][] = out.q;

      int slen = S.length;
      boolean visited[] = new boolean[ slen ];
      int numPerMotif = slen / cnt;
      for ( int mot = 0; mot < cnt; mot ++ ) {
	 double counts[][] = c[ mot ];
	 for ( int i = 0; i < numPerMotif; i ++ ) {
	    int ind = -1;
	    if ( mot < cnt - 1 ) ind = IntUtils.RandChoose( 0, slen - 1 );
	    else ind = BoolUtils.IndexOf( visited, false );
	    if ( ind < 0 || visited[ ind ] ) { i --; continue; }
	    short seq[] = S[ ind ].GetResidues();
	    int llen = seq.length, max = llen - W + 1;
	    for ( int loc = 0; loc < max; loc ++ )
	       for ( int j = 0; j < W; j ++ )
		  counts[ j ][ (int) seq[ loc + j ] ] ++;
	    visited[ ind ] = true;
	    i ++;
	 }
	 out.FillQ( mot );
      }
      return out;
   }

   /*public static void main( String args[] ) {
      Sequence S[] = null;
      try { S = Sequence.ReadSequences( args[ 0 ] ); }
      catch( Exception e ) { e.printStackTrace(); System.exit( 0 ); }
      
      PSSMMotifModel out = GenerateFromSequences( S, Integer.parseInt( args[ 1 ] ), 9 );
      out.Print( System.out, true, false, null );
      }*/

   /*public static void main( String args[] ) {
      PSSMMotifModel out = GenerateFromConsensi( args[ 0 ], 10 );

      String mots[] = MyUtils.Tokenize( args[ 0 ], ", " );
      out.Print( System.out, true, false, mots );

      for ( int i = 0; i < mots.length; i ++ ) {
	 for ( int j = 0; j < 10; j ++ ) 
	    System.out.println( new Sequence( out.GenerateSequence( i ) ) );
	 System.out.println();
      }
   }
   */
}
