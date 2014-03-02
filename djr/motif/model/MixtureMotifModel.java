package djr.motif.model;

import djr.motif.*;
import djr.motif.model.prior.*;
import djr.motif.model.*;
import djr.motif.sampler.sites.*;
import djr.util.array.*;
import djr.util.bio.Sequence;

/**
 * Describe class <code>MixtureMotifModel</code> here.
 *
 * @author <a href="mailto:dreiss@systemsbiology.org">David Reiss</a>
 * @version 1.9978 (Fri Nov 07 05:56:26 PST 2003)
 */
public class MixtureMotifModel extends AlignmentMotifModel {
   protected static String superModelID = "SUPERMODEL";

   double pseudo;
   boolean hasSuperModel;
   double cmix[][][] = null, c0[][] = null; // Supermodel
   int intPtrs[][], sh3Ptrs[][];

   public MixtureMotifModel( int nMotifs, Sites sites, int W,
			     int intPtrs[][], int sh3Ptrs[][], double pseud ) {
      super( nMotifs, sites, W );
      this.intPtrs = intPtrs;
      this.sh3Ptrs = sh3Ptrs;
      int NS = sites.GetSequenceCount();
      cmix = DoubleUtils.New( NS, W, J );
      q = DoubleUtils.New( NS, W, J );
      changed = new boolean[ NS ];
      this.pseudo = pseud;
      hasSuperModel = ( pseudo != 0 );
      if ( hasSuperModel ) c0 = DoubleUtils.New( W, J ); // supermodel
   }

   public MixtureMotifModel( MixtureMotifModel m ) {
      this( m.nMotifs, m.sites, m.W, m.intPtrs, m.sh3Ptrs, m.pseudo );
      CopyFrom( m );
   }

   public PSSMMotifModel Duplicate() {
      return new MixtureMotifModel( this );
   }

   public void Initialize() {
      super.Initialize();
      if ( hasSuperModel ) DoubleUtils.Zero( c0 );
      DoubleUtils.Zero( cmix );
   }

   public double[][] GetSuperModelCounts() { return c0; }
   public double GetSuperModelPseudoCount() { return pseudo; }

   public void Copy( PSSMMotifModel to ) {
      super.Copy( to );
      MixtureMotifModel mto = (MixtureMotifModel) to;
      mto.pseudo = pseudo;
      mto.hasSuperModel = hasSuperModel;
      if ( hasSuperModel ) DoubleUtils.Copy( mto.c0, c0 );
      DoubleUtils.Copy( mto.cmix, cmix );
      mto.intPtrs = intPtrs;
      mto.sh3Ptrs = sh3Ptrs;
      mto.W = W;
   }

   public void FillTheArrays() {
      Initialize();
      int NS = sites.GetSequenceCount();
      for ( int i = 0; i < NS; i ++ ) AddRemoveFromCountsForSequence( true, i );
   }

   public void AddRemoveFromCountsForSequence( boolean add, int seqNum ) {
      int sh3p[] = sh3Ptrs[ seqNum ];
      if ( sh3p == null ) return;
      int sh3len = sh3p.length;
      int ind = sites.GetSites( seqNum )[ 0 ];
      if ( ind < 0 ) return;
      short seq[] = sites.GetSequenceResidues( seqNum );
      int llen = seq.length;

      if ( hasSuperModel ) AddRemoveFromCounts( add, seq, ind, c0, pseudo );

      for ( int m = 0; m < sh3len; m ++ ) {
	 int mot = sh3p[ m ];
	 AddRemoveFromCounts( add, seq, ind, c[ mot ], 1.0 );

	 int intp[] = intPtrs[ mot ];
	 if ( intp == null ) continue;
	 int intlen = intp.length;
	 for ( int i = 0; i < intlen; i ++ ) {
	    int seqN = intp[ i ];
	    AddRemoveFromCounts( add, seq, ind, cmix[ seqN ], 1.0 );
	 }
      }
      changed[ seqNum ] = true;
   }

   public void FillQ() {
      int NS = sites.GetSequenceCount();
      for ( int i = 0; i < NS; i ++ ) FillQ( i );
   }

   public void FillQ( int mot ) {
      for ( int i = 0; i < W; i ++ ) FillQColumn( mot, i );      
   }

   public void FillQColumn( int seqNum, int col ) {
      double qq[] = q[ seqNum ][ col ];
      DoubleUtils.Copy( qq, cmix[ seqNum ][ col ] );
      // Add the supermodel into the mix
      if ( hasSuperModel ) DoubleUtils.Add( qq, c0[ col ] ); 
      fgmodel.GetProbs( qq, qq );
   }

   public double[][] GetMotifModel( int seqNum ) {
      if ( changed[ seqNum ] ) {
	 FillQ( seqNum );
	 changed[ seqNum ] = false;
      }
      return q[ seqNum ]; 
   }

   public void Print( java.io.PrintStream out, boolean ansi, boolean html, 
		      String names[], int mot ) {
      // Makes sure that the counts for each motif get printed out
      Print( out, ansi, html, names[ mot ], GetMotifCounts( mot ) ); }

   public int ComputeProbsForSequence( int seqNum, double[] qx, boolean normalize ) {
      return ComputeProbsForSequence( sites.GetSequenceResidues( seqNum ), GetMotifModel( seqNum ), qx, 
				      normalize );
   }

   public double ComputeMAPScore() {
      double F = 0.0; // Compute the total probability that all sequences' sites matches their mixture models.
      int NS = sites.GetSequenceCount();
      for ( int seqNum = 0; seqNum < NS; seqNum ++ ) {
	 int loc = sites.GetSites( seqNum )[ 0 ];
	 if ( loc < 0 ) continue;
	 double prob = ComputeProbForSite( sites.GetSequenceResidues( seqNum ), loc, 
					   GetMotifModel( seqNum ) );
	 if ( prob <= 0.0 ) continue; 
	 F += DoubleUtils.Log2( prob );
      }
      return F;
   }

   public boolean PrintSuperModel( java.io.PrintStream pout, boolean ansi, boolean html ) {
      if ( hasSuperModel ) {
	 double cc[][] = DoubleUtils.New( c0 );
	 DoubleUtils.Divide( cc, pseudo );
	 PSSMMotifModel.Print( pout, ansi, html, superModelID, cc, " %4d" );
	 pout.println();
      }
      return hasSuperModel;
   }
}
