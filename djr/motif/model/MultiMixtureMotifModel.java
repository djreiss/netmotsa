package djr.motif.model;

import djr.motif.*;
import djr.motif.model.prior.*;
import djr.motif.model.*;
import djr.motif.sampler.sites.*;
import djr.util.array.*;
import djr.util.bio.Sequence;

/**
 * Describe class <code>MultiMixtureMotifModel</code> here.
 *
 * @author <a href="mailto:reiss@uw.edu">David Reiss</a>
 * @version 1.9978 (Fri Nov 07 05:56:26 PST 2003)
 */
public class MultiMixtureMotifModel extends MixtureMotifModel {
   double pseudo2;
   boolean hasMixtureModels = false;
   transient double q2[][][][] = null;
   transient boolean changed2[][] = null;

   public MultiMixtureMotifModel( int nMotifs, Sites sites, int W,
				  int intPtrs[][], int sh3Ptrs[][], 
				  double pseud, double pseud2 ) {
      super( nMotifs, sites, W, intPtrs, sh3Ptrs, pseud );
      this.pseudo2 = pseud2;
      hasMixtureModels = ( pseudo2 != 0 );
      changed = null;
      q = null;
      InitializeQ2();
   }

   public MultiMixtureMotifModel( MultiMixtureMotifModel m ) {
      this( m.nMotifs, m.sites, m.W, m.intPtrs, m.sh3Ptrs, m.pseudo, m.pseudo2 );
      CopyFrom( m );
   }

   private void readObject( java.io.ObjectInputStream s ) throws java.io.IOException, 
								 ClassNotFoundException {
      s.defaultReadObject();
      InitializeQ2();
      BoolUtils.Set( changed2, true );
   }

   public PSSMMotifModel Duplicate() {
      return new MultiMixtureMotifModel( this );
   }

   public void Initialize() {
      super.Initialize();
      if ( q2 == null ) InitializeQ2();
      else {
	 DoubleUtils.Zero( q2 );
	 BoolUtils.Set( changed2, true );
      }
   }

   protected void InitializeQ2() {
      if ( q2 == null ) {
	 int NS = sites.GetSequenceCount();
	 q2 = new double[ NS ][ nMotifs ][][];
	 for ( int i = 0; i < NS; i ++ ) {
	    int ptrs[] = sh3Ptrs[ i ];
	    if ( ptrs == null ) continue;
	    for ( int j = 0; j < ptrs.length; j ++ )
	       q2[ i ][ ptrs[ j ] ] = new double[ W ][ J ];
	 }
	 changed2 = new boolean[ NS ][ nMotifs ];
      }
   }

   public void Copy( PSSMMotifModel to ) {
      super.Copy( to );
      MultiMixtureMotifModel mto = (MultiMixtureMotifModel) to;
      mto.pseudo2 = pseudo2;
      mto.hasMixtureModels = hasMixtureModels;
      DoubleUtils.Copy( mto.q2, q2 );
      BoolUtils.Copy( mto.changed2, changed2 );
      //System.err.println("HERE: "+ComputeMAPScore()+" "+mto.ComputeMAPScore());
   }

   public void FillTheArrays() {
      Initialize();
      for ( int mot = 0; mot < nMotifs; mot ++ ) {
	 int intp[] = intPtrs[ mot ];
	 if ( intp == null ) continue;
	 for ( int i = 0, intlen = intp.length; i < intlen; i ++ )
	    AddRemoveFromCounts( true, intp[ i ], mot );
      }
   }

   public void AddRemoveFromCounts( boolean add, int seqNum, int mot ) {
      if ( ! IntUtils.Contains( sh3Ptrs[ seqNum ], mot ) ) return;
      int ind = sites.GetSiteForMot( seqNum, mot );
      if ( ind < 0 ) return;

      short seq[] = sites.GetSequenceResidues( seqNum );

      AddRemoveFromCounts( add, seq, ind, c[ mot ], 1.0 );

      if ( hasMixtureModels ) 
	 AddRemoveFromCounts( add, seq, ind, cmix[ seqNum ], pseudo2 );

      if ( hasSuperModel ) AddRemoveFromCounts( add, seq, ind, c0, pseudo );

      changed2[ seqNum ][ mot ] = true;
   }

   public void FillQColumn( int seqNum, int mot, int col ) {
      double qq[] = q2[ seqNum ][ mot ][ col ];
      DoubleUtils.Copy( qq, c[ mot ][ col ] );
      // Add the mixture model into the mix
      if ( hasMixtureModels ) DoubleUtils.Add( qq, cmix[ seqNum ][ col ] );
      // Add the supermodel into the mix
      if ( hasSuperModel ) DoubleUtils.Add( qq, c0[ col ] ); 
      fgmodel.GetProbs( qq, qq );
   }

   public void FillQ() {
      int NS = sites.GetSequenceCount();
      for ( int mot = 0; mot < nMotifs; mot ++ ) 
	 for ( int i = 0; i < NS; i ++ ) FillQForSequenceAndMotif( i, mot );
   }

   public double[][] GetMotifModel( int seqNum ) { return null; }
   public void FillQ( int mot ) { };

   public void FillQForSequenceAndMotif( int seqNum, int mot ) {
      for ( int i = 0; i < W; i ++ ) FillQColumn( seqNum, mot, i ); 
   }

   public int ComputeProbsForSequence( int seqNum, int mot, double[] qx, 
				       boolean normalize ) {
      return ComputeProbsForSequence( sites.GetSequenceResidues( seqNum ), 
				      GetMotifModel( seqNum, mot ), qx, normalize );
   }

   public double[][] GetMotifModel( int seqNum, int mot ) {
      if ( changed2[ seqNum ][ mot ] ) {
	 FillQForSequenceAndMotif( seqNum, mot );
	 changed2[ seqNum ][ mot ] = false;
      }
      return q2[ seqNum ][ mot ]; 
   }

   public double ComputeMAPScore() {
      double F = 0.0;
      int NS = sites.GetSequenceCount();
      for ( int seqNum = 0; seqNum < NS; seqNum ++ ) {
	 short seq[] = sites.GetSequenceResidues( seqNum );
	 int sh3p[] = sh3Ptrs[ seqNum ];
	 if ( sh3p == null ) continue;
	 for ( int j = 0, nsh3p = sh3p.length; j < nsh3p; j ++ ) {
	    int mot = sh3p[ j ];
	    int loc = sites.GetSiteForMot( seqNum, mot );
	    if ( loc < 0 ) continue;
	    double prob = ComputeProbForSite( seq, loc, GetMotifModel( seqNum, mot ) );
	    if ( prob <= 0.0 ) continue; 
	    F += DoubleUtils.Log2( prob );
	 }
      }
      return F;
   }
}
