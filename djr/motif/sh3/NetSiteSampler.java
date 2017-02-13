package djr.motif.sh3;

import java.io.*;
import java.util.*;
import cern.jet.random.*;

import djr.nr.*;
import djr.util.*;
import djr.util.bio.*;
import djr.util.array.*;
import djr.motif.sampler.*;
import djr.motif.model.prior.*;
import djr.motif.sampler.sites.*;
import djr.motif.model.*;

/**
 * Class <code>NetSiteSampler</code>
 *
 * @author <a href="mailto:reiss@uw.edu">David Reiss</a>
 * @version 1.9978 (Fri Nov 07 05:56:26 PST 2003)
 */
public class NetSiteSampler extends SiteSampler {
   protected Map sh3s, ints;
   protected transient Map seqsHash, culledIntsHash;
   protected ObjVector sh3List;
   protected int nSh3, nInts, intPtrs[][], nonIntPtrs[][], sh3Ptrs[][], nonSh3Ptrs[][];

   protected double nonDiscrimF = 0, discrimF = 0, discrimPseudo = 0.5; /* <== pseudofrequency */
   protected double discrimLimit = 0.1;
   protected boolean implantPxxP = false;
   protected transient double ttestData[][][];
   protected transient double qtemp[][][];
   protected transient TTest tTest = null;

   public NetSiteSampler() { };
   public NetSiteSampler( String argv[] ) { Initialize( argv ); }
   public NetSiteSampler( String args ) { this( MyUtils.Tokenize( args, " " ) ); }

   protected boolean Initialize( String argv[] ) {
      argProc = new ArgProcessor();
      if ( ! argProc.ProcessArgs( this, argv ) ) MyUtils.Exit( 0 );
      sh3List = getSh3List( argProc.getArg( "sh3" ) );

      String yeastFile = "sequences/yeast_culled_nosh3.fst.gz";
      Sequence allYeast[];
      try { allYeast = Sequence.ReadSequences( argProc.getArg( "bgf", yeastFile ) ); }
      catch( Exception e ) { allYeast = null; System.err.println( yeastFile + ":" );
      e.printStackTrace(); }

      // Remove all sequences in yeast file that are similar to our SH3 interactor seqs.
      String non = argProc.getArg( "non", null );
      if ( ! MyUtils.IsNullString( non ) && ! yeastFile.equals( non ) ) {
	 ObjVector newYeast = new ObjVector();
	 for ( int j = 0, ss = allYeast.length; j < ss; j ++ ) {
	    allYeast[ j ].SetHeader( Sequence.GetGeneNameFromYeastGenome( 
									 allYeast[ j ].GetHeader() ).toUpperCase() );
	    boolean similar = false;
	    for ( int i = 0, s = S.length; i < s; i ++ )
	       if ( seqsAreEqual( S[ i ], allYeast[ j ] ) ) { similar = true; break; }
	    if ( ! similar ) newYeast.addElement( allYeast[ j ] );
	 }
	 allYeast = new Sequence[ newYeast.size() ];
	 for ( int i = 0, s = newYeast.size(); i < s; i ++ ) 
	    allYeast[ i ] = (Sequence) newYeast.elementAt( i );
	 newYeast = null;
	 Sequence.WriteSequences( yeastFile, allYeast );
	 argProc.putArg( "bgf", yeastFile );
      }

      initialized = super.Initialize( argv );
      if ( ! initialized ) return false;
      initialized = false;

      implantPxxP = argProc.getBooleanArg( "pxxp" );

      int i = 0;
      if ( ! ranking ) {
	 //this.alignThresh = 0;

	 seqsHash = Sequence.SequenceArrayToHashtable( S );
	 Map sif = readSif( argProc.getArg( "sif" ), true );
	 sh3s = getSh3Interactors( sif );
	 ints = getInteractorSh3s( sh3s );
	 nSh3 = sh3s.size();
	 nInts = ints.size();

	 modelIDs = new String[ nSh3 ];
	 Iterator e = sh3s.keySet().iterator(); 
	 while( e.hasNext() ) modelIDs[ i ++ ] = (String) e.next();

	 fillIntPtrs();
	 fillSh3Ptrs();

	 if ( verbose ) {
	    println( "A total of " + nSh3 + " SH3 proteins; " + nInts + 
		     " interactors; " + IntUtils.NElem( intPtrs ) + " interactions." );
	    double mean = 0; 
	    for ( i = 0; i < nInts; i ++ ) 
	       if ( sh3Ptrs[ i ] != null ) mean += sh3Ptrs[ i ].length;
	    printf( "An average of %.2f SH3 partners per interactor; ", 
		    mean / (double) nInts );

	    for ( i = 0, mean = 0; i < nSh3; i ++ ) 
	       if ( intPtrs[ i ] != null ) mean += intPtrs[ i ].length;
	    printf( "%.2f interactors per SH3 protein.\n", mean / (double) nSh3 );
	 }
      }

      Iterator e = sh3s.keySet().iterator(); i = 0;
      while( e.hasNext() ) modelIDs[ i ++ ] = (String) e.next();

      this.nMotifs = nSh3;      

      double pseudo2 = argProc.getFloatArg( "P2" );

      Sites sites = new Sites( S, W, nMotifs, minPerSeq, maxPerSeq );

      mmodel = new MixtureMotifModel( nMotifs, sites, W, intPtrs, sh3Ptrs, pseudo2 );
      mmodel.SetBackgroundAndForegroundModels( bgfile, bgorder, pseudo, fgType, 
					       dirichName, matrix );

      discrimPseudo = argProc.getFloatArg( "PD" );

      initialized = true;
      return initialized;
   }

   protected boolean seqsAreEqual( Sequence s1, Sequence s2 ) {
      return ( ( s1.GetHeader().equals( s2.GetHeader() ) ) ||
	       s1.IsHomologous( s2 ) );
   }

   protected boolean doDiscrimination() {
      return ( discrimPseudo > 0.0 && ( /*iter1 > maxIter1*9/10 ||*/ iter2 > maxIter2*3/4 || nworse > 0 ) );
   }

   protected void FillTheArrays() { 
      super.FillTheArrays(); 
      MixtureMotifModel mmod = (MixtureMotifModel) mmodel;
      if ( implantPxxP ) {
	 for ( int q = 0; q < nMotifs; q ++ ) 
	    implantPxxP( mmod.GetMotifCounts( q ) );
	 implantPxxP( mmod.GetSuperModelCounts() );
      }
   }

   protected void IterateSampler( int niter ) {
      int itermax = niter == -1 ? NS : niter;
      MixtureMotifModel mmod = (MixtureMotifModel) mmodel;
      for ( int i = 0; i < itermax; i ++ ) {
	 int seqNum = niter == -1 ? i : IntUtils.RandChoose( 0, NS-1 );
	 mmod.AddRemoveFromCountsForSequence( false, seqNum );

	 int max = ComputeProbsForSequence( seqNum, qx, true );
	 if ( ! noSamp ) {
	    int mmax = DoubleUtils.Sample( qx, len[ seqNum ] - W + 1 );
	    if ( mmax != -1 ) max = mmax;
	 }
	 if ( max != -1 /*&& IsSiteValid( max, seqNum )*/ ) {
	    GetSites().RemoveAllSites( seqNum );
	    GetSites().AddSite( seqNum, max, 0, true );
	 }
	 mmod.AddRemoveFromCountsForSequence( true, seqNum );
      }
   }   

   protected int ComputeProbsForSequence( int seqNum, double qx[], boolean normalize ) {
      MixtureMotifModel mmod = (MixtureMotifModel) mmodel;
      int max = mmod.ComputeProbsForSequence( seqNum, qx, normalize );
      if ( power != 1.0 ) DoubleUtils.Pow( qx, power );
      if ( doDiscrimination() ) max = AdjustProbsUsingDiscrimination( qx, seqNum, normalize, max );
      return max;
   }

   protected void initializeStuff() {
      if ( doDiscrimination() && ttestData == null ) {
	 ttestData = new double[ NS ][ 2 ][];
	 for ( int i = 0; i < NS; i ++ ) {
	    ttestData[ i ][ 0 ] = DoubleUtils.New( sh3Ptrs[ i ].length );
	    ttestData[ i ][ 1 ] = DoubleUtils.New( nonSh3Ptrs[ i ].length );
	 }
	 qtemp = DoubleUtils.New( nMotifs, W, J );
	 tTest = new TTest();
      }
   }

   public int AdjustProbsUsingDiscrimination( double qx[], int seqNum, boolean normalize,
					      int max ) {
      initializeStuff();

      // Compute motif models for the individual motif counts arrays 
      // (not using the mixture arrays that the mmodel would give us 
      // if we called mmodel.GetMotifModel(mot)).
      for ( int mot = 0; mot < nSh3; mot ++ )
	 ( (AlignmentMotifModel) mmodel ).
	    FillQ( mmodel.GetMotifCounts( mot ), qtemp[ mot ] );

      double ttest[][] = ttestData[ seqNum ];
      int sh3p[] = sh3Ptrs[ seqNum ], nonSh3p[] = nonSh3Ptrs[ seqNum ];
      int llen = len[ seqNum ] - W + 1;
      double oldMax = qx[ max ], limit = oldMax * discrimLimit;
      double newMax = -Double.MAX_VALUE; int newMaxSite = -1;
      short seq[] = GetSites().GetSequenceResidues( seqNum );
      for ( int site = 0; site < llen; site ++ ) {
	 if ( qx[ site ] >= limit ) { // Only do it for high-scoring sites (to speed it up)
	    double score = ComputeDiscriminationScore( seq, site, sh3p, nonSh3p, ttest );
	    //if ( site == max ) System.err.println("HERE: "+seqNum+" "+score+" "+sh3p.length);
	    qx[ site ] *= score;
	 }
	 if ( normalize && qx[ site ] > newMax ) { 
	    newMax = qx[ site ]; 
	    newMaxSite = site;
	 }
      }
      if ( normalize )
	 for ( int site = 0; site < llen; site ++ ) qx[ site ] /= newMax;
      return newMaxSite;
   }

   protected double ComputeDiscriminationScore( short seq[], int site, 
						int sh3p[], int nonSh3p[], 
						double ttest[][] ) {
      double score = 1.0;
      for ( int j = 0, sz = ttest[ 0 ].length; j < sz; j ++ ) 
	 ttest[ 0 ][ j ] = DoubleUtils.Log( mmodel.ComputeProbForSite( seq, site, qtemp[ sh3p[ j ] ] ) );
      for ( int j = 0, sz = ttest[ 1 ].length; j < sz; j ++ ) 
	 ttest[ 1 ][ j ] = DoubleUtils.Log( mmodel.ComputeProbForSite( seq, site, qtemp[ nonSh3p[ j ] ] ) );

      try {
	 tTest.recompute( ttest[ 1 ], ttest[ 0 ] );
	 score = ( tTest.getTValue() > 0 ? ( 1.0 - tTest.getProb() ) : 0 ) + discrimPseudo; 
      } catch( Exception e ) { score = 1.0; e.printStackTrace(); }
      return score;
   }

   public double ComputeMAPScore( Map state ) {
      if ( state == null ) state = selfs;
      double F = super.ComputeMAPScore( state );
      nonDiscrimF = F;

      // Only do discrimination if iters are high (doing it all iters is too slow)
      if ( doDiscrimination() ) {
	 initializeStuff();
	 discrimF = 0;

	 AlignmentMotifModel mm = (AlignmentMotifModel) state.get( "mm" );
	 Sites sites = mm.GetSites();

	 for ( int mot = 0; mot < nSh3; mot ++ )
	    mm.FillQ( mm.GetMotifCounts( mot ), qtemp[ mot ] );

	 for ( int mot = 0; mot < nSh3; mot ++ ) {
	    int intp[] = intPtrs[ mot ], intpl = intp.length;
	    if ( intpl <= 0 ) continue;

	    for ( int i = 0; i < intpl; i ++ ) {
	       int seqNum = intp[ i ];
	       double ttest[][] = ttestData[ seqNum ];
	       int sh3p[] = sh3Ptrs[ seqNum ], nonSh3p[] = nonSh3Ptrs[ seqNum ];
	       short seq[] = sites.GetSequenceResidues( seqNum );
	       ObjVector sitesList = new ObjVector();
	       sites.GetAlignmentsForSequence( seqNum, mot, sitesList );
	       if ( sitesList.size() > 0 ) {
		  int site = ( (int[]) sitesList.elementAt( 0 ) )[ 1 ];
		  double score = ComputeDiscriminationScore( seq, site, sh3p, nonSh3p, ttest );
		  F += DoubleUtils.Log2( score );
		  discrimF += DoubleUtils.Log2( score );
	       }
	    }
	 }
      }
      return F;
   }

   public int[][][] GetAlignments() {
      int out[][][] = new int[ nMotifs ][][];
      for ( int mot = 0; mot < nMotifs; mot ++ ) {
	 ObjVector tempout = new ObjVector( 10, 10 );
	 for ( int i = 0; i < NS; i ++ ) {
	    int partners[] = sh3Ptrs[ i ];
	    if ( partners == null || ! IntUtils.Contains( partners, mot ) ) continue;
	    GetSites().GetAlignmentsForSequence( i, tempout );
	 }
	 out[ mot ] = new int[ tempout.size() ][];
	 for ( int i = 0; i < tempout.size(); i ++ ) 
	    out[ mot ][ i ] = (int[]) tempout.elementAt( i );
      }
      return out;
   }

   public void RestoreStateFromFile( String saveStateFile ) {
      super.RestoreStateFromFile( saveStateFile );
      nworse = 1; // Make it so it displays the discrimination MAP score.
   }

   protected void PrintSummary( Map state, boolean printAlignments,
				boolean printGarbage ) {
      SamplerOutput.PrintSummary( this, state, printAlignments, printGarbage ); 
      MixtureMotifModel mm = (MixtureMotifModel) state.get( "mm" );
      boolean hasSuper = mm.PrintSuperModel( pout, ansi, html );
      if ( hasSuper && vverbose ) {
	 if ( "text".equals( logos ) ) SamplerOutput.PrintLogos( this );
	 else if ( "gui".equals( logos ) ) DisplayLogos();
      }
      if ( doDiscrimination() ) {
	 DoubleUtils.Printf( "MAP SCORES (non-discrim; discrim; total): ", "%.3f ", 
			     new double[] { nonDiscrimF, discrimF, nonDiscrimF + discrimF } );
	 println();
      }
   }

   /*protected void InsertMotifsOLD( Sequence SS[], String insert ) {
     super.InsertMotifs( SS, insert );
     String saveseqs = argProc.getArg( "saveseqs" );
     if ( MyUtils.IsNullString( saveseqs ) ) return;
     String savefname = saveseqs;
     int period = savefname.lastIndexOf( "." );
     if ( period > 0 ) savefname = savefname.substring( 0, period );
     String mots[] = MyUtils.Tokenize( insert, "," );
     if ( mots.length <= 1 ) return;
     int start = 0; if ( mots[ 0 ].length() <= 2 ) start = 1;
     String outSH3s[] = new String[ mots.length ];
     ObjVector outInts[] = new ObjVector[ mots.length ];
     for ( int i = start; i < mots.length; i ++ ) {
     outSH3s[ i ] = Sequence.CreateFakeGeneName();
     outInts[ i ] = new ObjVector();
     }
     for ( int i = 0, s = SS.length; i < s; i ++ ) {
     String seq = SS[ i ].GetSequence();
     for ( int j = start, ss = mots.length; j < ss; j ++ ) {
     if ( seq.indexOf( mots[ j ] ) >= 0 ) outInts[ j ].addElement( SS[ i ].GetHeader() );
     }
     }
     try {
     PrintStream ps = new PrintStream( new FileOutputStream( savefname + ".isSH3Protein.txt" ) );
     for ( int i = start; i < mots.length; i ++ ) ps.println( outSH3s[ i ] + " true" );
     ps = new PrintStream( new FileOutputStream( savefname + ".sif" ) );
     for ( int i = start; i < mots.length; i ++ ) {
     for ( int j = 0, s = outInts[ i ].size(); j < s; j ++ ) {
     ps.println( outSH3s[ i ] + " pp " + outInts[ i ].elementAt( j ) );
     }
     }
     ps = null;
     } catch( Exception e ) { };
     }
   */

   /*protected boolean OkayToCull( Sequence s ) {
     return ! sh3List.contains( s.GetHeader().toUpperCase() ); }*/

   // sh3Ptrs holds, for each interactor, the indexes (in c2[] or intPtrs[]) of its 
   // sh3 partners. NonSh3Ptrs holds, for each interactor, the indexes of all sh3 
   // proteins that it doesnt interact with
   protected void fillSh3Ptrs() {
      sh3Ptrs = new int[ NS ][];
      nonSh3Ptrs = new int[ NS ][];
      for ( int i = 0; i < NS; i ++ ) {
	 IntVector v = new IntVector(), vv = new IntVector();
	 for ( int j = 0; j < nSh3; j ++ ) {
	    int ip[] = intPtrs[ j ]; boolean found = false;
	    for ( int k = 0, size = ip.length; k < size; k ++ ) {
	       if ( ip[ k ] == i ) { v.addElement( j ); found = true; }
	    }
	    if ( ! found ) vv.addElement( j );
	 }
	 sh3Ptrs[ i ] = v.data();
	 nonSh3Ptrs[ i ] = vv.data();
      }
   }

   // intPtrs holds, for each sh3, the indexes (in seqs[] or sh3Ptrs[]) of its interactors
   // nonIntPtrs holds, for each sh3, the indexes of all sequences it doesn't interact with
   protected void fillIntPtrs() {
      intPtrs = new int[ nSh3 ][]; 
      nonIntPtrs = new int[ nSh3 ][]; 
      for ( int i = 0; i < nSh3; i ++ ) {
	 ObjVector v = (ObjVector) sh3s.get( modelIDs[ i ] );
	 IntVector vv = new IntVector(), vv2 = new IntVector();
	 for ( int j = 0, size = v.size(); j < size; j ++ ) {
	    String intr = (String) v.elementAt( j );
	    for ( int ii = 0; ii < NS; ii ++ ) {
	       if ( intr.equals( S[ ii ].GetHeader().toUpperCase() ) ) {
		  vv.addElement( ii ); break; }
	    }
	 }
	 intPtrs[ i ] = vv.data();
	 for ( int ii = 0; ii < NS; ii ++ )
	    if ( ! IntUtils.Contains( intPtrs[ i ], ii ) ) vv2.addElement( ii );
	 nonIntPtrs[ i ] = vv2.data();
      }
   }

   /** Get the a hashtable of vectors of sh3 partner names for each sh3 interactor */
   protected Map getInteractorSh3s( Map sh3s ) {
      Map out = new java.util.HashMap();
      Iterator keys = sh3s.keySet().iterator();
      while( keys.hasNext() ) {
	 String sh3 = (String) keys.next();
	 ObjVector v = (ObjVector) sh3s.get( sh3 );
	 for ( int i = 0, size = v.size(); i < size; i ++ ) {
	    String intr = (String) v.elementAt( i );
	    ObjVector vv = (ObjVector) out.get( intr );
	    if ( vv == null ) vv = new ObjVector();
	    if ( ! vv.contains( sh3 ) ) vv.addElement( sh3 );
	    out.put( intr, vv );
	 }
      }
      return out;
   }

   protected ObjVector getSh3List( String sh3ListName ) {
      ObjVector out = new ObjVector();
      try {
	 Enumeration lines = MyUtils.ReadFileLines( sh3ListName ).elements();
         while( lines.hasMoreElements() ) {
            String line = (String) lines.nextElement();
            String toks[] = MyUtils.Tokenize( line, " " );
	    toks[ 0 ] = replaceSynonym( toks[ 0 ] );
	    if ( "true".equals( toks[ 1 ] ) ) out.addElement( toks[ 0 ] );
	 }
      } catch( Exception e ) { };
      return out;
   }

   /** Get the a hashtable of vectors of interactor names for each SH3 */
   protected Map getSh3Interactors( Map sif ) {
      Map out = new java.util.HashMap();
      try {
	 for ( int i = 0, size = sh3List.size(); i < size; i ++ ) {
	    String sh3 = replaceSynonym( (String) sh3List.elementAt( i ) );
	    ObjVector v = (ObjVector) sif.get( sh3 );
	    if ( v != null ) out.put( sh3, v );
	    else System.err.println( "SH3 protein " + sh3 + " doesn't have any interactors." );
	 }
      } catch( Exception e ) { };
      return out;
   }

   protected Map readSif( String fname, boolean bothWays ) {
      Map interaxns = new java.util.HashMap();
      try {
	 ObjVector vlines = MyUtils.ReadFileLines( fname );
         Enumeration lines = vlines.elements();
	 int cullInds[] = null;
	 if ( argProc.getIntArg( "cull-ints" ) > 0 ) {
	    if ( culledIntsHash == null ) culledIntsHash = new java.util.HashMap();
	    int nlines = vlines.size();
	    boolean flags[] = BoolUtils.New( nlines );
	    int cullInts = argProc.getIntArg( "cull-ints" );
	    cullInds = IntUtils.New( cullInts );
	    for ( int i = 0; i < cullInts; i ++ ) {
	       int ind = -1;
	       while( ind == -1 || flags[ ind ] ) ind = IntUtils.Random( nlines );
	       cullInds[ i ] = ind;
	    }
	 }
	 int lnum = -1;
         while( lines.hasMoreElements() ) {
	    lnum ++;
            String line = (String) lines.nextElement();
            String toks[] = MyUtils.Tokenize( line, " " );
            String p1 = toks[ 0 ].toUpperCase();
	    //String p2 = "pp".equals( toks[ 1 ] ) ? toks[ 2 ].toUpperCase() : toks[ 1 ].toUpperCase();
	    String p2 = toks[ 2 ].toUpperCase();
	    p1 = replaceSynonym( p1 );
	    p1 = p1.toUpperCase();
	    p2 = replaceSynonym( p2 );
	    p2 = p2.toUpperCase();
	    if ( cullInds != null && IntUtils.Contains( cullInds, lnum ) ) {
	       if ( verbose ) println( "Culling interaction " + p1 + " -> " + p2 + 
				       " from the dataset." );
	       culledIntsHash.put( p1, p2 );
	       continue;
	    }
	    if ( culledHash != null && 
		 ( culledHash.containsKey( p1 ) || culledHash.containsKey( p2 ) ) ) {
	       if ( culledIntsHash == null ) culledIntsHash = new java.util.HashMap();
	       culledIntsHash.put( p1, p2 );
	       continue;
	    } 
	    if ( ! seqsHash.containsKey( p1 ) ) {
	       System.err.println( "WARNING: no sequence for protein " + p1 + ". Skipping." );
	       if ( ! seqsHash.containsKey( p2 ) )
		  System.err.println( "WARNING: no sequence for protein " + p2 + ". Skipping." );
	       continue;
	    } else if ( ! seqsHash.containsKey( p2 ) ) {
	       System.err.println( "WARNING: no sequence for protein " + p2 + ". Skipping." );
	       continue;
	    }
            ObjVector v1 = (ObjVector) interaxns.get( p1 );
            if ( v1 == null ) v1 = new ObjVector();
            if ( ! v1.contains( p2 ) ) v1.addElement( p2 );
            interaxns.put( p1, v1 );
	    if ( bothWays ) {
	       ObjVector v2 = (ObjVector) interaxns.get( p2 );
	       if ( v2 == null ) v2 = new ObjVector();
	       if ( ! v2.contains( p1 ) ) v2.addElement( p1 );
	       interaxns.put( p2, v2 );
	    }
	 }
      } catch( Exception e ) { e.printStackTrace(); }
      return interaxns;
   }

   protected void InsertMotifs( Sequence SS[], String insert ) {
      String mots[] = MyUtils.Tokenize( insert, ", " );
      if ( mots.length <= 1 ) return;
      boolean isNum = true; int num = 1;
      try { num = Integer.parseInt( mots[ 0 ] ); }
      catch( Exception e ) { num = 1; isNum = false; }
      if ( isNum ) {
	 insert = MyUtils.Join( mots, ",", 1, mots.length-1 );
	 mots = MyUtils.Tokenize( insert, ", " );
      }

      int specificity = argProc.getIntArg( "spec", 10 );
      PSSMMotifModel model = PSSMMotifModel.GenerateFromConsensi( insert, specificity );

      double sims[][] = DoubleUtils.New( mots.length, mots.length );
      for ( int i = 0; i < mots.length; i ++ ) {
	 for ( int j = i+1; j < mots.length; j ++ ) {
	    double simScore = model.ComputeSimilarityScore( i, j );
	    sims[ i ][ j ] = sims[ j ][ i ] = simScore * simScore * simScore;
	 }
	 DoubleUtils.MaxNorm( sims );
      }

      String saveseqs = argProc.getArg( "saveseqs" );
      if ( MyUtils.IsNullString( saveseqs ) ) return;
      String savefname = saveseqs;
      int period = savefname.lastIndexOf( "." );
      if ( period > 0 ) savefname = savefname.substring( 0, period );

      String outSH3s[] = new String[ mots.length ];
      ObjVector outInts[] = new ObjVector[ mots.length ];
      boolean isSH3[] = BoolUtils.New( SS.length );
      for ( int i = 0; i < mots.length; i ++ ) {
	 int ind = -1;
	 while( ind < 0 || isSH3[ ind ] ) ind = IntUtils.RandChoose( 0, S.length - 1 );
	 isSH3[ ind ] = true;
	 outSH3s[ i ] = S[ ind ].GetHeader();
	 outInts[ i ] = new ObjVector();
      }

      PrintStream infops = null;
      try {
	 infops = new PrintStream( new FileOutputStream( savefname + ".info" ) );
	 infops.println( "SIMILARITY SCORES:\n" ); 

	 for ( int i = 0; i < mots.length; i ++ ) {
	    for ( int j = i+1; j < mots.length; j ++ ) {
	       infops.println( MyUtils.SPrintf( "%16s ", mots[ i ] ) +
			       MyUtils.SPrintf( "%16s  ", mots[ j ] ) +
			       MyUtils.SPrintf( "%4s ", outSH3s[ i ] ) +
			       MyUtils.SPrintf( "%4s ", outSH3s[ j ] ) +
			       DoubleUtils.SPrintf( "%.3f ", sims[ i ][ j ] ) +
			       DoubleUtils.SPrintf( "%.3f", 0.1/sims[ i ][ j ] ) );
	    }
	 }
      } catch( Exception e ) { e.printStackTrace(); }

      // Set up the distrib of SH3 parners per interactor -- scale free distrib: k^(-2)
      double sh3partners_per_interactor[] = DoubleUtils.New( mots.length + 1 );
      for ( int i = 1; i <= mots.length; i ++ ) 
	 sh3partners_per_interactor[ i ] = Math.pow( (double) i, -2.0 );
      // We'll use a flat distribution for the distrib of interactors per SH3,
      // since that's what we see.

      boolean USE_MIXTURE = true; // Use one mixture motif for all sh3 interactors (true), 
      // or else insert a different motif at a different site for each one

      double nints[] = DoubleUtils.New( mots.length );
      double nints2[] = DoubleUtils.New( mots.length );
      DoubleUtils.Random( nints );
      DoubleUtils.MaxNorm( nints );

      ObjVector insertedInfo = new ObjVector();
      boolean chosen[] = BoolUtils.New( mots.length );
      for ( int i = 0, s = SS.length; i < s; i ++ ) {
	 Sequence ss = SS[ i ];
	 String seq = ss.GetSequence();
	 int slen = ss.GetLength();
	 String head = ss.GetHeader();
	 boolean indchosen[] = BoolUtils.New( slen );
	 int nint = 0, lastMot = -1;
	 while( nint == 0 ) nint = DoubleUtils.Sample( sh3partners_per_interactor );
	 BoolUtils.Zero( chosen );
	 IntVector ints = null;
	 if ( USE_MIXTURE ) ints = new IntVector();
	 for ( int j = 0; j < nint; j ++ ) {
	    int mot = -1, ind = -1, tries = 0;
	    DoubleUtils.Copy( nints2, nints );
	    //System.out.print("HERE1 ");DoubleUtils.Printf("%.3f ",nints2);
	    if ( lastMot >= 0 ) {
	       DoubleUtils.Mult( nints2, sims[ lastMot ] );
	       DoubleUtils.MaxNorm( nints2 );
	       //System.out.print("HERE2a ");DoubleUtils.Printf("%.3f ",sims[lastMot]);
	       //System.out.print("HERE2 "+lastMot+" ");DoubleUtils.Printf("%.3f ",nints2);
	       //System.exit(0);
	    }
	    while( tries ++ < 1000 && mot < 0 || chosen[ mot ] || 
		   ( isSH3[ i ] && head.equals( outSH3s[ mot ] ) ) ) // No self-interaxns
	       //mot = IntUtils.RandChoose( 0, mots.length - 1 );
	       mot = DoubleUtils.Sample( nints );
	    if ( tries >= 1000 ) continue;
	    chosen[ mot ] = true;
	    lastMot = mot;
	    if ( USE_MIXTURE ) {
	       ints.add( mot );
	    } else {
	       short res[] = model.GenerateSequence( mot );
	       while( ind < 0 || indchosen[ ind ] )
		  ind = IntUtils.RandChoose( 0, slen - res.length );
	       for ( int ii = ind - mots[ mot ].length(); 
		     ii < ind + mots[ mot ].length(); ii ++ )
		  if ( ii >= 0 && ii < indchosen.length ) indchosen[ ii ] = true;
	       Sequence news = new Sequence( res );
	       ss.Replace( news, ind );
	       outInts[ mot ].addElement( head );
	       insertedInfo.add( head + " " + ind + " " + 
				 news.GetSequence() + "\t" + outSH3s[ mot ] + " " + 
				 mots[ mot ] );
	    }
	 }
	 if ( USE_MIXTURE ) {
	    int data[] = ints.data();
	    short res[] = model.GenerateSequence( data );
	    int ind = -1;
	    while( ind < 0 || indchosen[ ind ] )
	       ind = IntUtils.RandChoose( 0, slen - res.length );
	    indchosen[ ind ] = true;
	    Sequence news = new Sequence( res );
	    ss.Replace( news, ind );
	    for ( int j = 0; j < data.length; j ++ ) {
	       int mot = data[ j ];
	       outInts[ mot ].addElement( head );
	       insertedInfo.add( head + " " + ind + " " + 
				 news.GetSequence() + "\t" + outSH3s[ mot ] + " " + 
				 mots[ mot ] );
	    }
	 }
      }

      try {
	 infops.println( "\nSH3s and their motifs:" );
	 for ( int i = 0; i < mots.length; i ++ ) 
	    infops.println( outSH3s[ i ] + " " + mots[ i ] + 
			    DoubleUtils.SPrintf( " %.3f ", nints[ i ] ) );
	 infops.println( "\nSequences and the inserted peptides:" );
	 for ( int i = 0; i < insertedInfo.size(); i ++ ) 
	    infops.println( insertedInfo.get( i ) );

	 PrintStream ps = 
	    new PrintStream( new FileOutputStream( savefname + ".isSH3Protein" ) );
	 for ( int i = 0; i < mots.length; i ++ ) ps.println( outSH3s[ i ] + " true" );

	 ps = new PrintStream( new FileOutputStream( savefname + ".sif" ) );
	 for ( int i = 0; i < mots.length; i ++ ) {
	    for ( int j = 0, s = outInts[ i ].size(); j < s; j ++ ) {
	       ps.println( outSH3s[ i ] + " pp " + outInts[ i ].elementAt( j ) );
	    }
	 }
	 ps = null;
      } catch( Exception e ) { };
   }

   protected static void implantPxxP( double counts[][] ) {
      int ww = counts.length, centre = ww / 2, start = centre - 2;
      double sum = DoubleUtils.Sum( counts[ start ] );
      counts[ start ][ 12 ] += sum;
      counts[ start + 3 ][ 12 ] += sum;
      counts[ start - 1 ][ 12 ] = counts[ start + 1 ][ 12 ] = 0;
      counts[ start + 2 ][ 12 ] = counts[ start + 4 ][ 12 ] = 0;
   }

   public static String replaceSynonym( String p1 ) {  // A hack to replace synonyms...
      if ( "YHR114W".equals( p1 ) ) return "BZZ1";
      //if ( "YHR016C".equals( p1 ) ) return "YSC84";
      if ( "YSC84".equals( p1 ) ) return "YHR016C";
      if ( "YBR098W".equals( p1 ) ) return "MMS4";
      //if ( "YHL027W".equals( p1 ) ) return "RIM101";
      return p1;
   }

   public void SetupArgs( ArgProcessor argProc ) {
      super.SetupArgs( argProc );
      argProc.AddArg( "NetSiteSampler Parameters" );
      argProc.AddArg( "sif", "<fname>", "", "SIF containing interactions" );
      argProc.AddArg( "sh3", "<fname>", "", "file with list of interaction domain nodes" );
      argProc.AddArg( "non", "<fname>", "", "FASTA with non-interacting sequences" );
      argProc.AddArg( "spec", "<int>", "10", "specificity for simulated motifs" );
      argProc.AddArg( "P2", "<frange:0:100>", "0.5", "supermodel pseudocount fraction" );
      argProc.AddArg( "PD", "<frange:0:100>", "0.5", "discrimination pseudofrequency" );
      argProc.AddArg( "pxxp", null, "false", "pre-implant a PxxP into each motif model" );
   }

   public static void main( String argv[] ) {
      NetSiteSampler samp = new NetSiteSampler( argv );
      if ( ! samp.IsInitialized() ) return;
      samp.Run();
   }
}
