package djr.motif.sampler;
import java.util.*;
import java.io.*;
// START GUI
import java.awt.*;
import java.awt.event.*;
import javax.swing.*;
import djr.motif.gui.*;
import djr.util.gui.*;
// END GUI
import cern.jet.random.Poisson;
import corejava.Format;

import djr.util.array.*;
import djr.util.*;
import djr.util.bio.*;
import djr.motif.model.prior.*;
import djr.motif.model.*;
import djr.motif.sampler.prior.*;
import djr.motif.sampler.sites.*;

/**
 * <code>Sampler</code> class.
 *
 * @author <a href="mailto:reiss@uw.edu">David Reiss</a>
 * @version 1.9978 (Fri Nov 07 05:56:26 PST 2003)
 */
public class Sampler implements ArgProcessorInterface, Runnable, Serializable 
// START GUI
				, ActionListener
// END GUI
{
   static final boolean MAPSUMSQUARED = true; // "true" seems to work better.
   public static final String PROGRESS_CHANGE = "PROGRESS_CHANGE";
   public static final String PROGRESS_COMPLETE = "PROGRESS_COMPLETE";
   public static final String SEQUENCE_HEADER_COMPLEMENT = "(COMPLEMENT)";
   public static final String SEQUENCE_HEADER_PLUS_COMPLEMENT = "(plus COMPLEMENT)";
   
   public String file, save, modelIDs[], localHostName;
   public String logos, postProcessors, sitesPriorsString;
   public double cull, smallIter, finalThresh = 0;
   public int doRankTest, onRankTest, smallIterSize;
   public boolean quiet, verbose, vverbose, initialized, noSamp, isRestored;
   public boolean ranking, ansi, gui, html, display, askGuiOpt, noExit, threaded;
   public boolean signalPause, signalPrintAlignments, signalCancel, inNearOptimumSampling;

   public long seed, startTime, myID = System.currentTimeMillis();
   public int iter1, iter2, nworse, maxIter1, maxIter2, nbad, testing, debug;
   public int nearOptimumIter1, nearOptimumIter2, runMultiples, onMultiple;

   public Sequence S[], culled[];
   public Map selfs = null, globals = null;
   //public double convPseudoFrac, hydroPrior;
   //public int convWidth, hydroWidth; 
   public ArgProcessor argProc = null;
   public Map culledHash = null;

   public AlignmentMotifModel mmodel;
   public ObjVector sitesPriors;

   public String fgType, dirichName, bgfile, matrix; //, Wstring;

   // Model-dependent params (should be moved to other model-based classes [see SiteSampler]):
   public double pseudo, power;
   public int nMotifs, W, bgorder, minPerSeq, maxPerSeq;

   public transient ObjVector savedStates = null;
   // START GUI
   transient ObjVector listeners = null;
   transient javax.swing.JFrame logoFrame = null;
   // END GUI

   public short J;
   public int NS, len[], maxLen;

   public double qx[]; // Temporary array kept here for speed

   public transient Thread mainThread; 
   public transient PrintStream pout = System.out;

   public Sampler() { }
   public Sampler( String argv[] ) { Initialize( argv ); }
   public Sampler( String args ) { this( MyUtils.Tokenize( args, " " ) ); }

   protected boolean Initialize( String args ) { 
      return Initialize( MyUtils.Tokenize( args, " " ) ); }

   protected boolean Initialize( String argv[] ) {
      int removed = 0, numculled = 0;
      if ( ! ranking && onMultiple <= 0 ) {

	 // START GUI
	 if ( listeners != null ) 
	    for ( int li = 0, size = listeners.size(); li < size; li ++ ) 
	       ( (ActionListener) listeners.elementAt( li ) ).
		  actionPerformed( new ActionEvent( this, 0, PROGRESS_CHANGE ) );
	 // END GUI

	 if ( argv != null ) {
	    if ( ! askGuiOpt ) {
	       argProc = new ArgProcessor();
	       if ( ! argProc.ProcessArgs( this, argv ) ) MyUtils.Exit( 0 );
	    }
	    // START GUI	
	    if ( askGuiOpt || argv == null || argv.length <= 0 ) {
	       argProc = new ArgProcessorGUI();
	       if ( ! argProc.ProcessArgs( this, argv ) ) MyUtils.Exit( 0 );
	       argv = ( (ArgProcessorGUI) argProc ).getNewArgs();
	    }
	    // END GUI
	 }
	 if ( vverbose ) verbose = true;

	 String tmpFl = argProc.getArg( "tmp" );
	 if ( ( ! MyUtils.IsNullString( tmpFl ) ) ) {
	    File tmpF = new File( tmpFl );
	    if ( tmpF.exists() ) System.exit( 1 );
	    tmpF.deleteOnExit();
	    try {
	       FileOutputStream tfo = new FileOutputStream( tmpF );
	       tfo.write( 0 ); tfo.flush(); tfo.close();
	    } catch ( Exception e ) {
	       System.err.println( "ERROR: Could not create temporary file " + tmpFl );
	    }
	 }

	 if ( gui ) { InitializeGUI(); }
	 else if ( html ) {
	    ansi = gui = false;
	    print( HTML.begin( GetHeaderString() ) );
	    print( HTML.font( "courier" ) );
	 }
	 String log = argProc.getArg( "log" );
	 if ( ! MyUtils.IsNullString( log ) ) {
	    try { pout = new LoggedSystemOut( log ); }
	    catch( Exception e ) { 
	       if ( ! quiet ) System.err.println( "Cannot open logfile " + log + " for output" );
	    }
	 }

	 DoubleUtils.SetSeed( seed );

	 String fake = argProc.getArg( "fake" );
	 if ( ! MyUtils.IsNullString( file ) ) {
	    try {
	       S = Sequence.ReadSequences( file );
	    } catch( Exception e ) {
	       if ( ! quiet ) { System.err.println( "File '" + file + "' not found. Exiting.\n" );
	       e.printStackTrace(); }
	       MyUtils.Exit( 0 );
	    }
	 } else if ( ! MyUtils.IsNullString( fake ) ) {
	    S = Sequence.CreateFakeSequences( fake, pseudo );
	 } 
      
	 if ( S == null ) {
	    if ( ! quiet ) System.err.println( "ERROR: Cannot open sequence file " + file + " to read." );
	    MyUtils.Exit( 1 );
	    return false;
	 }

	 NS = S.length;
	 if ( argProc.getBooleanArg( "rd" ) ) {
	    S = Sequence.RemoveDups( S, vverbose );
	    removed = NS - S.length;
	 }
	 NS = S.length;
	 if ( cull != 0.0 ) {
	    Sequence[][] out = CullSeqs( S );
	    S = out[ 0 ];
	    culled = out[ 1 ];
	    numculled = culled.length;
	    if ( verbose ) {
	       for ( int i = 0 ; i < culled.length; i ++ )
		  println( "Culled " + culled[ i ].GetHeader() + " from the sample." );
	    }
	 }
	 NS = S.length;

	 if ( argProc.getBooleanArg( "shuffle" ) ) ShuffleSeqs( S );
	 String insert = argProc.getArg( "insert" );
	 if ( ! MyUtils.IsNullString( insert ) ) InsertMotifs( S, insert );
	 PostProcessSeqs( S );

	 if ( argProc.getBooleanArg( "wc" ) ) {
	    Sequence newS[] = new Sequence[ NS * 2 ];
	    for ( int i = 0, j = 0; j < NS; j ++ ) {
	       newS[ i ++ ] = S[ j ];
	       newS[ i ++ ] = new Sequence( S[ j ].GetComplement(), S[ j ].GetHeader() + 
					    " " + SEQUENCE_HEADER_COMPLEMENT );
	    }
	    S = newS;
	    NS = S.length;
	 } else if ( argProc.getBooleanArg( "worc" ) ) {
	    for ( int j = 0; j < NS; j ++ ) {
	       S[ j ].Append( new Sequence( S[ j ].GetComplement(), S[ j ].GetHeader() + 
					    " " + SEQUENCE_HEADER_PLUS_COMPLEMENT ) );
	    }
	 }

	 if ( argProc.getBooleanArg( "prot" ) ) 
	    for ( int j = 0; j < NS; j ++ ) S[ j ].ToProtein();

	 String saveseqs = argProc.getArg( "saveseqs" );
	 if ( ! MyUtils.IsNullString( saveseqs ) ) {
	    Sequence.WriteSequences( saveseqs, S ); MyUtils.Exit( 0 ); }
      } else if ( initialized && ranking ) {
	 SetupSeqsForRanking();
      }

      //if ( W == null ) W = IntUtils.FTokenize( Wstring, "," );
      nMotifs = 1; //W.length;

      J = S[ 0 ].GetAlphabetSize();
      for ( int i = 0; i < NS; i ++ ) {
	 if ( S[ i ].GetAlphabetSize() != J ) {
	    if ( ! quiet ) System.err.println( "ERROR: Sequence number " + (i+1) + " (" +
					       S[ i ].GetHeader() + ") is a " + S[ i ].GetTypeName() +
					       " but other sequences are not.\nExiting ungracefully." );
	    MyUtils.Exit( 1 );
	 }
      }
      len = IntUtils.New( NS );
      maxLen = Integer.MIN_VALUE;
      for ( int i = 0; i < NS; i ++ ) {
	 int llen = S[ i ].GetLength();
	 len[ i ] = llen;
	 if ( llen > maxLen ) maxLen = llen;
      }

      Sites sites = new Sites( S, W, nMotifs, minPerSeq, maxPerSeq );
      qx = DoubleUtils.New( maxLen );

      mmodel = new AlignmentMotifModel( nMotifs, sites, W );
      mmodel.SetBackgroundAndForegroundModels( bgfile, bgorder, pseudo, fgType, 
					       dirichName, matrix );

      LoadSitesPriors( sitesPriorsString );

      if ( verbose && ! quiet ) {
	 if ( html ) println( HTML.name( "top" ) );
	 if ( ! ranking && vverbose ) {
	    for ( int i = 0; i < NS; i ++ ) {
	       if ( S[ i ] == null ) continue;
	       println( "SEQUENCE " + (i+1) + " " + S[ i ].GetHeader() + " (l=" + S[ i ].GetLength() 
			+ "): " + S[ i ].GetSequence() );
	    }	
	 } 
	 if ( ranking ) {
	    if ( html ) println( HTML.BR );
	    if ( ! quiet ) println( "\nPerforming Wilcoxon rank siginificance test #" + onRankTest +
				    " (" + maxIter1 + " iterations)..." );
	 } else if ( onMultiple <= 0 ) {
	    println( "Read " + (NS+removed+(culled!=null?culled.length:0)) + " sequences." );
	    if ( removed > 0 ) println( "Removed " + removed + " duplicate sequences." );
	    if ( numculled > 0 ) println( "Culled an additional " + numculled + " sequences." );
	    println( "Running on " + NS + " remaining sequences." );
	    println( "Alphabet size is " + J + "." );
	    println( "A total of " + IntUtils.Sum( len ) + " residues." );
	    printf( "Average sequence length: %.1f ", IntUtils.Mean( len ) );
	    printf( "+/- %.1f\n", IntUtils.Stddev( len ) );
	 }
      }

      modelIDs = new String[ nMotifs ];
      for ( int i = 0; i < nMotifs; i ++ ) modelIDs[ i ] = (i+1) + "";

      if ( ! ranking ) {
	 Date d = new Date( System.currentTimeMillis() );
	 if ( ! quiet && onMultiple <= 0 ) {
	    println( ReflectUtils.getClassName( this ) + " version " + 
		     DoubleUtils.version + "\n" );
	    println( d + "\n\nCOMMAND LINE: java " + 
		     ReflectUtils.getFullClassName( this ) + " " + 
		     argProc.getFinalArgs() );
	    println( "\nINPUT PARAMETERS:\n" + argProc.getFinalParams() );
	    println( "ISEED = " + seed ); }
      }
      startTime = System.currentTimeMillis();
      localHostName = MyUtils.GetLocalHostName();
      initialized = true;
      return true;
   }

   public void SetSequences( Sequence SS[] ) {
      this.S = SS;
      Initialize( (String[]) null );
   }

   public void SetSequences( String argString, Sequence SS[] ) {
      this.S = SS;
      Initialize( argString );
   }

   public void Run() {
      run();
   }

   public void run() { // Runnable interface method
      if ( ! threaded ) RunMultiples( runMultiples );
      else {
	 Thread thr = new Thread( this, "sampler" );
	 threaded = false;
	 thr.start();
      }
      if ( ! noExit && ! display && ! "gui".equals( logos ) ) MyUtils.Exit( 0 );
   }

   public void SetParameters( String args ) {
      Initialize( MyUtils.Tokenize( args, " " ) ); }

   protected void IterateSampler( int niter ) { };

   protected void AddRemoveFromCounts( boolean add, int seqNum, int mot ) {
      mmodel.AddRemoveFromCounts( add, seqNum, mot );
   }

   public Sites GetSites() { return mmodel.GetSites(); }

   // Make sure all un-serialized bits are given valid pointers
   private void readObject( ObjectInputStream s ) throws IOException, 
							 ClassNotFoundException {
      s.defaultReadObject();
      GetSites().SetSequences( this.S );
   }

   protected Map InitializeState( boolean self ) {
      Map out = new java.util.HashMap();
      if ( self ) {
	 GetSites().SetSequences( S );
	 mmodel.Initialize();
	 out.put( "mm", mmodel );
	 out.put( "raw", new Boolean( true ) );
      } else {
	 out.put( "mm", mmodel.Duplicate() );
      }
      return out;
   }

   protected void CopyState( Map to, Map from ) {
      ( (AlignmentMotifModel) to.get( "mm" ) ).
	 CopyFrom( (AlignmentMotifModel) from.get( "mm" ) );
      //to.put( "mm", MyUtils.DeepCopy( from.get( "mm" ) ) );
   }   

   public void SetState( Map newState ) {
      CopyState( selfs, newState );
   }

   public boolean IsInitialized() { return initialized; }
   public AlignmentMotifModel GetMotifModel() { return mmodel; }
   public ObjVector GetSavedStates() { return savedStates; }
   public ArgProcessor GetArgProcessor() { return argProc; }
   public void SetNoExit( boolean noe ) { noExit = noe; MyUtils.SetNoExit( noe ); }
   public Map GetState() { return selfs; }
   public void SetNoExit() { noExit = true; }
   public PrintStream GetOutput() { return pout; }
   public void SetNearOptimumIters( int ni1, int ni2 ) {
      nearOptimumIter1 = ni1; nearOptimumIter2 = ni2; }
   public void SetVerbosity( boolean v, boolean vv, boolean q ) {
      verbose = v; vverbose = vv; quiet = q; }

   protected void RunSampler() {
      if ( ! initialized ) return;
      mainThread = Thread.currentThread();

      selfs = InitializeState( true );
      globals = InitializeState( false );

      Map locals = InitializeState( false );

      iter1 = 0;
      while( iter1 ++ < maxIter1 ) {
	 if ( signalCancel ) return;
	 if ( signalPause ) {
	    println( "Pausing. Press \"Pause\" button again to resume." );
	    while ( signalPause ) {
	       pout.flush();
	       try { mainThread.sleep( 1000 ); }
	       catch( Exception e ) { };
	    }
	    println( "Resuming..." );
	 }
	 PrintIters( iter1 + "" );

	 ShuffleSites();
	 FillTheArrays();

	 locals.put( "raw", new Boolean( true ) );
	 InnerIter( locals );

	 double sum1 = ComputeMAPScore( locals );
	 double sum2 = ComputeMAPScore( globals );
	 //System.err.println( "\nHERE1: old=" + sum2 + " new=" + sum1 );
	 if ( sum1 != 0.0 && ( Double.isInfinite( sum2 ) || sum1 > sum2 ) ) {
	    CopyState( globals, locals );
	    if ( verbose ) {
	       SetState( globals );
	       PrintSummary( globals, vverbose );
	    }
	    //SendResults( globals );
	 }

	 // START GUI
	 if ( listeners != null ) {
	    int level = onMultiple * maxIter1 + iter1 + 1;
	    level = (int) ( ( 100.0 * level ) / 
			    (double) ( maxIter1 * ( runMultiples + 1 ) ) );
	    for ( int li = 0, size = listeners.size(); li < size; li ++ ) 
	       ( (ActionListener) listeners.elementAt( li ) ).
		  actionPerformed( new ActionEvent( this, level, PROGRESS_CHANGE ) );
	 }
	 // END GUI
      }

      if ( savedStates == null ) savedStates = new ObjVector();
      savedStates.addElement( globals );

      this.iter1 = this.iter2 = 0;
      PerformNearOptimumSampling();
      this.nearOptimumIter1 = this.nearOptimumIter2 = 0;
      PrintResults();
   }

   protected void RunMultiples( int nRepeats ) {
      Sequence newS[] = new Sequence[ S.length ];
      for ( int i = 0; i < S.length; i ++ ) newS[ i ] = new MaskedSequence( S[ i ] );
      S = newS;

      RunSampler();      
      for ( int i = 0; i < nRepeats; i ++ ) {
	 MaskAlignment( globals );
	 //for ( int j = 0; j < S.length; j ++ ) ( (MaskedSequence) S[ j ] ).ApplyMask();
	 onMultiple = i + 1;
	 boolean out = Initialize( MyUtils.Tokenize( argProc.getFinalArgs(), " " ) );
	 RunSampler();
      }

      PrintFinalResults( false, false );

      if ( ! quiet && ! ranking ) {
	 println( "ISEED WAS " + seed );
	 double secs = (double) ( System.currentTimeMillis() - startTime ) / 1000.0;
	 printf( "Completed in %.2f seconds.\n", secs );
      }

      // START GUI
      if ( ! signalCancel && listeners != null ) 
	 for ( int li = 0, size = listeners.size(); li < size; li ++ ) 
	    ( (ActionListener) listeners.elementAt( li ) ).
	       actionPerformed( new ActionEvent( this, 0, PROGRESS_COMPLETE ) );
      // END GUI

      CloseUpShop();
   }

   public void CloseUpShop() {
      if ( doRankTest > 0 ) PerformRankTest();
      if ( ! MyUtils.IsNullString( postProcessors ) ) RunPostProcessors( postProcessors );
   }

   public void PrintIters( String iters ) {
      if ( ! quiet && ! ranking && ( verbose || ranking ) ) {
	 if ( ansi ) print( ANSI.leftCursor( 99999 ) + 
			    ANSI.clearLine() + "ITER = " + iters );
	 else if ( ! ansi || ! gui ) println( "ITER = " + iters );
	 // START GUI
	 if ( gui && html ) ( (JPrintStreamPanel) pout ).setStatus( "ITER = " + iters );
	 else if ( gui ) ( (PrintStreamPanel) pout ).setStatus( "ITER = " + iters );
	 // END GUI
      }
   }

   public void PrintResults() {
      if ( ranking || quiet ) return;
      SetState( globals );
      if ( ! quiet ) PrintSummary( globals, vverbose || runMultiples <= 0 );
   }

   protected void PrintSummary( Map state, boolean printAlignments ) {
      SamplerOutput.PrintSummary( this, state, printAlignments, true ); }

   public void PrintBestMotif( Map state, String name/*int mot*/ ) {
      println(); if ( html ) print( HTML.HR + HTML.BR );
      println( "Info for motif " + name + ":" ); println();
      SetState( state );
      SamplerOutput.PrintSummary( this, state, ! ranking, false );

      String gffName = argProc.getArg( "gff" );
      if ( ! MyUtils.IsNullString( gffName ) )
	 SamplerOutput.SaveAlignmentToGFF( this, gffName, true );
   }

   public void PrintBestMotifs( ObjVector states ) {
      for ( int j = 0, size = states.size(); j < size; j ++ ) {
	 Map state = (Map) states.elementAt( j );
	 PrintBestMotif( state, "#"+( j + 1 ) );
      }
      if ( ! quiet ) println( "\n" + new Date( System.currentTimeMillis() ) );
   }

   public void PrintFinalResults( boolean force, boolean forceDisplay ) {
      if ( force ) { quiet = false; verbose = vverbose = true; }
      if ( forceDisplay ) display = true;
      if ( ! quiet ) println( "\n\nDONE...  " );
      // START GUI
      if ( gui && html ) ( (JPrintStreamPanel) pout ).setStatus( "DONE!" );
      else if ( gui ) ( (PrintStreamPanel) pout ).setStatus( "DONE!" );
      // END GUI

      String gffName = argProc.getArg( "gff" );
      if ( ! MyUtils.IsNullString( gffName ) ) {
	 try {
	    File gff = new File( gffName );
	    if ( gff.exists() ) gff.delete();
	 } catch( Exception e ) { };
      }

      if ( ! quiet ) PrintBestMotifs( savedStates );

      // START GUI
      if ( display )
	 ( new MotifDisplayPanel( "Best motifs", savedStates, this ) ).putInFrame().QuitOnClose();
      // END GUI

      String saveStateFile = force ? null : argProc.getArg( "S" );
      if ( ! isRestored && ! MyUtils.IsNullString( saveStateFile ) ) 
	 SaveStateToFile( saveStateFile );
   }

   public int ComputeProbsForSequence( short seq[], double q[][], double[] qx, 
				       boolean normalize ) {
      int whereMax = mmodel.ComputeProbsForSequence( seq, q, qx, normalize );
      int whereMax2 = AdjustProbsForSequence( qx, seq, -1, normalize );
      if ( whereMax2 >= 0 ) whereMax = whereMax2;
      return whereMax;
   }

   public int ComputeProbsForSequence( int seqNum, int mot, double[] qx, 
				       boolean normalize ) {
      short seq[] = GetSites().GetSequenceResidues( seqNum );
      int whereMax = mmodel.ComputeProbsForSequence( seq, mot, qx, normalize );
      int whereMax2 = AdjustProbsForSequence( qx, seq, mot, normalize );
      if ( whereMax2 >= 0 ) whereMax = whereMax2;
      return whereMax;
   }

   public int AdjustProbsForSequence( double qx[], short seq[], int mot, 
				      boolean normalize ) {
      int whereMax = -1;
      if ( sitesPriors != null || power != 1.0 ) {
	 double max = -Double.MAX_VALUE;
	 int maxl = seq.length; 
	 for ( int i = 0; i < maxl; i ++ ) {
	    qx[ i ] = AdjustProbForSite( qx[ i ], seq, i, mot );
	    if ( qx[ i ] > max ) { whereMax = i; max = qx[ i ]; }
	 }
	 if ( normalize ) DoubleUtils.Divide( qx, max );
      }
      return whereMax;
   }

   public int AdjustProbsForSequence( double qx[], int seqNum, int mot, 
				      boolean normalize ) {
      return AdjustProbsForSequence( qx, GetSites().GetSequenceResidues( seqNum ), mot, normalize );
   }

   public double AdjustProbForSite( double prob, short seq[], int site, int mot ) {
      if ( power != 1.0 ) prob = Math.pow( prob, power );
      if ( sitesPriors != null && sitesPriors.size() > 0 ) {
	 for ( int i = 0, sz = sitesPriors.size(); i < sz; i ++ ) {
	    prob *= ( (SitesPrior) sitesPriors.get( i ) ).
	       GetPriorValue( GetSites(), seq, site, 0 );
	 }
      }
      return prob;
   }

   public double AdjustProbForSite( double prob, int seqNum, int site, int mot ) {
      return AdjustProbForSite( prob, GetSites().GetSequenceResidues( seqNum ), site, mot );
   }

   protected void FillTheArrays() { mmodel.FillTheArrays(); }

   protected void ShuffleSites() { GetSites().ShuffleSites(); }

   protected void ShiftSites( int amount ) { GetSites().ShiftSites( amount ); }

   protected boolean IsSiteValid( int site, int seq ) {
      return GetSites().IsSiteValid( seq, site );
   }

   /** Format of output alignments aa[][][] is mot x num x (alignment length + 1) where
       num is the number of alignments attributed to motif number mot in all sequences.
       The alignment is in the last dimension. First index is the seq. number, and indices
       1->length give the indexes in that sequence that are part of the alignment. (Sorry!) **/
   public int[][][] GetAlignments() {
      if ( finalThresh == 0 ) return GetSites().GetAlignments();
      else return (int[][][]) GetSites().LocateAlignments( this, minPerSeq, maxPerSeq, finalThresh ).get( "a" );
   }

   public double ComputeMAPScore() {
      return ComputeMAPScore( null );
   }

   public double ComputeMAPScore( Map state ) {
      if ( state == null ) state = selfs;
      if ( state.get( "raw" ) != null ) { state.remove( "raw" ); return -Double.MAX_VALUE; }
      AlignmentMotifModel mmod = (AlignmentMotifModel) state.get( "mm" );
      double F = mmod.ComputeMAPScore();
      if ( sitesPriors != null && sitesPriors.size() > 0 ) {
	 ObjVector sitesList = new ObjVector();
	 for ( int ii = 0, sz = sitesPriors.size(); ii < sz; ii ++ ) {
	    SitesPrior sp = (SitesPrior) sitesPriors.get( ii );
	    for ( int mot = 0; mot < nMotifs; mot ++ ) {
	       for ( int i = 0; i < NS; i ++ ) {
		  sitesList.removeAll();
		  short seq[] = mmod.GetSites().GetSequenceResidues( i );
		  mmod.GetSites().GetAlignmentsForSequence( i, mot, sitesList );
		  for ( int j = 0, s = sitesList.size(); j < s; j ++ )
		     F += DoubleUtils.Log2( sp.GetPriorValue( mmod.GetSites(), seq, 
							      ( (int[]) sitesList.elementAt( j ) )[1], mot ) );
	       }
	    }
	 }
      }
      return F;
   }

   protected void InnerIter( Map locals ) {
      iter2 = 0; nworse = 0;
      int sis = smallIterSize; if ( sis <= 0 ) sis = NS / 10;
      while( iter2 ++ < maxIter2 && nworse <= nbad ) {
	 boolean smallIt = smallIter != 1.0 && DoubleUtils.Random() <= ( 1.0 - smallIter );
	 if ( smallIt ) IterateSampler( -1 );
	 else { IterateSampler( sis ); iter2 --; if ( nworse >= nbad ) nworse --; }
	 double sum1 = ComputeMAPScore( selfs );
	 double sum2 = ComputeMAPScore( locals );
	 //System.err.println( "\nHERE2: old=" + sum2 + " new=" + sum1 );
	 if ( sum1 != 0.0 && ( sum1 > sum2 || Double.isInfinite( sum2 ) ) ) {
	    CopyState( locals, selfs );
	    nworse = 0;
	 } else {
	    if ( nworse > nbad ) break;
	    nworse ++;
	 }
      }
   }

   protected void PerformNearOptimumSampling() {
      inNearOptimumSampling = true;
      if ( nearOptimumIter1 > 0 ) {
	 if ( verbose ) println( "\n\nPerforming near-optimum sampling..." );
	 for ( int i = 0; i < savedStates.size(); i ++ ) {
	    Map state = (Map) savedStates.get( i );
	    Map newState = DoNearOptimumSampling( state );
	    double sum1 = ComputeMAPScore( state );
	    double sum2 = ComputeMAPScore( newState );
	    //System.err.println("\nHERE3: "+sum1+" "+sum2);
	    if ( sum2 != 0.0 && ( Double.isInfinite( sum1 ) || sum2 > sum1 ) ) {
	       /*boolean similar = true;
		 if ( simThresh > 0 ) {
		 AlignmentMotifModel m1 = (AlignmentMotifModel) state.get( "mm" );
		 AlignmentMotifModel m2 = (AlignmentMotifModel) newState.get( "mm" );
		 for ( int mot = 0; mot < nMotifs; mot ++ )
		 if ( ! m1.IsSimilarTo( m2, mot, simThresh ) ) similar = false;
		 }
		 if ( similar )*/ CopyState( state, newState );
	    }
	 }
      }
      inNearOptimumSampling = false;
      // return states;
   }

   protected Map DoNearOptimumSampling( Map state ) {
      Map glob = InitializeState( false ), loc = InitializeState( false );
      CopyState( glob, state );

      int noiter1 = 0, saveIter2 = maxIter2;
      maxIter2 = nearOptimumIter2;
      while( noiter1 ++ < nearOptimumIter1 ) {
	 SetState( state );
	 CopyState( loc, state );
	 CopyState( selfs, state );
	 if ( noiter1 >= nearOptimumIter1 ) noSamp = true;
	 loc.put( "raw", new Boolean( true ) );
	 InnerIter( loc );
	 //if ( mmodel.GetMotifModel( 0 ) != null && DoubleUtils.Sum( mmodel.GetMotifModel( 0 ) ) <= 0 )
	 // while ( DoubleUtils.Sum( mmodel.GetMotifModel( 0 ) ) <= 0 ) InnerIter( loc );
	 double sum1 = ComputeMAPScore( loc );
	 double sum2 = ComputeMAPScore( glob );
	 //System.err.println("\nHERE4: "+sum1+" "+sum2);
	 PrintIters( noiter1 + " of " + nearOptimumIter1 );
	 if ( sum1 != 0.0 && ( Double.isInfinite( sum2 ) || sum1 > sum2 ) )
	    CopyState( glob, loc );
      }
      maxIter2 = saveIter2;
      return glob;
   }

   public ObjVector GetBestMotifStates() {
      return savedStates;
   }

   public void SaveStateToFile( String saveStateFile ) {
      if ( ! MyUtils.IsNullString( saveStateFile ) ) {
	 if ( vverbose ) println( "Saving state of sampler to file " + saveStateFile +
				  ".gz" );
	 String args = argProc.getFinalArgs();
	 if ( args.indexOf( "-iseed" ) < 0 ) args += " -iseed " + seed;
	 Map out = new HashMap();
	 out.put( "args", args );
	 out.put( "classname", ReflectUtils.getFullClassName( this ) );
	 out.put( "states", savedStates );
	 MyUtils.SaveObject( out, saveStateFile + ".gz" );
      }
   }

   public void RestoreStateFromFile( String saveStateFile ) {
      if ( ! MyUtils.IsNullString( saveStateFile ) ) {
	 if ( ! saveStateFile.toUpperCase().endsWith( ".GZ" ) ) saveStateFile += ".gz";
	 if ( vverbose ) println( "Restoring state of sampler from file " + 
				  saveStateFile );
	 Map in = (Map) MyUtils.ReadObject( saveStateFile );
	 savedStates = (ObjVector) in.get( "states" );
	 if ( selfs == null ) selfs = InitializeState( true );
	 if ( globals == null ) globals = InitializeState( false );
	 Map state = (Map) savedStates.get( 0 );
	 SetState( state );
	 CopyState( globals, state );
	 isRestored = true;
      }
   }

   protected Sequence[][] CullSeq( Sequence SS[], String toCull ) {
      culledHash = new java.util.HashMap();
      int indToCull = -1; boolean isString = false;
      try { indToCull = Integer.parseInt( toCull ) - 1; }
      catch( Exception e ) { indToCull = -1; isString = true; }
      if ( indToCull < 0 ) {
	 toCull = toCull.toUpperCase();
	 for ( int i = 0, ll = SS.length; i < ll; i ++ ) 
	    if ( SS[ i ].GetHeader().toUpperCase().equals( toCull ) &&
		 OkayToCull( SS[ i ] ) ) { indToCull = i; break; }
      } //else if ( ! OkayToCull( SS[ indToCull ] ) ) indToCull = -1;
      Sequence out1[] = SS, out2[] = null;
      if ( indToCull >= SS.length ) {
	 println( "Sequence requested to cull, #" + indToCull + 
		  " does not exist. Culling nothing." );
      } else if ( indToCull >= 0 ) {
	 int realInd = indToCull;
	 if ( ! isString ) {
	    realInd = 0;
	    for ( int j = 0; j <= indToCull; realInd ++ ) {
	       if ( realInd >= SS.length || realInd < 0 ) { realInd = -1; break; }
	       if ( realInd < SS.length && OkayToCull( SS[ realInd ] ) ) j ++;
	    }
	    realInd --;
	 }
	 if ( realInd < SS.length && realInd >= 0 ) {
	    String head = SS[ realInd ].GetHeader().toUpperCase();
	    println( "Culling requested sequence #" + realInd + ": " + head );
	    culledHash.put( head, new Boolean( true ) );
	    out1 = new Sequence[ SS.length - 1 ]; out2 = new Sequence[ 1 ];
	    out2[ 0 ] = SS[ realInd ];
	    for ( int i = 0, j = 0, ll = SS.length; i < ll; i ++ ) 
	       if ( i != realInd ) out1[ j ++ ] = SS[ i ];
	 }
      }
      return new Sequence[][] { out1, out2 };
   }      

   protected Sequence[][] CullSeqs( Sequence SS[] ) {
      culledHash = new java.util.HashMap();
      int num = 0, ll = SS.length;
      if ( cull < 1.0 ) num = (int) Math.round( (1.0-cull) * ll );
      else num = (int) Math.round( cull );
      boolean flags[] = BoolUtils.New( ll );
      for ( int i = 0; i < num; i ++ ) {
	 int ind = IntUtils.RandChoose( 0, ll - 1 );
	 if ( flags[ ind ] || ! OkayToCull( SS[ ind ] ) ) { i --; continue; }
	 flags[ ind ] = true;
      }
      Sequence out2[] = new Sequence[ num ], out1[] = new Sequence[ SS.length - num ];
      for ( int i = 0, j = 0, k = 0; i < SS.length; i ++ ) {
	 if ( flags[ i ] ) {
	    out2[ j ++ ] = SS[ i ];
	    culledHash.put( SS[ i ].GetHeader().toUpperCase(), new Boolean( true ) );
	 } else out1[ k ++ ] = SS[ i ];
      }
      return new Sequence[][] { out1, out2 };
   }

   protected boolean OkayToCull( Sequence s ) {
      return true;
   }

   protected void ShuffleSeqs( Sequence SS[] ) {
      for ( int i = 0; i < NS; i ++ ) if ( SS[ i ] != null ) SS[ i ].Shuffle();
   }

   protected void InsertMotifs( Sequence SS[], String insert ) {
      int where;
      for ( int i = 0; i < NS; i ++ ) {
	 if ( SS[ i ] == null ) continue;
	 ObjVector inserted = SS[ i ].InsertMotifs( insert );
	 if ( verbose ) {
	    for ( int j = 0; j < inserted.size(); j += 2 ) {
	       println( "INSERTED MOTIF " + inserted.elementAt( j ) + " IN SEQUENCE " + 
			(i+1) + " AT POSITION " + inserted.elementAt( j+1 ) );
	    }
	 }
      }
   }

   protected void SetupSeqsForRanking() {
      if ( ! initialized ) return;
      BackgroundModel bgmod = new BackgroundModel( S, (short) 3, pseudo );
      for ( int i = 0; i < NS; i ++ ) {
	 if ( S[ i ] == null ) continue;
	 short randseq[] = bgmod.GenerateSequence( S[ i ].GetLength() - W/*[ 0 ]*/ );
	 Sequence temp = new Sequence( S[ i ] );
	 temp.Append( new Sequence( randseq ) );
	 S[ i ] = temp;
      }
   }

   protected void PostProcessSeqs( Sequence SS[] ) {
      // For ROWAN (peroxisomes) -- search for " (C)" at end of header (assumes all
      // sequences in fasta file are given as W)
      if ( ! argProc.getBooleanArg( "wc" ) && 
	   ! argProc.getBooleanArg( "worc" ) ) { // Skip if already searching both directions
	 for ( int i = 0; i < NS; i ++ ) {
	    if ( SS[ i ] == null ) continue;
	    String head = SS[ i ].GetHeader();
	    if ( head.endsWith( " (C)" ) ) {
	       SS[ i ] = new Sequence( SS[ i ].GetComplement(), head );
	       if ( vverbose ) println( "Using reverse-complement of sequence " + head );
	    }
	 }
      }

      if ( ! "".equals( argProc.getArg( "cull-seq" ) ) ) {
	 Sequence[][] out = CullSeq( SS, argProc.getArg( "cull-seq" ) );
	 S = out[ 0 ];
	 culled = out[ 1 ];
	 NS = S.length;
      }
   }

   protected void MaskAlignment( Map state ) {
      if ( state == null ) return;
      SetState( state );
      GetSites().MaskAlignment( S, GetAlignments() );
   }

   protected void PerformRankTest() {
      ranking = true;
      maxIter1 *= 2; maxIter2 *= 2;
      nearOptimumIter1 *= 2; nearOptimumIter2 *= 2;
      minPerSeq *= 2;
      Sequence saveS[] = new Sequence[ S.length ];
      System.arraycopy( S, 0, saveS, 0, S.length );
      onRankTest = 0;
      int numRankTests = doRankTest;
      doRankTest = 0;
      while( onRankTest ++ < numRankTests ) {
	 System.arraycopy( saveS, 0, S, 0, S.length );
	 verbose = true;
	 boolean out = Initialize( MyUtils.Tokenize( argProc.getFinalArgs(), " " ) );
	 verbose = vverbose = false;
	 logos = null;
	 RunSampler();

	 ObjVector states = GetBestMotifStates();
	 if ( states != null ) {
	    Map state = (Map) states.elementAt( states.size() - 1 );
	    SetState( state );
	 }

	 double probs[] = DoubleUtils.New( nMotifs );
	 int aa[][][] = GetAlignments();
	 double pvalues[][] = mmodel.ComputePolyScores( aa );

	 int gsum = 0, gokay = 0, gbad = 0;
	 for ( int mot = 0; mot < nMotifs; mot ++ ) {
	    int ns = aa[ mot ].length;
	    double[] los = DoubleUtils.New( ns );
	    //mmodel.ComputeLogOddsScores( GetSites(), mot, los, aa, pvalues[ mot ], alignThresh );
	    for ( int i = 0; i < ns; i ++ ) {
	       double pp = pvalues[ mot ][ i ];
	       los[ i ] = -DoubleUtils.Log10( pp == 0.0 ? 1.0 / Double.MAX_VALUE : pp );
	    }
	    double[] arr = DoubleUtils.New( los );
	    int[] indx = IntUtils.New( ns );
	    int[] irank = IntUtils.New( ns );
	    DoubleUtils.Mult( arr, -1.0 );
	    DoubleUtils.Indexx( arr, indx );
	    DoubleUtils.Mult( arr, -1.0 );
	    //indx = IntUtils.Reverse( indx );
	    int ngood = 0, nbad = 0;
	    int nOkay = ns - DoubleUtils.NEqualTo( los, 0.0 ), icount = 0;
	    for ( int i = 0; i < ns; i ++ ) {
	       int indd = indx[ i ];
	       irank[ indd ] = 0;
	       if ( los[ indd ] == 0.0 ) continue;
	       int seq = aa[ mot ][ indd ][ 0 ];
	       if ( seq < 0 || S[ seq ] == null ) continue;
	       int loc = aa[ mot ][ indd ][ 1 ];
	       irank[ indd ] = nOkay/*ns*/ - icount/*i*/; // - 1;
	       if ( loc > len[ seq ] / 2 - W/*W*/ ) {
		  irank[ indd ] = -irank[ indd ];
		  nbad ++;
	       } else ngood ++;
	       //println(seq+" "+loc+" "+indd+" "+irank[indd]+" "+los[indd]);
	       icount ++;
	    }
	    //double sig = Math.sqrt( (double) ( ns * ( ns + 1 ) * ( 2 * ns + 1 ) / 6 ) );
	    double sig = Math.sqrt( (double) ( nOkay * ( nOkay + 1 ) * ( 2 * nOkay + 1 ) / 6 ) );
	    int isum = IntUtils.Sum( irank );
	    double nsig = ( isum - 0.5 ) / sig;
	    probs[ mot ] = nsig > 0 ? com.imsl.math.Sfun.erf( nsig ) : 0.0;
	    
	    if ( ! quiet ) {
	       print( "FINAL SIGNIFICANCE (motif " );
	       println( modelIDs[ mot ] + "):" );
	       printf( "\t%.3f% of sites found on bogus sequences (" + 
		       nbad + "/" + nOkay +").\n", 
		       ( (double) nbad / (double) nOkay/*ns*/ ) * 100.0 );
	       printf( "\tOverall P-value is %.3f\n\n", probs[ mot ] );
	    }

	    gsum += isum; gokay += nOkay; gbad += nbad;
	 }

	 if ( nMotifs > 1 ) { 
	    double gsig = Math.sqrt( (double) ( gokay * ( gokay + 1 ) * ( 2 * gokay + 1 ) / 6 ) );
	    double gnsig = ( gsum - 0.5 ) / gsig;
	    double gprob = gnsig > 0 ? com.imsl.math.Sfun.erf( gnsig ) : 0.0;
	    
	    if ( ! quiet ) {
	       println( "FINAL GLOBAL SIGNIFICANCE:" );
	       printf( "\t%.3f% of sites found on bogus sequences (" + 
		       gbad + "/" + gokay + ").\n", 
		       ( (double) gbad / (double) gokay ) * 100.0 );
	       printf( "\tOverall P-value is %.3f\n\n", gprob );
	    }
	 }
      }
   }

   protected void LoadSitesPriors( String sitesPriorsString ) {
      if ( MyUtils.IsNullString( sitesPriorsString ) ) return;
      sitesPriors = new ObjVector();
      String sps[] = MyUtils.Tokenize( sitesPriorsString, ":" );
      try {
	 for ( int i = 0; i < sps.length; i ++ ) {
	    SitesPrior proc = (SitesPrior) Class.forName( sps[ i ] ).
	       getConstructor( new Class[] { Sampler.class } ).
	       newInstance( new Object[] { this } );
	 }
      } catch( Exception e ) { e.printStackTrace(); }
   }

   protected void RunPostProcessors( String postProcessors ) {
      String pps[] = MyUtils.Tokenize( postProcessors, ":" );
      try {
	 for ( int i = 0; i < pps.length; i ++ ) {
	    PostProcessor proc = (PostProcessor) Class.forName( pps[ i ] ).
	       getConstructor( new Class[] { Sampler.class } ).
	       newInstance( new Object[] { this } );
	 }
      } catch( Exception e ) { e.printStackTrace(); }
   }

   protected void finalize() {
      if ( ! quiet && html ) print( HTML.end() );
      pout.flush(); 
   }

   public void SetOut( PrintStream out ) {
      pout = out;
      try { System.setOut( pout ); } catch( Exception e ) { }; // Catch for applet
   }

   // START GUI
   public PrintStreamPanel CreatePrintStreamPanel() {
      PrintStreamPanel psw = new PrintStreamPanel(); 
      psw.addButton( "Pause", this );
      if ( ! vverbose && MyUtils.IsNullString( save ) ) psw.addButton( "Print Alignments", this );
      psw.setFont( new java.awt.Font( "Courier", java.awt.Font.PLAIN, 12 ) );
      psw.setForeground( java.awt.Color.white );
      psw.setBackground( java.awt.Color.black );
      SetOut( psw );
      return psw;
   }

   public JPrintStreamPanel CreateJPrintStreamPanel() {
      JPrintStreamPanel psw = new JPrintStreamPanel();
      psw.setContentHTML();
      psw.addButton( "Pause", this );
      if ( ! vverbose && MyUtils.IsNullString( save ) ) psw.addButton( "Print Alignments", this );
      SetOut( psw );
      return psw;
   }
   // END GUI

   protected void InitializeGUI() {
      ansi = false;
      // START GUI
      if ( html ) {
	 JPrintStreamPanel psw = CreateJPrintStreamPanel();
	 psw.putInFrame( GetHeaderString(), 800, 500 );
	 psw.setStatus( "INITIALIZING..." );
      } else {
	 PrintStreamPanel psw = CreatePrintStreamPanel();
	 psw.putInFrame( GetHeaderString(), 800, 500 );
	 psw.setStatus( "INITIALIZING..." );
      }
      // END GUI
   }

   public void DisplayLogos() {
      // START GUI
      if ( ! vverbose ) return;

      String title = "Motif Logos: ";

      String savef = null;
      if ( ! MyUtils.IsNullString( save ) ) savef = IntUtils.SPrintf( save + "_%05d", iter1 );
      title += "Iteration " + iter1 + ( savef != null ? "; " + savef : "" );

      logoFrame = LogoPanel.showLogos( title, mmodel, modelIDs, logoFrame );
      // END GUI
   }

   protected String GetHeaderString() {
      if ( ! file.startsWith( "/tmp/" ) && ! file.endsWith( ".tmp" ) &&
	   file.indexOf( "readseq" ) < 0 ) {
	 return "Results of run of " + ReflectUtils.getClassName( this ) + 
	    " on " + file;
      } else {
	 return "Results of run of " + ReflectUtils.getClassName( this );
      }
   }

   // START GUI
   public JPrintStreamPanel SetUseHTMLOutput() {
      JPrintStreamPanel psw = CreateJPrintStreamPanel();
      quiet = ansi = false;
      verbose = gui = html = true;
      //logos = "text";
      psw.putInFrame( GetHeaderString(), 800, 500 );
      print( HTML.begin( GetHeaderString() ) );
      print( HTML.font( "courier" ) );
      return psw;
   }

   public PrintStreamPanel SetUseTextGuiOutput() {
      PrintStreamPanel psw = CreatePrintStreamPanel();
      psw.putInFrame( GetHeaderString(), 800, 500 );
      quiet = html = false;
      verbose = gui = true;
      //logos = "text";
      return psw;
   }

   public void actionPerformed( ActionEvent evt ) {
      if ( "Print Alignments".equals( evt.getActionCommand() ) ) signalPrintAlignments = true;
      else if ( "Cancel".equals( evt.getActionCommand() ) ) signalCancel = true;
      else if ( "Pause".equals( evt.getActionCommand() ) ) signalPause = ! signalPause;
      else if ( "Close".equals( evt.getActionCommand() ) ) signalPause = ! signalPause;
   }
   // END GUI

   // START GUI
   public void AddActionListener( Object sl ) {
      if ( listeners == null ) listeners = new ObjVector();
      if ( ! listeners.contains( sl ) ) listeners.addElement( sl );
   }
   // END GUI

   // Functions to facilitate simultaneous printing and logging of output
   public void println() {
      pout.println();
      if ( html ) pout.println( HTML.BR ); }
   public void println( Object str ) {
      pout.println( str.toString() );
      if ( html ) pout.println( HTML.BR ); }
   public void print( char str ) {
      pout.print( str ); } 
   public void print( Object str ) {
      pout.print( str.toString() ); } 
   public void printf( String fmt, double val ) {
      if ( Double.isInfinite( val ) ) val = -Double.MAX_VALUE;
      Format.print( pout, fmt, val ); 
      if ( html && fmt.endsWith( "\n" ) ) pout.println( HTML.BR ); }
   public void printf( String fmt, int val ) {
      Format.print( pout, fmt, val ); 
      if ( html && fmt.endsWith( "\n" ) ) pout.println( HTML.BR ); }
   public void printf( String fmt, String val ) {
      Format.print( pout, fmt, val ); 
      if ( html && fmt.endsWith( "\n" ) ) pout.println( HTML.BR ); }
   protected void DEBUG( Object str ) { if ( debug > 0 ) System.err.println( str ); }

   public void SetupArgs( ArgProcessor argProc ) {
      argProc.AddArg( "Primary Sampling Arguments" );
      argProc.AddArg( "f", "<fname>", "", "input sequence filename or URL" );
      argProc.AddArg( "iseed", "<int>", "10", "random seed" );
      argProc.AddArg( "pow", "<frange:0.01:100>", "1", "power (higher->faster convergence)" );
      argProc.AddArg( "nosamp", null, "false", "deterministic (no sampling->faster)" );
      argProc.AddArg( "i1", "<irange:1:10000>", "500", "number of outer iterations" );
      argProc.AddArg( "i2", "<irange:1:1000>", "50", "number of inner iterations" );
      argProc.AddArg( "nbad", "<irange:0:200>", "10", "number of bad inner iterations allowed" );
      argProc.AddArg( "si", "<frange:0:1>", "0.1", "fraction of small iterations" );
      argProc.AddArg( "sis", "<irange:1:1000>", "1", "number of sequences to run per small iteration" );
      argProc.AddArg( "ni1", "<irange:0:1000>", "20", "number of near-optimum sampling outer iterations" );
      argProc.AddArg( "ni2", "<irange:1:100>", "5", "number of near-optimum sampling inner iterations" );

      argProc.AddArg( "Sequence Pre-Processing Options" );
      argProc.AddArg( "rd", null, "false", "remove duplicate sequences" );
      argProc.AddArg( "cull", "<float>", "0", "randomly cull this fraction/number of sequences from the dataset" );
      argProc.AddArg( "cull-seq", "<string>", "", "cull this named sequence (and its interactions) from the dataset" );
      argProc.AddArg( "wc", null, "false", "search both W and C strands" );
      argProc.AddArg( "worc", null, "false", "search either W or C strands" );
      argProc.AddArg( "prot", null, "false", "convert nucleotide seqs to protein seqs" );

      argProc.AddArg( "Model and Prior Parameters" );
      argProc.AddArg( "W", "<int>", "6", "motif width" );
      argProc.AddArg( "min", "<irange:0:100>", "1", "min. number of motif sites per sequence" );
      argProc.AddArg( "max", "<irange:1:100>", "1", "max. number of motif sites per sequence" );
      argProc.AddArg( "N", "<irange:0:20>", "0", "number of masked repeats to run" );
      argProc.AddArg( "T", "<float>", "0", "negative log of p-value threshold for final motif search" );
      argProc.AddArg( "bgo", "<irange:0:3>", "0", "order of background model" );
      argProc.AddArg( "bgf", "<fname>", "", "file or URL of sequences to generate background model" );
      argProc.AddArg( "fg", "<pseudo|avgscore|dirichlet>", "pseudo", "type of foreground prior model" );
      argProc.AddArg( "P", "<frange:0:1>", "0.001", "pseudocount fraction" );
      argProc.AddArg( "mat", "<fname>", "matrix/BLOSUM62", "scoring matrix filename or URL" );
      argProc.AddArg( "dirich", "<fname>", "dirichlet/hydro-cons.3comp", "dirichlet mixture filename or URL" );
      argProc.AddArg( "p", "<string>", "", "\":\"-separated list of site-based prior class plugins" );

      argProc.AddArg( "Status, Reporting, Misc. Parameters" );
      argProc.AddArg( "rank", "<irange:0:100>", "0", "perform ranking significance tests" );
      argProc.AddArg( "pp", "<string>", "", "class name(s) of sampler post-processor to run before exit (separated by \":\")" );
      argProc.AddArg( "v", null, "false", "be verbose" );
      argProc.AddArg( "vv", null, "false", "be very verbose" );
      argProc.AddArg( "log", "<fname>", "", "save output to given log filename" );
      argProc.AddArg( "tmp", "<fname>", "", "create a temporary file to say \"I'm running!\" (good for clustering)" );
      argProc.AddArg( "S", "<fname>", "", "save final sampler state to file" );
      argProc.AddArg( "save", "<fname>", "", "save alignments to FASTA file" );
      argProc.AddArg( "gff", "<fname>", "", "save alignments to GFF format file" );
      argProc.AddArg( "quiet", null, "false", "be vewwy vewwy quiet" );
      argProc.AddArg( "ansi", null, "false", "enable colored ANSI terminal output" );
      argProc.AddArg( "thread", null, "false", "run in a background thread?" );
      argProc.AddArg( "html", null, "false", "enable colored HTML output" );
      argProc.AddArg( "logo", "<none|gui|text>", "none", "motif logo display mode" );
      argProc.AddArg( "test", "<irange:0:5>", "0", "testing level" );
      argProc.AddArg( "debug", "<irange:0:5>", "0", "debugging level" );
      argProc.AddArg( "h", null, "false", "display this help" );
      argProc.AddArg( "?", null, "false", "display this help" );

      argProc.AddArg( "Simulation Parameters" );
      argProc.AddArg( "insert", "<string>", "", "insert motif(s) into each sequence; <[nper:[gapfreq:]]mot>" );
      argProc.AddArg( "fake", "<string>", "", "generate fake sequences; <J,NS,len|[order:]file>" );
      argProc.AddArg( "saveseqs", "<fname>", "", "save preprocessed seqs to file and exit" );
      argProc.AddArg( "shuffle", null, "false", "shuffle the sequences" );

      // Option flag, type of arg, default, description
      argProc.ModifyDefaultArg( "iseed", String.valueOf( System.currentTimeMillis() ) );
      // START GUI
      argProc.AddArg( "GUI Options" );
      argProc.AddArg( "display", null, "false", "display motif alignment window" );
      argProc.AddArg( "gui", null, "false", "run in Java-powered terminal" );
      argProc.AddArg( "guiopt", null, "false", "run options dialog" );
      // END GUI
   }

   public void SetArgs( ArgProcessor proc ) {
      file = proc.getArg( "f" );
      save = proc.getArg( "save" );
      W = proc.getIntArg( "W" );
      cull = proc.getFloatArg( "cull" );
      seed = proc.getLongArg( "iseed" );
      if ( seed == 0 ) seed = System.currentTimeMillis();
      power = proc.getFloatArg( "pow" );
      noSamp = proc.getBooleanArg( "nosamp" );
      maxIter1 = proc.getIntArg( "i1" );
      maxIter2 = proc.getIntArg( "i2" );
      nbad = proc.getIntArg( "nbad" );
      nearOptimumIter1 = proc.getIntArg( "ni1" );
      nearOptimumIter2 = proc.getIntArg( "ni2" );
      smallIter = proc.getFloatArg( "si" );
      smallIterSize = proc.getIntArg( "sis" );
      bgorder = proc.getIntArg( "bgo" );
      bgfile = proc.getArg( "bgf" );
      fgType = proc.getArg( "fg" );
      matrix = proc.getArg( "mat" );
      dirichName = proc.getArg( "dirich" );
      runMultiples = proc.getIntArg( "N" );
      finalThresh = proc.getFloatArg( "T" );

      minPerSeq = proc.getIntArg( "min" );
      maxPerSeq = proc.getIntArg( "max" );
      pseudo = proc.getFloatArg( "P" ); 
      doRankTest = proc.getIntArg( "rank" );
      sitesPriorsString = proc.getArg( "p" );
      postProcessors = proc.getArg( "pp" );
      verbose = proc.getBooleanArg( "v" );
      vverbose = proc.getBooleanArg( "vv" );
      quiet = proc.getBooleanArg( "quiet" );
      if ( quiet ) { verbose = vverbose = html = ansi = false; debug = 0; }
      threaded = proc.getBooleanArg( "thread" );
      // START GUI
      display = proc.getBooleanArg( "display" );
      gui = proc.getBooleanArg( "gui" );
      askGuiOpt = proc.getBooleanArg( "guiopt" );
      // END GUI
      logos = proc.getArg( "logo" );
      ansi = proc.getBooleanArg( "ansi" );
      html = proc.getBooleanArg( "html" );
      testing = proc.getIntArg( "test" );
      debug = proc.getIntArg( "debug" ); 
   }

   public static Sampler GetSamplerByClassName( String className, Object args ) {
      if ( ! className.startsWith( "djr.motif." ) ) 
	 className = "djr.motif.sampler." + className;
      Object theArgs[] = ( args instanceof String ) ? MyUtils.Tokenize( (String) args, " " ) : 
	 (Object[]) args;
      if ( className.indexOf( "Gapped" ) < 0 ) {
	 boolean gapped = false;
	 for ( int i = 1; i < theArgs.length; i ++ ) {
	    if ( "-gp".equals( theArgs[ i ] ) || "-ge".equals( theArgs[ i ] ) ) { 
	       gapped = true; break; }
	 }
	 if ( gapped ) className += "Gapped";
      }
      Sampler samp = (Sampler) ReflectUtils.CallConstructor( className, 
							     new Object[] { MyUtils.Join( (String[]) theArgs, " " ) } );
      return samp;
   }

   public static Sampler RestoreSamplerFromSavedStateFile( String fname ) {
      Sampler sampler = null;
      if ( fname.toUpperCase().endsWith( ".GZ" ) ) fname = fname.substring( 0, fname.lastIndexOf( '.' ) );
      if ( ! MyUtils.IsNullString( fname ) ) {
	 try { 
	    Map state = (Map) MyUtils.ReadObject( fname + ".gz" );
	    sampler = GetSamplerByClassName( (String) state.get( "classname" ), 
					     (String) state.get( "args" ) );
	    if ( sampler != null ) {
	       sampler.pout = System.out;
	       sampler.RestoreStateFromFile( fname );
	    }
	 } catch( Exception ee ) { sampler = null; }
      }
      return sampler;
   }

   public static void main( String args[] ) {
      String className = args[ 0 ];
      if ( ! className.startsWith( "-" ) ) {
	 if ( ! className.startsWith( "djr.motif.sampler." ) ) 
	    className = "djr.motif.sampler." + className;
	 if ( className.indexOf( "Gapped" ) < 0 ) {
	    boolean gapped = false;
	    for ( int i = 1; i < args.length; i ++ ) {
	       if ( "-gp".equals( args[ i ] ) || "-ge".equals( args[ i ] ) ) { 
		  gapped = true; break; }
	    }
	    if ( gapped ) className += "Gapped";
	 }
	 if ( ReflectUtils.ExecuteStaticMethod( className, "public static void " + 
						className + ".main", 
						new Object[] { args } ) ) return;
      } else if ( "-restore".equals( args[ 0 ] ) ) {
	 // This option reads in a serialized sampler state object from a file and 
	 // re-creates the sampler.
	 Sampler samp = RestoreSamplerFromSavedStateFile( args[ 1 ] );
	 samp.Run();
      } else if ( "-printresults".equals( args[ 0 ] ) ) {
	 // This option reads in a serialized sampler object from a file
	 // and prints out the final output from it.
	 if ( ! MyUtils.IsNullString( args[ 1 ] ) ) {
	    Sampler sampler = null;
	    try {
	       if ( ! args[ 1 ].toUpperCase().endsWith( ".GZ" ) ) args[ 1 ] += ".gz";
	       Map info = (Map) MyUtils.ReadObject( args[ 1 ] );
	       sampler = GetSamplerByClassName( (String) info.get( "classname" ), 
						(String) info.get( "args" ) );
	       if ( sampler != null ) {
		  sampler.pout = System.out;
		  sampler.RestoreStateFromFile( args[ 1 ] );
		  sampler.PrintFinalResults( true, true );
	       }
	    } catch( Exception ee ) { sampler = null; }
	 }
      }
   }
}

