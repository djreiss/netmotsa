package djr.motif.sampler;

import java.io.*;
import java.util.*;
import corejava.*;
import djr.util.*;
import djr.util.bio.*;
import djr.util.array.*;

/**
 * Class <code>SamplerOutput</code>
 *
 * @author <a href="mailto:reiss@uw.edu">David Reiss</a>
 * @version 1.9978 (Fri Nov 07 05:56:26 PST 2003)
 */
public class SamplerOutput {
   final static String GAP = "-";

   public static void PrintLogos( Sampler samp ) {
      boolean ansi = samp.ansi, html = samp.html;
      short J = samp.J;
      for ( int mot = 0; mot < samp.nMotifs; mot ++ ) {
	 double logo[][] = samp.GetMotifModel().ComputeLogoScores( mot );
	 double nsteps = J == 4 ? 10.0 : 20.0;
	 double max = -Double.MAX_VALUE;
	 int WW = samp.W;
	 for ( int i = 0; i < WW; i ++ ) max = DoubleUtils.Max( max, DoubleUtils.Sum( logo[ i ] ) );
	 if ( max <= 0.01 ) return;
	 double step = max / nsteps;
	 int temp[][] = IntUtils.New( WW, (int) nsteps ); IntUtils.Set( temp, -1 );
	 for ( int ii = 0; ii < WW; ii ++ ) {
	    double li[] = logo[ ii ]; int ti[] = temp[ ii ];
	    int indexx[] = DoubleUtils.Indexx( li, null );
	    for ( int j = 0, ind = 0; j < J; j ++ ) {
	       int indx = indexx[ j ];
	       for ( double val = step; val <= li[ indx ]; val += step ) ti[ ind ++ ] = indx;
	    }
	 }

	 if ( html ) samp.println();
	 samp.println( "Sequence Logo (motif " + samp.modelIDs[ mot ] + "):" );
	 int ind = (int) nsteps - 1;
	 for ( double i = max; i > 0.0001; i -= step, ind -- ) {
	    if ( html ) samp.print( HTML.cellStart() );
	    samp.printf( "%5.2f ", i );
	    if ( html ) samp.print( HTML.cellEnd() );
	    for ( int ii = 0; ii < WW; ii ++ ) {
	       int ti[] = temp[ ii ];
	       if ( ti[ ind ] < 0 ) {
		  if ( html ) samp.print( HTML.cellStart() );
		  if ( ! html ) samp.print( "    " );
		  if ( html ) samp.print( HTML.cellEnd() );
	       } else {
		  char ntide = 0; 
		  if ( J == 4 ) ntide = (char) Sequence.SEQ_TABLE[ ti[ ind ] ];
		  else if ( J != 4 ) ntide = (char) Sequence.PROT_TABLE[ ti[ ind ] ];
		  if ( ansi ) samp.print( ANSI.setAttr( Sequence.getANSIColorBG( ntide, J ) ) );
		  else if ( html ) {
		     String colorbg = Sequence.getHTMLColorBG( ntide, J );
		     samp.print( HTML.cellStart( colorbg ) + HTML.BOLD );
		  }
		  samp.printf( "%3s ", ntide+"" );
		  if ( ansi ) samp.print( ANSI.setAttr( ANSI.RESET ) );
		  else if ( html ) samp.print( HTML.cellEnd() );
	       }
	    }
	    if ( html ) samp.print( HTML.NEW_ROW ); else samp.println(); 
	 }
	 if ( html ) samp.print( HTML.cellStart() );
	 samp.print( " BITS " ); 
	 if ( html ) samp.print( HTML.cellEnd() );
	 for ( int ii = 1; ii <= WW; ii ++ ) {
	    if ( html ) samp.print( HTML.cellStart() );
	    samp.printf( "%3d ", ii ); 
	    if ( html ) samp.print( HTML.cellEnd() );
	 }
	 if ( ! html ) samp.println(); else samp.print( HTML.reset() );
	 if ( ! html ) samp.println(); 
      }
   }

   public static void PrintSummary( Sampler samp, Map state, boolean printAlignments,
				    boolean printGarbage ) {
      if ( samp.html ) samp.print( HTML.HR );
      if ( printGarbage ) {
	 if ( ! samp.ansi ) samp.println( "ITER = " + samp.iter1 );
	 else if ( ! samp.html ) samp.println();
	 samp.println( new java.util.Date( System.currentTimeMillis() ) + "\n" );
      }
      int[][][] alignment = null;
      double pvalues[][] = null;
      if ( samp.vverbose || printAlignments ) {
	 alignment = samp.GetAlignments();
	 pvalues = samp.GetMotifModel().ComputePolyScores( alignment );
      }
      if ( printAlignments ) PrintAlignment( samp, alignment, pvalues, true );
      samp.printf( "MAP SCORE = %.3f\n", samp.ComputeMAPScore( state ) );
      int mi[] = samp.GetSites().GetMotifCounts();
      samp.printf( "SITE COUNTS (sum=%d): ", IntUtils.Sum( mi ) );
      for ( int mmm = 0; mmm <= samp.nMotifs; mmm ++ ) samp.printf( "%d ", mi[ mmm ] );
      samp.println();
      samp.print( "WIDTHS: " );
      for ( int mmm = 0; mmm < samp.nMotifs; mmm ++ ) samp.printf( "%5d ", samp.W );
      samp.println();
      if ( samp.vverbose || printAlignments ) {
	 samp.print( "SP SCORES: " );
	 double SP[] = samp.GetMotifModel().ComputeSPScores( alignment );
	 for ( int mot = 0; mot < samp.nMotifs; mot ++ ) samp.printf( "%.3f ", SP[ mot ] );
	 samp.println();
	 samp.print( "E-VALUES: " );
	 double EV[] = samp.GetMotifModel().ComputeEValues( alignment );
	 for ( int mot = 0; mot < samp.nMotifs; mot ++ ) samp.printf( "%.3f ", EV[ mot ] );
	 samp.println();
      }
      samp.println();
      samp.GetMotifModel().Print( samp, samp.ansi, samp.html, false, samp.modelIDs );
      if ( ! MyUtils.IsNullString( samp.logos ) ) {
	 if ( "text".equals( samp.logos ) ) PrintLogos( samp );
	 else if ( "gui".equals( samp.logos ) && 
		   ( samp.vverbose || samp.iter1 >= samp.maxIter1 ) ) samp.DisplayLogos();
      }
      if ( samp.html ) samp.print( HTML.name( "iter" + samp.iter1 ) + "\n" );
      //if ( printAlignments ) SamplerStats.ComputePolyScores( samp, alignment );
   }

   public static void PrintGarbage( Sampler samp, Map state ) {
      samp.printf( "\nMAP SCORE = %.3f\n", samp.ComputeMAPScore( state ) );
      int mi[] = samp.GetSites().GetMotifCounts();
      samp.printf( "SITE COUNTS (sum=%d): ", IntUtils.Sum( mi ) );
      for ( int mmm = 0; mmm <= samp.nMotifs; mmm ++ ) samp.printf( "%d ", mi[ mmm ] );
      samp.println();
      samp.print( "WIDTHS: " );
      for ( int mmm = 0; mmm < samp.nMotifs; mmm ++ ) samp.printf( "%5d ", samp.W/*[ mmm ]*/ );
      samp.println();
   }

   public static void PrintAlignment( Sampler samp, int aa[][][], double pvalues[][], 
				      boolean printAll ) {
      if ( aa == null ) return;
      final int maxw = 20, maxhw = 10;
      boolean saving = ! MyUtils.IsNullString( samp.save );
      boolean html = samp.html, ansi = samp.ansi;
      boolean searchWorC = samp.GetArgProcessor().getBooleanArg( "worc" );
      boolean searchBothWC = samp.GetArgProcessor().getBooleanArg( "wc" );
      Sequence S[] = samp.S;
      short J = samp.J;
      //samp.println( "\nBEST ALIGNMENTS: " );
      if ( html ) samp.print( HTML.BR + HTML.tableStart() );
      String savef = IntUtils.SPrintf( samp.save + "_%05d", samp.iter1 );
      //PrintStream oldpout = samp.pout;
      //samp.pout = new PrintStream( new BufferedOutputStream( oldpout ) );
      if ( html ) samp.println( "<td align=center>Seq</td><td align=center>Mot</td>" +
				"<td align=center>Name</td><td align=center>Strand</td>" +
				"<td>-log(P-<br>value)" +
				"</td><td>Loc</td><td align=right>Alignment</td></tr>" );
      for ( int mot = 0; mot < samp.nMotifs; mot ++ ) {
	 PrintStream sout = null;
	 if ( saving ) {
	    try { sout = new PrintStream( new FileOutputStream( savef + IntUtils.SPrintf( ".%03d", mot+1 ) ) ); }
	    catch( Exception e ) { System.err.println( "Cannot open save file " + samp.save + 
						       (mot+1) + " for output" ); }
	 }
	 if ( aa[ mot ] == null ) continue;
	 for ( int i = 0; i < aa[ mot ].length; i ++ ) {
	    int seq = aa[ mot ][ i ][ 0 ];
	    if ( seq < 0 || S[ seq ] == null ) continue;
	    int locind = 1, loc = aa[ mot ][ i ][ locind ];
	    while ( loc < 0 && locind ++ < aa[ mot ][ i ].length - 1 ) 
	       loc = aa[ mot ][ i ][ locind ];
	    if ( loc < 0 ) continue;
	    boolean inRevComp = ( searchWorC && loc > samp.len[ seq ] / 2 );
	    boolean inRevComp2 = ( searchBothWC && seq % 2 == 1 );
	    char[] sseq = S[ seq ].GetSequence().toCharArray();
	    String header = S[ seq ].GetHeader();
	    String wc = inRevComp || inRevComp2 ? "W" : "C";
	    String temp = maxhw < header.length() ? header.substring( 0, maxhw-1 ) : header;
	    double pp = pvalues != null && pvalues[ mot ] != null ? pvalues[ mot ][ i ] : 
	       samp.GetMotifModel().ComputeMatchScore( aa, i, mot );
	    pp = DoubleUtils.Log10( pp == 0.0 ? 1.0 / Double.MAX_VALUE : pp );
	    if ( ! printAll && pvalues != null && pvalues[ mot ] != null ) continue;
	    if ( sout != null ) sout.println( "> MOTIF = " + (mot+1) + "; " + S[ seq ].GetHeader() );
	    if ( html ) { samp.print( HTML.cellStartA( "center" ) );
	    samp.printf( "%3d", (seq+1) ); }
	    else samp.printf( "%3d ", (seq+1) );
	    if ( html ) { samp.print( HTML.cellEnd() + HTML.cellStartA( "center" ) );
	    samp.print( samp.modelIDs[ mot ] ); }
	    else samp.print( samp.modelIDs[ mot ] + " " );
	    if ( html ) { samp.print( HTML.cellEnd() + HTML.cellStartA( "center" ) );
	    samp.printf( "%" + maxhw + "s", temp ); }
	    else samp.printf( "\"%10s\" ", temp );
	    if ( html ) samp.print( HTML.cellEnd() + HTML.cellStartA( "center" ) + wc );
	    else samp.print( wc + " " );
	    if ( html ) { samp.print( HTML.cellEnd() + HTML.cellStartA( "center" ) );
	    samp.printf( "%5.1f", -pp ); }
	    else samp.printf( "%5.1f ", -pp );
	    if ( html ) { samp.print( HTML.cellEnd() + HTML.cellStartA( "center" ) );
	    samp.printf( "%4d", ( inRevComp ? (loc-samp.len[ seq ]/2+1) : (loc+1) ) ); }
	    else samp.printf( "%4d ", ( inRevComp ? (loc-samp.len[ seq ]/2+1) : (loc+1) ) );
	    if ( html ) samp.print( HTML.cellEnd() + HTML.cellStartA( "right" ) );
	    for ( int j = loc - maxw; j < loc; j ++ ) {
	       if ( ( inRevComp && j < samp.len[seq]/2 ) || ( ! inRevComp && j < 0 ) ) samp.print( " " );
	       else { 
		  if ( ansi ) samp.print( Sequence.setANSIColorFG( sseq[ j ], J ) );
		  else if ( html ) samp.print( Sequence.setHTMLColorFG( sseq[ j ], J ) );
		  else samp.print( sseq[ j ] );
	       }
	       if ( sout != null ) { if ( j < 0 ) sout.print( GAP ); else sout.print( sseq[ j ] ); }
	    }
	    samp.print( " " );
	    if ( html ) samp.print( HTML.cellEnd() );
	    int lastind = -1;
	    for ( int j = 1; j < aa[ mot ][ i ].length; j ++ ) {
	       int jj = aa[ mot ][ i ][ j ];
	       if ( ( ! inRevComp && jj < 0 ) || ( inRevComp && jj < samp.len[seq]/2 ) ) {
		  if ( html ) samp.print( HTML.cellStartA( "center" ) );
		  samp.print( GAP );
		  if ( html ) samp.print( HTML.cellEnd() );
	       } else {
		  if ( ansi ) samp.print( Sequence.setANSIColorBG( sseq[ jj ], J ) );
		  else if ( html ) samp.print( Sequence.setHTMLColorBG( sseq[ jj ], J ) );
		  else samp.print( sseq[ jj ] );
		  if ( jj > 0 ) lastind = jj;
	       }
	       if ( sout != null ) { 
		  if ( ( ! inRevComp && jj < 0 ) || ( inRevComp && jj < samp.len[seq]/2 ) ) sout.print( GAP );
		  else sout.print( sseq[ jj ] );
	       }
	    }
	    if ( ! html ) samp.print( " " );
	    int llen = samp.len[ seq ];
	    int start = lastind + 1;
	    if ( html ) { samp.print( HTML.cellEnd() ); samp.print( HTML.cellStartA( "left" ) ); }
	    for ( int j = start; j < start + maxw; j ++ ) {
	       if ( ( ! inRevComp && j < llen/2 ) || ( ( ! searchWorC || inRevComp ) && j < llen ) ) {
		  if ( ansi ) samp.print( Sequence.setANSIColorFG( sseq[ j ], J ) );
		  else if ( html ) samp.print( Sequence.setHTMLColorFG( sseq[ j ], J ) );
		  else samp.print( sseq[ j ] );
	       } else samp.print( " " );
	       if ( sout != null ) { 
		  if ( ( ! inRevComp && j < llen/2 ) || ( ( ! searchWorC || inRevComp ) && j < llen ) ) sout.print( sseq[ j ] ); 
		  else sout.print( GAP ); 
	       }
	    }
	    if ( html ) samp.print( HTML.cellEnd() + HTML.NEW_ROW );
	    else samp.println();
	    if ( sout != null ) sout.println();
	 }
	 if ( html ) {
	    samp.print( HTML.cellEnd() + HTML.NEW_ROW + HTML.cellStartA( "left" ) + 
			HTML.SPACE + HTML.cellEnd() + HTML.NEW_ROW );
	 } else samp.println();
	 if ( sout != null ) sout.close();
	 if ( saving ) samp.println( "Alignment for motif " + (mot+1) + " saved to file " + 
				     savef + IntUtils.SPrintf( ".%03d", mot+1 ) );
      }
      if ( html ) samp.print( HTML.reset() );
      samp.println();
      //samp.pout.flush();
      //samp.pout = oldpout;
      //System.gc();
   }

   // Format of GFF is a tab-limited of the following info:
   // <seqname> <source> <feature> <start> <end> <score> <strand> <frame> 
   // [attributes] [comments]
   // Feature should be "enhancer" for regulatory element or binding site
   // For source, use name of software
   // Attribute field is optional; can use the actual motif sequence for reference
   public static void SaveAlignmentToGFF( Sampler samp, String gffName, boolean printAll ) {
      int aa[][][] = samp.GetAlignments();
      double pvalues[][] = samp.GetMotifModel().ComputePolyScores( aa );
      if ( aa == null ) return;
      boolean searchWorC = samp.GetArgProcessor().getBooleanArg( "worc" );
      boolean searchBothWC = samp.GetArgProcessor().getBooleanArg( "wc" );
      String cname = ReflectUtils.getClassName( samp );
      Sequence S[] = samp.S;
      short J = samp.J;
      PrintStream gffout = null;
      try { gffout = new PrintStream( MyUtils.OpenOutputFile( gffName, true ) ); }
      catch( Exception e ) { System.err.println( "Cannot open save file " + gffName + 
						 " for output" ); }
      if ( gffout == null ) return;
      for ( int mot = 0; mot < samp.nMotifs; mot ++ ) {
	 if ( aa[ mot ] == null ) continue;
	 for ( int i = 0; i < aa[ mot ].length; i ++ ) {
	    int seq = aa[ mot ][ i ][ 0 ];
	    if ( seq < 0 || S[ seq ] == null ) continue;
	    int locind = 1, loc = aa[ mot ][ i ][ locind ];
	    while ( loc < 0 && locind ++ < aa[ mot ][ i ].length - 1 ) 
	       loc = aa[ mot ][ i ][ locind ];
	    if ( loc < 0 ) continue;
	    boolean inRevComp = ( searchWorC && loc > samp.len[ seq ] / 2 );
	    boolean inRevComp2 = ( searchBothWC && seq % 2 == 1 );
	    char[] sseq = S[ seq ].GetSequence().toCharArray();
	    String wc = inRevComp || inRevComp2 ? "-" : "+";
	    double pp = pvalues != null && pvalues[ mot ] != null ? pvalues[ mot ][ i ] : 
	       samp.GetMotifModel().ComputeMatchScore( aa, i, mot );
	    pp = DoubleUtils.Log10( pp == 0.0 ? 1.0 / Double.MAX_VALUE : pp );
	    if ( ! printAll && pvalues != null && pvalues[ mot ] != null ) continue;

	    // strip out the "(COMPLEMENT)" and "(plus COMPLEMENT)" occurrences
	    String header = S[ seq ].GetHeader();
	    int index = header.indexOf( Sampler.SEQUENCE_HEADER_COMPLEMENT );
	    if ( -1 != index ) header = header.substring( 0, index ).trim();
	    else {
	       index = header.indexOf( Sampler.SEQUENCE_HEADER_PLUS_COMPLEMENT );
	       if ( -1 != index ) header = header.substring( 0, index ).trim();
	    }
	    gffout.print( header + "\t" + cname + "\t" );
	    String feature = samp.modelIDs[ mot ];
	    if ( feature.equals( "" + ( mot + 1 ) ) && J == 4 ) feature = "enhancer";
	    gffout.print( feature + "\t" );
	    int start = inRevComp ? loc-samp.len[ seq ]/2+1 : loc+1;
	    int end = inRevComp ? start - samp.W + 1 : start + samp.W - 1;
	    gffout.print( start+ "\t" + end + "\t" );
	    gffout.print( DoubleUtils.SPrintf( "%5.1f\t", -pp ) );
	    gffout.print( wc + "\t.\t" );

	    for ( int j = 1; j < aa[ mot ][ i ].length; j ++ ) {
	       int jj = aa[ mot ][ i ][ j ];
	       if ( ( ! inRevComp && jj < 0 ) || ( inRevComp && jj < samp.len[seq]/2 ) )
		  gffout.print( GAP );
	       else gffout.print( sseq[ jj ] );
	    }
	    gffout.println();
	 }
      }
      gffout.flush(); gffout.close();
   }
}

