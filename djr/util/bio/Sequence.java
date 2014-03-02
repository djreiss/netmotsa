package djr.util.bio;
import java.util.*;
import java.io.*;

import iubio.bioseq.*;
import iubio.readseq.*;
import djr.util.*;
import djr.util.array.*;
//import djr.util.bio.align.*;
import djr.motif.model.BackgroundModel;

/**
 * <code>Sequence</code> class.
 *
 * @author <a href="mailto:dreiss@systemsbiology.org">David Reiss</a>
 * @version 1.9978 (Fri Nov 07 05:56:26 PST 2003)
 */
public class Sequence implements Serializable {
   public static final short PROTEIN = 1;
   public static final short DNA = 2;
   public static final short RNA = 3;

   public static char NULL_CHARACTER = '-';
   public static String NUCLEOTIDES = "GATC";
   public static String ACIDS = "ACDEFGHIKLMNPQRSTVWY";
   public static char SEQ_TABLE[] = NUCLEOTIDES.toCharArray();
   public static char PROT_TABLE[] = ACIDS.toCharArray();
   public static char PROT_HYDROPHOBIC[] = "ILVCAGMFYWHKT".toCharArray();
   public static char PROT_AROMATIC[] = "FYWH".toCharArray();
   public static char PROT_CHARGED_POS[] = "RHK".toCharArray();
   public static char PROT_CHARGED_NEG[] = "DE".toCharArray();
   public static char PROT_POLAR[] = "YWHKREQDNST".toCharArray();
   public static short REV_SEQ_TABLE[] = new short[ (int) 'Z' + 1 ];
   public static short REV_PROT_TABLE[] = new short[ (int) 'Z' + 1 ];
   public static char VALID_DNA[] = "GATCUXN[]".toCharArray();
   public static String CODON_TABLE[][] = { 
      { "TCA", "S" }, { "TCG", "S" }, { "TCC", "S" }, { "TCT", "S" }, { "TTT", "F" }, 
      { "TTC", "F" }, { "TTA", "L" }, { "TTG", "L" }, { "TAT", "Y" }, { "TAC", "Y" }, 
      { "TAA", "*" }, { "TAG", "*" }, { "TGT", "C" }, { "TGC", "C" }, { "TGA", "*" }, 
      { "TGG", "W" }, { "CTA", "L" }, { "CTG", "L" }, { "CTC", "L" }, { "CTT", "L" }, 
      { "CCA", "P" }, { "CCG", "P" }, { "CCC", "P" }, { "CCT", "P" }, { "CAT", "H" }, 
      { "CAC", "H" }, { "CAA", "Q" }, { "CAG", "Q" }, { "CGA", "R" }, { "CGG", "R" }, 
      { "CGC", "R" }, { "CGT", "R" }, { "ATT", "I" }, { "ATC", "I" }, { "ATA", "I" }, 
      { "ATG", "M" }, { "ACA", "T" }, { "ACG", "T" }, { "ACC", "T" }, { "ACT", "T" }, 
      { "AAT", "N" }, { "AAC", "N" }, { "AAA", "K" }, { "AAG", "K" }, { "AGT", "S" }, 
      { "AGC", "S" }, { "AGA", "R" }, { "AGG", "R" }, { "GTA", "V" }, { "GTG", "V" }, 
      { "GTC", "V" }, { "GTT", "V" }, { "GCA", "A" }, { "GCG", "A" }, { "GCC", "A" }, 
      { "GCT", "A" }, { "GAT", "D" }, { "GAC", "D" }, { "GAA", "E" }, { "GAG", "E" }, 
      { "GGA", "G" }, { "GGG", "G" }, { "GGC", "G" }, { "GGT", "G" } };
   public static Map CODONS = new java.util.HashMap();

   public static Map HYDROPHOBICITIES = new java.util.HashMap();

   public static String PRINT_FORMAT = "ANSI"; // "PLAIN", "ANSI" or "HTML"
   public static String FG_N_HTML_COLORS[] = { "#E6E600", "#E60A0A", "#0F820F", "#00DCDC" };
   public static String BG_N_HTML_COLORS[] = FG_N_HTML_COLORS;

   public static String FG_N_ANSI_COLORS[] = { ANSI.FGYELLOW, ANSI.FGRED, ANSI.FGGREEN, 
					       ANSI.FGCYAN };
   public static String BG_N_ANSI_COLORS[] = { ANSI.BGYELLOW, ANSI.BGRED, ANSI.BGGREEN, 
					       ANSI.BGCYAN };

   // Acid colors taken from http://info.bio.cmu.edu/Courses/BiochemMols/RasFrames/SHAPELY.HTM
   public static String FG_AA_ANSI_COLORS[] = { ANSI.FGWHITE + ";" + ANSI.UNDERSCORE, // grey
						ANSI.FGYELLOW, ANSI.FGRED, ANSI.FGRED,
						ANSI.FGBLUE, ANSI.FGWHITE, 
						ANSI.FGBLUE + ";" + ANSI.UNDERSCORE,
						ANSI.FGGREEN, ANSI.FGBLUE, ANSI.FGGREEN, 
						ANSI.FGYELLOW, ANSI.FGCYAN, ANSI.FGMAGENTA,
						ANSI.FGCYAN, ANSI.FGBLUE, ANSI.FGYELLOW, 
						ANSI.FGMAGENTA + ";" + ANSI.UNDERSCORE,
						ANSI.FGGREEN, ANSI.FGMAGENTA, 
						ANSI.FGBLUE + ";" + ANSI.UNDERSCORE };
   public static String BG_AA_ANSI_COLORS[] = { ANSI.BGWHITE + ";" + ANSI.UNDERSCORE, // grey
						ANSI.BGYELLOW, ANSI.BGRED, ANSI.BGRED,
						ANSI.BGBLUE, ANSI.BGWHITE, 
						ANSI.BGBLUE + ";" + ANSI.UNDERSCORE,
						ANSI.BGGREEN, ANSI.BGBLUE, ANSI.BGGREEN, 
						ANSI.BGYELLOW, ANSI.BGCYAN, ANSI.BGMAGENTA,
						ANSI.BGCYAN, ANSI.BGBLUE, ANSI.BGYELLOW,
						ANSI.BGMAGENTA + ";" + ANSI.UNDERSCORE,
						ANSI.BGGREEN, ANSI.BGMAGENTA, 
						ANSI.BGBLUE + ";" + ANSI.UNDERSCORE };
   public static String FG_AA_HTML_COLORS[] = { "#C8C8C8", "#E6E600", "#E60A0A", "#E60A0A",
						"#3232AA", "#ABABAB", "#8282D2", "#0F820F", 
						"#145AFF", "#0F820F", "#E6E600", "#00DCDC", 
						"#DC9682", "#00DCDC", "#145AFF", "#FA9600", 
						"#FA9600", "#0F820F", "#B45AB4", "#3232AA" };
   public static String BG_AA_HTML_COLORS[] = FG_AA_HTML_COLORS;

   public static final String getHTMLColorFG( short res, short J ) { 
      return J == 4 ? FG_N_HTML_COLORS[ res ] : FG_AA_HTML_COLORS[ res ]; }
   public static final String getHTMLColorFG( char res, short J ) { 
      return J == 4 ? FG_N_HTML_COLORS[ REV_SEQ_TABLE[ res ] ] : 
	 FG_AA_HTML_COLORS[ REV_PROT_TABLE[ res ] ]; }
   public static final String getHTMLColorBG( short res, short J ) { 
      return J == 4 ? BG_N_HTML_COLORS[ res ] : BG_AA_HTML_COLORS[ res ]; }
   public static final String getHTMLColorBG( char res, short J ) { 
      return J == 4 ? BG_N_HTML_COLORS[ REV_SEQ_TABLE[ res ] ] : 
	 BG_AA_HTML_COLORS[ REV_PROT_TABLE[ res ] ]; }

   public static final String getANSIColorFG( short res, short J ) { 
      return J == 4 ? FG_N_ANSI_COLORS[ res ] : FG_AA_ANSI_COLORS[ res ]; }
   public static final String getANSIColorFG( char res, short J ) { 
      return J == 4 ? FG_N_ANSI_COLORS[ REV_SEQ_TABLE[ res ] ] : 
	 FG_AA_ANSI_COLORS[ REV_PROT_TABLE[ res ] ]; }
   public static final String getANSIColorBG( short res, short J ) { 
      return ( J == 4 ? BG_N_ANSI_COLORS[ res ] : BG_AA_ANSI_COLORS[ res ] ) + ";" + ANSI.FGBLACK; }
   public static final String getANSIColorBG( char res, short J ) { 
      return ( J == 4 ? BG_N_ANSI_COLORS[ REV_SEQ_TABLE[ res ] ] : 
	       BG_AA_ANSI_COLORS[ REV_PROT_TABLE[ res ] ] ) + ";" + ANSI.FGBLACK; }

   public static final String setHTMLColorFG( short res, short J ) {
      return //HTML.cellStartA( "center" ) + 
	 HTML.colorStart( Sequence.getHTMLColorFG( res, J ) ) + 
	 ( J == 4 ? SEQ_TABLE[ res ] : PROT_TABLE[ res ] ) + 
	 HTML.colorStop() /*+ HTML.cellEnd()*/; }
   public static final String setHTMLColorFG( char res, short J ) {
      return //HTML.cellStartA( "center" ) + 
	 HTML.colorStart( Sequence.getHTMLColorFG( res, J ) ) + res + 
	 HTML.colorStop() /*+ HTML.cellEnd()*/; }
   public static final String setHTMLColorBG( short res, short J ) {
      return HTML.cellStart( getHTMLColorBG( res, J ), "center" ) + HTML.BOLD +
	 ( J == 4 ? SEQ_TABLE[ res ] : PROT_TABLE[ res ] ) + 
	 HTML.cellEnd(); }
   public static final String setHTMLColorBG( char res, short J ) {
      return HTML.cellStart( getHTMLColorBG( res, J ), "center" ) + HTML.BOLD + res +
	 HTML.cellEnd(); }

   public static final String setANSIColorFG( short res, short J ) {
      return ANSI.setAttr( getANSIColorFG( res, J ) ) + 
	 ( J == 4 ? SEQ_TABLE[ res ] : PROT_TABLE[ res ] ) + 
	 ANSI.setAttr( ANSI.RESET ); }
   public static final String setANSIColorFG( char res, short J ) {
      return ANSI.setAttr( getANSIColorFG( res, J ) ) + res + 
	 ANSI.setAttr( ANSI.RESET ); }
   public static final String setANSIColorBG( short res, short J ) {
      return ANSI.setAttr( getANSIColorBG( res, J ) ) +
	 ( J == 4 ? SEQ_TABLE[ res ] : PROT_TABLE[ res ] ) + 
	 ANSI.setAttr( ANSI.RESET ); }
   public static final String setANSIColorBG( char res, short J ) {
      return ANSI.setAttr( getANSIColorBG( res, J ) ) + res + 
	 ANSI.setAttr( ANSI.RESET ); }

   static {
      for ( short i = 0; i < 4; i ++ ) REV_SEQ_TABLE[ (int) SEQ_TABLE[ i ] ] = i;
      REV_SEQ_TABLE[ (int) 'U' ] = 2;
      for ( short i = 0; i < 20; i ++ ) REV_PROT_TABLE[ (short) PROT_TABLE[ i ] ] = i;
      for ( int i = 0; i < CODON_TABLE.length; i ++ ) 
	 CODONS.put( CODON_TABLE[ i ][ 0 ], CODON_TABLE[ i ][ 1 ] );
      HYDROPHOBICITIES.put( "Eisenberg", 
			    new double[] { 0.4, -0.6, -4.2, -3.7, 2.2, 0.0, -2.7, 2.7, 
					   -6.0, 1.8, 0.5, -3.8, -1.1, -4.1, -9.2, -2.0, 
					   -1.6, 1.8, 1.0, -0.7 } );
      HYDROPHOBICITIES.put( "Cornette", 
			    new double[] {0.2, 4.1, -3.1, -1.8, 4.4, 0.0, 0.5, 4.8, -3.1, 
					  5.7, 4.2, -0.5, -2.2, -2.8, 1.4, -0.5, -1.9, 4.7,
					  1.0, 3.2 } );
      HYDROPHOBICITIES.put( "Boyko",
			    new double[] { -0.4, -0.5, 6.2, 3.6, -8.8, 0.0, 2.4, -9.5, 6.1,
					   -10.6, -8.5, 1.2, 3.0, 6.4, -1.8, 0.9, 3.8, 
					   -9.7, -2.0, -7.4 } );
      HYDROPHOBICITIES.put( "ParkerHPLC", 
			    new double[] { 2.1, 1.4, 10.0, 7.8, -9.2, 5.7, 2.1, -8.0, 5.7, 
					   -9.2, -4.2, 7.0, 2.0, 6.0, 4.2, 6.5, 5.2, -3.7, 
					   -10.0, -1.9 } );
      HYDROPHOBICITIES.put( "KyteDoolittle", 
			    new double[] { 1.8, 2.5, -3.5, -3.5, 2.8, -0.4, -3.2, 4.5, 
					   -3.9, 3.8, 1.9, -3.5, -1.6, -3.5, -4.5, -0.8, 
					   -0.7, 4.2, -0.9, -1.3 } );
   }

   String seq, header;
   short[] residues;
   short type;
   int length;

   public String GetSequence() { return seq; }
   public String GetHeader() { return header; }
   public void SetHeader( String newH ) { header = newH; }
   public String GetName() { return header; }
   public short GetType() { return type; }
   public String GetTypeName() { return ( type == PROTEIN ? "protein" : ( type == DNA ? "DNA" : "RNA" ) ); }
   public int GetLength() { return length > 0 ? length : 0; }
   public short[] GetResidues() { return residues; }
   public short[] GetData() { return residues; }
   public short GetAlphabetSize() { return GetType() == PROTEIN ? (short) 20 : (short) 4; }

   public boolean IsMasked( int ind ) { return GetResidues()[ ind ] < 0; }
   public boolean IsMaskedInRange( int start, int end ) { 
      short r[] = GetResidues();
      if ( start < 0 ) start = 0;
      if ( end >= GetLength() ) end = GetLength() - 1;
      for ( int i = start; i <= end; i ++ ) if ( r[ i ] < 0 ) return true;
      return false;
   }

   // These are destructive!!!
   public void SetMask( int ind ) { GetResidues()[ ind ] = -1; }
   public void MaskRange( int start, int end, int val ) { 
      ShortUtils.Set( GetResidues(), start, end + 1, (short) -1 ); }

   public Sequence( String s, String h ) {
      seq = header = null;
      residues = null;
      Initialize( s, h );
   }

   public Sequence( String s ) {
      seq = header = null;
      residues = null;
      Initialize( s, null );
   }

   public Sequence( short[] res ) {
      this( res, null ); 
   }

   public Sequence( short[] res, String head ) {
      seq = null;
      this.header = head;
      residues = null;

      //length = len;
      length = res.length;
      residues = ShortUtils.New( length );
      ShortUtils.Copy( residues, res );
      type = DNA;
      for ( int i = 0; i < length; i ++ ) {
	 if ( type == DNA && residues[ i ] > 3 ) { type = PROTEIN; break; }
      }
      if ( header == null ) header = "NO HEADER";
      seq = stringFromResidues( residues, type );
   }

   public Sequence( Sequence in ) {
      this( in.GetSequence(), in.GetHeader() ); 
   }

   protected static String stringFromResidues( short residues[], short type ) {
      int len = residues.length;
      StringBuffer buf = new StringBuffer( len );
      switch( type ) {
      case DNA: case RNA: {
	 for ( int i = 0; i < len; i ++ ) {
	    short ind = residues[ i ];
	    if ( ind < SEQ_TABLE.length ) buf.append( SEQ_TABLE[ ind ] );
	 }
	 break;
      }
      case PROTEIN: {
	 for ( int i = 0; i < len; i ++ ) {
	    short ind = residues[ i ];
	    if ( ind == -1 ) buf.append( NULL_CHARACTER );
	    else if ( ind < PROT_TABLE.length ) buf.append( PROT_TABLE[ ind ] );
	 }	
	 break;
      }
      default: { break; }
      }
      return new String( buf );
   }

   protected void Initialize( String s, String h ) {
      length = 0;
      if ( s != null ) {
	 seq = new String( s );
	 length = seq.length();
	 char[] buf = seq.toCharArray();
	 int len = length;
	 for ( int i = 0; i < len; i ++ ) 
	    if ( buf[ i ] != '-' && ! Character.isLetter( buf[ i ] ) ) length --;
	 residues = ShortUtils.New( length );
	 if ( type != PROTEIN ) type = DNA;
	 for ( int i = 0, j = 0; i < len; i ++ ) {
	    if ( buf[ i ] == '-' || buf[ i ] == 'x' || buf[ i ] == 'X' ) {
	       residues[ j ++ ] = -1;
	       continue;
	    }
	    if ( ! Character.isLetter( buf[ i ] ) ) continue;
	    buf[ i ] = Character.toUpperCase( buf[ i ] );
	    residues[ j ++ ] = (short) buf[ i ];
	    if ( type == DNA && ( buf[ i ] != 'G' && buf[ i ] != 'A' && buf[ i ] != 'T' &&
				  buf[ i ] != 'C' && buf[ i ] != 'U' ) ) {
	       type = PROTEIN;
	    } else if ( type == DNA && buf[ i ] == 'U' ) type = RNA;
	 }
      }
      if ( h != null ) header = new String( h );
      else header = "NO HEADER";

      if ( header.toLowerCase().indexOf( "[saccharomyces cerevisiae]" ) >= 0 ) {
	 header = GetGeneNameFromYeastGenome( header ).toUpperCase();
	 if ( header.endsWith( "P" ) ) header = header.substring( 0, header.length() - 2 );
	 int ind = header.indexOf( 'Y' );
	 if ( ind != 0 ) {
	    try {
	       ObjVector lines = MyUtils.ReadFileLines( "synonyms_all.noa" );
	       for ( int i = 1, sz = lines.size(); i < sz; i ++ ) {
		  String toks[] = MyUtils.Tokenize( (String) lines.get( i ), " " );
		  if ( toks[ 0 ].equals( header ) ) { header = toks[ 2 ]; break; }
	       }
	    } catch ( Exception e ) { };//e.printStackTrace(); }
	 }
	 if ( header.indexOf( '-' ) >= 0 ) header = header.substring( 0, header.indexOf( '-' ) - 1 );
      }

      switch( type ) {
      case DNA: case RNA: {
	 for ( int i = 0; i < length; i ++ ) {
	    short ind = residues[ i ];
	    residues[ i ] = REV_SEQ_TABLE[ ind ];
	 }
	 break;
      }
      case PROTEIN: {
	 for ( int i = 0; i < length; i ++ ) {
	    short ind = residues[ i ];
	    if ( ind == -1 ) continue;
	    residues[ i ] -= 'A';
	    //if ( ind >= 'Z' ) residues[ i ] -= 6; // No acids BJOUXZ
	    /*else*/ if ( ind >= 'X' ) residues[ i ] -= 5;
	    else if ( ind >= 'U' ) residues[ i ] -= 4;
	    else if ( ind >= 'O' ) residues[ i ] -= 3;
	    else if ( ind >= 'J' ) residues[ i ] -= 2;
	    else if ( ind >= 'B' ) residues[ i ] --;
	 }
	 break;
      }
      default: { break; }
      }
      seq = stringFromResidues( residues, type );
   }

   public String GetReverse() {
      return new String( ( new StringBuffer( seq ) ).reverse() );
   }

   public int[] GetResiduesAsInts() { 
      int out[] = IntUtils.New( residues.length );
      for ( int i = 0, s = residues.length; i < s; i ++ ) out[ i ] = (int) residues[ i ];
      return out;
   }

   public String GetComplement() {
      if ( type == PROTEIN ) return new String( ( new StringBuffer( seq ) ).reverse() );
      StringBuffer comp = new StringBuffer( length );
      char[] buf = seq.toCharArray();
      for ( int i = 0; i < length; i ++ ) {
	 char residue = buf[ i ];
	 switch( residue ) {
	 case 'G': comp.append( 'C' ); break;
	 case 'C': comp.append( 'G' ); break;
	 case 'A': comp.append( 'T' ); break;
	 case 'T': comp.append( 'A' ); break;
	 case 'U': comp.append( 'A' ); break;
	 }
      }
      return new String( comp.reverse() );
   }

   public String toString() {
      //String t = GetHeader() + ": " + ( type == DNA ? "DNA" : ( type == PROTEIN ? "PROTEIN" : 
      //				 ( type == RNA ? "RNA" : "Unknown" ) ) );
      //t += " (" + t + "): " + seq + " " + GetLength() + "\n[";
      //for ( int i = 0; i < GetLength(); ++ i ) t += " " + residues[ i ];
      //t += " ]\n";
      //return t;
      return "ANSI".equals( PRINT_FORMAT ) ? ToANSIString() : 
	 ( "HTML".equals( PRINT_FORMAT ) ? ToHTMLString( true ) : GetSequence() );
   }

   public String ToANSIString() {
      char c[] = GetSequence().toCharArray();
      short J = GetAlphabetSize();
      StringBuffer out = new StringBuffer();
      for ( int i = 0, s = c.length; i < s; i ++ )
	 out.append( setANSIColorFG( c[ i ], J ) );
      return out.toString();
   }

   public String ToHTMLString( boolean tableWrapper ) {
      char c[] = GetSequence().toCharArray();
      short J = GetAlphabetSize();
      StringBuffer out = new StringBuffer();
      if ( tableWrapper ) out.append( HTML.tableStart() );
      for ( int i = 0, s = c.length; i < s; i ++ ) {
	 out.append( HTML.cellStartA( "center" ) );
	 out.append( setHTMLColorFG( c[ i ], J ) );
	 out.append( HTML.cellEnd() );
      }
      if ( tableWrapper ) out.append( HTML.reset() );
      return out.toString();
   }

   public boolean Equals( Sequence other ) {
      return GetLength() == other.GetLength() &&
	 ShortUtils.Equals( other.GetResidues(), this.GetResidues() );
   }

   //static Match.Blosum50 blosum = null;
   //static Match.SW sw = null;

   /*
   public double GetMatchScore( Sequence other ) {
      return GetMatchScore( GetSequence(), other.GetSequence() );
   }

   public static double GetMatchScore( String seq1, String seq2 ) {
      if ( blosum == null ) {
	 blosum = new Match.Blosum50();
	 sw = new Match.SW( blosum, 8, seq1, seq2 );
      } else sw.reset( seq1, seq2, false );
      return ( (double) sw.getScore() / (double) ( seq1.length() + seq2.length() ) );
   }
   */

   public boolean IsHomologous( Sequence other ) {
      if ( this.Equals( other ) ) return true;
      Sequence s1 = GetLength() >= other.GetLength() ? this : other;
      Sequence s2 = s1 == this ? other : this;
      String ss1 = s1.GetSequence(), ss2 = s2.GetSequence();
      if ( ss1.indexOf( ss2 ) >= 0 ) return true;
      int len = s2.GetLength(), frac = 2, half = len / frac;
      while( half > 50 ) {
	 if ( ss1.indexOf( ss2.substring( len - half ) ) >= 0 ) return true;
	 for ( int ind = 0; ind < len - half - 6; ind += half )
	    if ( ss1.indexOf( ss2.substring( ind, ind + half ) ) >= 0 ) return true;
	 frac ++; half = len / frac;
      }
      //double score = GetMatchScore( other );
      //return score > 1.0;
      return false;
   }

   public static int GetMaxLength( Sequence seqs[] ) {
      int out = -1;
      for ( int i = 0, ns = seqs.length; i < ns; i ++ ) 
	 if ( seqs[ i ].GetLength() > out ) out = seqs[ i ].GetLength();
      return out;
   }

   public static String GetGeneNameFromYeastGenome( String header ) {
      String s[] = MyUtils.Tokenize( header, " []" );
      int ind = 0;
      while( ind < s.length && ! s[ ind ].toUpperCase().startsWith( "SACCHAROMYCES" ) ) ind ++;
      if ( ind >= s.length ) return s[ 0 ];
      return s[ ind - 1 ].substring( 0, s[ ind - 1 ].length() - 1 );
   }

   public static int GetOrfStartFromYeastGenome( String header ) {
      if ( header.indexOf( '|' ) < 0 ) return -1;
      String s[] = MyUtils.Tokenize( header, "|" );
      return Integer.parseInt( s[ 1 ] );
   }

   public static void RemoveHomologousSequences( Map seqs, boolean verbose ) {
      Iterator it1 = seqs.keySet().iterator();
      ObjVector removed = new ObjVector();
      int i = 1;
      while( it1.hasNext() ) {
	 String p1 = (String) it1.next();
	 if ( removed.contains( p1 ) ) continue;
	 String pp1 = Sequence.GetGeneNameFromYeastGenome( p1 );
	 if ( verbose ) System.err.println( "TESTING " + i + " " + pp1 );
	 Sequence s1 = (Sequence) seqs.get( p1 );
	 Iterator it2 = seqs.keySet().iterator();
	 while( it2.hasNext() ) {
	    String p2 = (String) it2.next();
	    if ( removed.contains( p2 ) ) continue;
	    String pp2 = Sequence.GetGeneNameFromYeastGenome( p2 );
	    Sequence s2 = (Sequence) seqs.get( p2 );
	    if ( p1.equals( p2 ) ) continue;
	    /*double score = s1.GetMatchScore( s2 );
	    if ( score > 1.0 ) {
	       Sequence toRemove = s1.GetLength() < s2.GetLength() ? s1 : s2;
	       String remove = s1.GetLength() < s2.GetLength() ? pp1 : pp2;
	       String keep = s1.GetLength() < s2.GetLength() ? pp2 : pp1;
	       if ( verbose ) 
		  System.err.println( "REMOVING " + remove + " (sim. to " + keep + 
				      DoubleUtils.SPrintf( "; score = %.3f)", score ) );
	       removed.addElement( s1.GetLength() < s2.GetLength() ? p1 : p2 );
	       if ( verbose ) System.err.println( "Removed: " + removed.size() + " of " +
						  seqs.size() );
	       if ( remove == pp2 ) continue;
	    }*/
	 }
      }
      if ( verbose ) System.err.println( "Removed " + removed.size() + " sequences:" );
      for ( i = 0; i < removed.size(); i ++ ) {
	 if ( verbose ) System.err.println( removed.elementAt( i ) );
	 seqs.remove( removed.elementAt( i ) );
      }
   }

   public static Map SequenceArrayToHashtable( Sequence S[] ) {
      Map out = new java.util.HashMap();
      for ( int i = 0, size = S.length; i < size; i ++ ) 
	 out.put( S[ i ].GetHeader().toUpperCase(), S[ i ] );
      return out;
   }

   public static Map ReadSequencesToHashtable( String fname ) throws Exception {
      Sequence S[] = ReadSequences( fname );
      return SequenceArrayToHashtable( S );
   }

   public static Sequence[] ReadSequences( String fname ) throws Exception {
      return ReadSequences( fname, Integer.MAX_VALUE );
   }

   public static Sequence[] ReadSequences( String fname, int max ) throws Exception {
      ObjVector vec = new ObjVector();
      InputStreamReader dis = null;
      try {
	 DataInputStream ddis = new DataInputStream( MyUtils.OpenFile( fname ) );
	 dis = new InputStreamReader( ddis );
      } catch( Exception e ) { e.printStackTrace(); dis = null; }
      if ( dis == null ) return null;
      Readseq rd = new Readseq();
      String seqname = rd.setInputObject( dis );
      int ii = 0;
      if ( rd.isKnownFormat() && rd.readInit() )  {
	 while ( rd.readNext() && ii ++ < max ) {
	    BioseqRecord seqrec = new BioseqRecord( rd.nextSeq() );
	    String title = ( seqrec.getdoc() != null && seqrec.getdoc().getTitle() != null ) ? seqrec.getdoc().getTitle() : seqrec.getID();
	    try { 
	       Sequence s = new Sequence( seqrec.getseq().toString(), title ); 
	       vec.addElement( s );
	    } catch( Throwable t ) { 
	       System.err.println( "Sequence " + title + " cannot be translated. Skipping." ); 
	       t.printStackTrace();
	    }
	 }
      } 
      dis.close();
      MyUtils.DeleteTempFiles( "readseq-", ".tmp" );
      Sequence[] out = new Sequence[ vec.size() ];
      for ( int i = 0, size = vec.size(); i < size; i ++ ) 
	 out[ i ] = (Sequence) vec.elementAt( i );
      return out;
   }

   public static int WriteSequences( String fname, Map seqs ) {
      Sequence S[] = new Sequence[ seqs.size() ];
      Iterator it = seqs.keySet().iterator();
      int i = 0;
      while( it.hasNext() ) S[ i ++ ] = (Sequence) seqs.get( it.next() );
      return WriteSequences( fname, S );
   }

   public static int WriteSequences( String fname, Sequence S[] ) {
      if ( S == null || S.length <= 0 ) return 0;
      try {
	 BioseqWriterIface writer = 
	    BioseqFormats.newWriter( BioseqFormats.formatFromName( "fasta" ) );
	 //writer.setOutput( new FileWriter( fname ) );
	 OutputStream ostream = MyUtils.OpenOutputFile( fname );
	 writer.setOutput( new OutputStreamWriter( ostream ) ); 
	 writer.writeHeader();
	 for ( int i = 0; i < S.length; i ++ ) {
	    Bioseq seq = new Bioseq( S[ i ].GetSequence() );
	    BasicBioseqDoc seqdoc = new BasicBioseqDoc( S[ i ].GetHeader() );
	    seqdoc.addSequenceStats( seq );
	    BioseqRecord seqrec = new BioseqRecord( seq, seqdoc );
	    if ( writer.setSeq( seqrec ) ) writer.writeSeqRecord();
	    else System.err.println( "Failed to write " + seqrec ); // or throw exception
	 }
	 ostream.flush();
	 writer.writeTrailer();
	 writer.close();
      } catch (Exception e) { e.printStackTrace(); }
      return 1;
   }

   public Sequence Shuffle() {
      int l = GetLength();
      short[] flags = ShortUtils.New( l );
      StringBuffer sseq = new StringBuffer( l );
      char[] old = GetSequence().toCharArray();
      int ind = 0;
      while( ind < l ) {
	 int pos = (int)( DoubleUtils.Random() * l );
	 if ( flags[ pos ] != 0 || old[ pos ] == 0 ) continue;
	 sseq.append( old[ pos ] );
	 flags[ pos ] = 1;
      }
      Initialize( new String( sseq ), new String( GetHeader() ) );
      return this;
   }

   public char GetRandomResidue() {
      if ( GetType() == DNA ) return SEQ_TABLE[ (int)( DoubleUtils.Random() * 4 ) ];
      return PROT_TABLE[ (int)( DoubleUtils.Random() * 20 ) ];
   }

   static ObjVector names = new ObjVector();
   static char alph[] = "ABCDEFGHIJKLMNOPQRSTUVWXYZ".toCharArray();

   public static String CreateFakeGeneName() {
      String name = "" + CharUtils.RandChoose( alph ) + CharUtils.RandChoose( alph ) +
	 CharUtils.RandChoose( alph ) + IntUtils.RandChoose( 1, 9 );
      while( names.contains( name ) ) {
	 name = "" + CharUtils.RandChoose( alph ) + CharUtils.RandChoose( alph ) +
	    CharUtils.RandChoose( alph ) + IntUtils.RandChoose( 1, 9 );
      }
      return name;
   }

   public static Sequence FakeSequence( int llen, int JJ ) {
      char[] seq = CharUtils.New( llen );
      for ( int j = 0; j < llen; j ++ ) {
	 if ( JJ == 4 ) seq[ j ] = CharUtils.RandChoose( Sequence.SEQ_TABLE );
	 else seq[ j ] = CharUtils.RandChoose( Sequence.PROT_TABLE );
      }
      return new Sequence( new String( seq ), CreateFakeGeneName() );
   }

   public static Sequence[] CreateFakeSequences( String fake, double pseudo ) {
      Sequence SS[] = null;
      if ( fake.indexOf( ',' ) >= 0 ) {
	 double[] f = DoubleUtils.FTokenize( fake, "," );
	 int JJ = (int) f[ 0 ];
	 int NS = (int) f[ 1 ];
	 int llen = (int) f[ 2 ];
	 SS = new Sequence[ NS ];
	 for ( int i = 0; i < NS; i ++ ) {
	    Sequence tempS = Sequence.FakeSequence( llen, JJ );
	    SS[ i ] = new Sequence( tempS.GetResidues(), Sequence.CreateFakeGeneName() );
	 }
      } else {
	 int order = 1;
	 try {
	    try { 
	       if ( fake.indexOf( ':' ) >= 0 ) {
		  order = Integer.parseInt( fake.substring( 0, fake.indexOf( ':' ) ) );
		  fake = fake.substring( fake.indexOf( ':' ) + 1 );
	       }
	    } catch( Exception e ) { order = 1; }
	    SS = Sequence.ReadSequences( fake );
	 } catch( Exception e ) {
	    System.err.println( "File '" + fake + "' not found for generating fake sequences.\n" );
	    e.printStackTrace();
	    return null;
	 }
	 if ( SS != null ) {
	    BackgroundModel bg = new BackgroundModel( SS, (short) order, pseudo );
	    for ( int i = 0; i < SS.length; i ++ ) {
	       SS[ i ] = new Sequence( bg.GenerateSequence( SS[ i ].GetLength() ), 
				       Sequence.CreateFakeGeneName() );
	    }
	 }
      }
      return SS;
   }

   public static Sequence[] RemoveDups( Sequence SS[], boolean verbose ) {
      int[] dups = IntUtils.New( SS.length ); IntUtils.Set( dups, -1 );
      int ndups = 0;
      for ( int i = 0; i < SS.length; i ++ ) {
	 if ( SS[ i ] == null ) continue;
	 if ( dups[ i ] >= 0 ) continue;
	 for ( int j = i + 1; j < SS.length; j ++ ) {
	    if ( SS[ j ] == null ) continue;
	    if ( dups[ j ] >= 0 ) continue;
	    if ( SS[ i ].GetSequence().equals( SS[ j ].GetSequence() ) ) {
	       dups[ j ] = i;
	       ndups ++;
	       if ( verbose ) System.err.println( SS[ i ].GetHeader() + " and " +
						  SS[ j ].GetHeader() + " are duplicates." );
	    } 
	 }
      }
      if ( ndups > 0 ) {
	 Sequence[] newS = new Sequence[ SS.length-ndups ];
	 int j = 0;
	 for ( int i = 0; i < SS.length; i ++ ) if ( dups[ i ] == -1 ) newS[ j ++ ] = SS[ i ];
	 return newS;
      }
      return SS;
   }

   public ObjVector InsertMotifs( String insert ) {
      String[] inserts = MyUtils.Tokenize( insert, "," );
      int ninsert = inserts.length;
      int nper = 1, ind = 0; double gapFreq = 0.0;
      ObjVector inserted = new ObjVector();
      try { 
	 nper = Integer.parseInt( inserts[ ind ] ); ind ++; ninsert --;
	 if ( nper <= 0 || nper > 10 ) nper = 1;
	 try { 
	    gapFreq = Double.parseDouble( inserts[ ind ] ) / 100.0; ind ++; ninsert --;
	    if ( gapFreq < 0.0 || gapFreq > 1.0 ) gapFreq = 0.0;
	 } catch( NumberFormatException e ) { gapFreq = 0.0; }
      } catch( NumberFormatException e ) { nper = 1; }

      int l = GetLength();
      char[] sseq = GetSequence().toCharArray();
      short[] flagged = null;
      if ( nper > 1 ) flagged = ShortUtils.New( l );
      int toInsert = ind, pos = 0, ii = 0;
      for ( int iii = 0; iii < nper; iii ++ ) {
	 if ( ninsert > 1 ) toInsert = IntUtils.RandChoose( ind, ninsert + ind - 1 );
	 int lllen = inserts[ toInsert ].length();

	 int sssum = 0;
	 if ( inserts[ toInsert ].indexOf( '[' ) >= 0 ) {
	    int xind = inserts[ toInsert ].indexOf( '[' ), 
	       xind2 = inserts[ toInsert ].indexOf( ']', xind );
	    sssum += ( xind2 - xind );
	    while ( inserts[ toInsert ].indexOf( '[', xind2 ) >= 0 ) {
	       xind = inserts[ toInsert ].indexOf( '[', xind2 );
	       xind2 = inserts[ toInsert ].indexOf( ']', xind );
	       sssum += ( xind2 - xind );
	    }
	 }
	 lllen -= sssum;
	 pos = IntUtils.RandChoose( 0, l - lllen );
	 if ( flagged != null )
	    while( flagged[ pos ] > 0 ) pos = IntUtils.RandChoose( 0, l - lllen );
	 lllen += sssum;

	 if ( pos < 0 ) pos = 0;
	 ii = pos;
	 
	 boolean letterOrHydrophobic = false;
	 char[] ttoInsert = inserts[ toInsert ].toCharArray();
	 for ( int j = 0; j < lllen; j ++ ) {
	    if ( ii >= l ) break;
	    if ( j >= lllen ) break;
	    if ( sseq[ ii ] == ' ' ) { ii ++; continue; }
	    if ( gapFreq > 0 && DoubleUtils.Random() < gapFreq ) {
	       if ( DoubleUtils.Random() < 0.5 ) { // Insert a random letter (a random insert)
		  if ( j > 0 && ii < l-1 ) {
		     sseq[ ii ++ ] = GetRandomResidue();
		     j --;
		     continue;
		  }
	       } else { // Skip this letter in the motif (a deletion)
		  j ++;
		  continue;
	       }
	    }
	    if ( letterOrHydrophobic ) {
	       if ( DoubleUtils.Random() < 0.5 ) sseq[ ii ++ ] = 
						    Character.toUpperCase( ttoInsert[ j ] );
	       else sseq[ ii ++ ] = CharUtils.RandChoose( PROT_HYDROPHOBIC );
	       letterOrHydrophobic = false;
	    } else if ( ttoInsert[ j ] == 'x' || ttoInsert[ j ] == 'X' ) {
	       sseq[ ii ++ ] = GetRandomResidue();
	    } else if ( ttoInsert[ j ] == '#' ) {
	       sseq[ ii ++ ] = CharUtils.RandChoose( PROT_HYDROPHOBIC );
	    } else if ( ttoInsert[ j ] == '@' ) {
	       sseq[ ii ++ ] = CharUtils.RandChoose( PROT_AROMATIC );
	    } else if ( ttoInsert[ j ] == '+' ) {
	       sseq[ ii ++ ] = CharUtils.RandChoose( PROT_CHARGED_POS );
	    } else if ( ttoInsert[ j ] == '-' ) {
	       sseq[ ii ++ ] = CharUtils.RandChoose( PROT_CHARGED_NEG );
	    } else if ( ttoInsert[ j ] == '[' ) {
	       j ++;
	       int jj = j; 
	       while( jj < inserts[ toInsert ].length() && ttoInsert[ jj ] != ']' ) jj ++;
	       int jlen = jj - j;
	       char letter = ttoInsert[ j + CharUtils.RandChoose( 0, jlen-1 ) ];
	       if ( letter >= 'a' && letter <= 'y' ) {
		  if ( DoubleUtils.Random() < 0.5 ) sseq[ ii ++ ] = Character.toUpperCase( letter );
		  else sseq[ ii ++ ] = GetRandomResidue(); 
	       } else { sseq[ ii ++ ] = letter; }
	       j = jj;
	    } else if ( ttoInsert[ j ] == 'Z' ) {
	       letterOrHydrophobic = true;
	    } else if ( ttoInsert[ j ] >= 'a' && ttoInsert[ j ] <= 'y' ) {
	       if ( DoubleUtils.Random() < 0.5 ) sseq[ ii ++ ] = 
						    Character.toUpperCase( ttoInsert[ j ] );
	       else sseq[ ii ++ ] = GetRandomResidue(); 
	    } else if ( ttoInsert[ j ] >= 'A' && ttoInsert[ j ] <= 'Y' ) {
	       if ( DoubleUtils.Random() < 0.9 ) sseq[ ii ++ ] = ttoInsert[ j ]; 
	       else sseq[ ii ++ ] = GetRandomResidue(); 
	    }
	 }
	 if ( flagged != null ) {
	    int llen = ii - pos;
	    for ( int iiii = pos - llen; iiii < ii; iiii ++ )
	       if ( iiii >= 0 && iiii < l ) flagged[ iiii ] = 1;
	 }
	 Initialize( new String( sseq ), GetHeader() );
	 try {
	    inserted.addElement( GetSequence().substring( pos, ii ) );
	    inserted.addElement( new Integer( pos + 1 ) );
	 } catch ( Exception e) { };
      }
      return inserted;
   }

   public void ForceProtein() {
      if ( type == PROTEIN ) return;
      type = PROTEIN;
      Initialize( GetSequence(), GetHeader() );
   }

   public Sequence Reverse() {
      Initialize( GetReverse(), new String( GetHeader() ) );
      return this;
   }

   public Sequence Complement() {
      if ( type == PROTEIN ) return this;
      Initialize( GetComplement(), new String( GetHeader() ) );
      return this;
   }

   public Sequence Append( Sequence s ) {
      Initialize( GetSequence() + s.GetSequence(), GetHeader() + " + " + s.GetHeader() );
      return this;
   }

   public Sequence Replace( Sequence s, int where ) {
      short inres[] = s.GetResidues(), outres[] = GetResidues();
      System.arraycopy( inres, 0, outres, where, s.GetLength() );
      String str = stringFromResidues( outres, GetType() );
      Initialize( str, GetHeader() );
      return this;
   }

   public String Subseq( int start, int end ) {
      if ( end != -999 && end < start ) { return ""; }
      if ( start > GetLength() || ( end != -999 && start > end ) ) return null;
      if ( start < 0 ) start = 0;
      if ( end < 0 ) end = GetLength() - 1;
      if ( end >= GetLength() ) end = GetLength() - 1;
      int newlen = end - start + 1;
      if ( newlen > GetLength() ) return null;
      return seq.substring( start, start+newlen );
   }

   public int[] GetCounts() {
      int[] out = IntUtils.New( GetAlphabetSize() );
      short[] res = GetResidues();
      for ( int i = 0; i < GetLength(); i ++ ) out[ res[ i ] ] ++;
      return out;
   }

   public static double[] GetHydrophobicities( short res[], String name, double into[] ) {
      if ( into == null ) into = DoubleUtils.New( res.length );
      double hydro[] = (double[]) HYDROPHOBICITIES.get( name );
      for ( int i = 0, s = res.length; i < s; i ++ ) into[ i ] = hydro[ res[ i ] ];
      return into;
   }

   public static double[] GetHydrophobicities( short res[], double into[] ) {
      return GetHydrophobicities( res, "KyteDoolittle", into );
   }

   public static double[] GetHydrophobicities( String name, short res[], int win, 
					       double into[] ) {
      double temp[] = DoubleUtils.New( res.length );
      double hydro[] = (double[]) HYDROPHOBICITIES.get( name );
      for ( int i = 0, s = res.length; i < s; i ++ ) temp[ i ] = hydro[ res[ i ] ];
      double kern[] = DoubleUtils.New( win * 4 + 1 );
      for ( int i = 0, s = kern.length; i < s; i ++ ) // Convolve with a 2*win-sigma gaussian 
	 kern[ i ] = cern.jet.stat.Probability.normal( (double) (win*2+1), (double) win, 
						       (double) i );
      if ( into == null ) into = DoubleUtils.New( res.length );
      DoubleUtils.Convolve( temp, kern, into );
      return into;
   }

   public static double[] GetHydrophobicities( short res[], int win, double into[] ) {
      return GetHydrophobicities( "KyteDoolittle", res, win, into );
   }

   public double[] GetHydrophobicities( double into[] ) {
      return GetHydrophobicities( "KyteDoolittle", into );
   }

   public double[] GetHydrophobicities( String name, double into[] ) {
      if ( type != PROTEIN ) return into;
      return GetHydrophobicities( GetResidues(), name, into );
   }

   public double[] GetHydrophobicities( int win, double into[] ) {
      return GetHydrophobicities( "KyteDoolittle", win, into );
   }

   public double[] GetHydrophobicities( String name, int win, double into[] ) {
      if ( type != PROTEIN ) return into;
      return GetHydrophobicities( name, GetResidues(), win, into );
   }

   public double[][] AddOrRemoveCountsAtSite( boolean add, int site, int width ) {
      double pssm[][] = DoubleUtils.New( width, GetAlphabetSize() );
      return AddOrRemoveCountsAtSite( add, site, pssm );
   }

   public double[][] AddOrRemoveCountsAtSite( boolean add, int site, double pssm[][] ) {
      int W = pssm.length, len = GetLength(), max = len - W + 1;
      double ad = add ? 1.0 : -1.0;
      short seq[] = GetResidues();
      for ( int i = 0; i < W; i ++ ) {
         if ( site + i >= len ) break;
         int j = (int) seq[ site + i ];
         double[] ci = pssm[ i ];
         ci[ j ] += ad;
         if ( ci[ j ] < 0 ) ci[ j ] = 0;
      }
      return pssm;
   }

   public static String dnaToProt( String dna ) {
      StringBuffer out = new StringBuffer(), in = new StringBuffer( dna );
      for ( int i = 0, len = in.length()-3; i <= len; i += 3 )
	 out.append( CODONS.get( in.substring( i, i + 3 ) ) );
      return out.toString();
   }

   public Sequence ToProtein() {
      if ( type == PROTEIN ) return this;
      Initialize( dnaToProt( GetSequence() ), GetHeader() + " (PROT-CONVERTED)" );
      return this;
   }
   
   public static Sequence[] RunTest( String argv[] ) {
      gnu.getopt.Getopt g = new gnu.getopt.Getopt( "SeqTest", argv, "f:o:i:" );
      String fname = "test.fst", out = null, insert = null;
      int c;
      while ( ( c = g.getopt() ) != -1 ) {
	 switch( c ) {
	 case 'f': fname = g.getOptarg(); break;
	 case 'o': out = g.getOptarg(); break;
	 case 'i': insert = g.getOptarg(); break;
	 case '?': break; // getopt() already printed an error
	 default: System.out.print( "getopt() returned " + c + "\n" );
	 }
      }
 
      System.out.println( "READING FILE " + fname );
      Sequence[] seqs = null;
      try {
	 seqs = ReadSequences( fname );
      } catch( Exception e ) {
	 System.err.println( "File not found...\n" + e );
	 MyUtils.Exit( 0 );
      }
      if ( seqs != null ) {
	 for ( int i = 0; i < seqs.length; i ++ ) {
	    System.out.println( (i+1) + " " + seqs[ i ] );
	    if ( insert != null ) seqs[ i ].InsertMotifs( insert );
	    System.out.println( (i+1) + " " + seqs[ i ] );
	 }
      }
      if ( out != null ) {
	 System.out.println( "WRITING TO FILE " + out );
	 WriteSequences( out, seqs );
      }

      // START GUI
      //double h1[] = seqs[ 0 ].GetHydrophobicities( null );
      //djr.util.gui.MyPlot.Plot( "Hydrophobicities", "Sequence Pos.", "Hydrophob.", h1 );
      //h1 = seqs[ 0 ].GetHydrophobicities( 6, h1 );
      //djr.util.gui.MyPlot.Plot( "Hydrophobicities", "Sequence Pos.", "Hydrophob.", h1 );
      double h1[] = seqs[ 0 ].GetHydrophobicities( 15, null );
      djr.util.gui.MyPlot.Plot( "Hydrophobicities", "Sequence Pos.", "Hydrophob.", h1 );
      h1 = seqs[ 0 ].GetHydrophobicities( 30, h1 );
      djr.util.gui.MyPlot.Plot( "Hydrophobicities", "Sequence Pos.", "Hydrophob.", h1 );
      // END GUI

      return seqs;
   }

   /*public static void main( String argv[] ) {
     RunTest( argv );
     }
   */
}
