package djr.motif.model;
import java.util.Enumeration;
import djr.motif.sampler.*;
import djr.util.array.*;
import djr.util.MyUtils;
import djr.util.bio.Sequence;

/**
 * <code>ScoringMatrix</code> class.
 *
 * @author <a href="mailto:dreiss@systemsbiology.org">David Reiss</a>
 * @version 1.9978 (Fri Nov 07 05:56:26 PST 2003)
 */
public class ScoringMatrix implements java.io.Serializable {
   // BG freqs (percent) used for PAM matrices
   final static double pamFreqs[] = { 8.69e-2, 3.30e-2, 4.70e-2, 5.00e-2, 4.00e-2, 
			       8.89e-2, 3.40e-2, 3.70e-2, 8.09e-2, 8.49e-2, 
			       1.50e-2, 4.00e-2, 5.09e-2, 3.80e-2, 4.10e-2, 
			       6.99e-2, 5.79e-2, 6.49e-2, 1.00e-2, 3.00e-2 };
   // BG freqs (percent) from SwissProt database used for BLOSUM matrices
   final static double spFreqs[] =  { 7.55e-2, 1.70e-2, 5.30e-2, 6.32e-2, 4.07e-2, 
			       6.86e-2, 2.23e-2, 5.73e-2, 5.95e-2, 9.32e-2, 
			       2.36e-2, 4.55e-2, 4.92e-2, 4.03e-2, 5.16e-2, 
			       7.19e-2, 5.77e-2, 6.52e-2, 1.26e-2, 3.21e-2 };

   boolean isDiagonal;
   short J;
   double[][] S; // The scoring matrix (raw log-odds score)
   double[][] Spow; // 2^(S/2) (S "un-logged")
   double[][] freq; // The scoring matrix converted to a mutation frequency
   double[][] logFreq; // The log2 of the frequency matrix
   double[][] normFreq; // The frequency table, normalized so each column is 1
   double[][] logNormFreq; // The frequency table, normalized so each column is 1, in log2
   Sequence residues;
   String matrix;

   public double[][] GetScoreMatrix() { return S; }
   public double[][] GetPowScoreMatrix() { return Spow; }
   public double[][] GetFrequencyMatrix() { return freq; }
   public double[][] GetLogFrequencyMatrix() { return logFreq; }
   public double[][] GetNormFrequencyMatrix() { return normFreq; }   
   public double[][] GetLogNormFrequencyMatrix() { return logNormFreq; }   
   public boolean IsDiagonal() { return isDiagonal; }
   public Sequence GetResidues() { return residues; }
   public void Print() { DoubleUtils.Printf( "%.3f ", normFreq ); }
   public static final double SwissProtFreq( int j ) { return spFreqs[ j ]; }
   public static final double PAMFreq( int j ) { return pamFreqs[ j ]; }
   public double GetScore( short i, short j ) { return normFreq[ i ][ j ]; }
   public double GetRawScore( short i, short j ) { return S[ i ][ j ]; }

   public ScoringMatrix( String fname, short j ) {
      J = j;
      if ( fname != null ) matrix = new String( fname );
      else matrix = "";
      if ( matrix.equals( "" ) ) {
	 if ( J == 4 ) matrix = "matrix/DNA.txt";
	 else if ( J == 20 ) matrix = "matrix/PROT.txt";
      }
      S = Spow = freq = logFreq = normFreq = logNormFreq = null;
      residues = null;
      readLUT();
      double sum = 0;
      for ( int i = 0; i < J; i ++ ) sum += normFreq[ i ][ i ];
      isDiagonal = ( sum == J );
   }

   void readLUT() {
      S = DoubleUtils.New( J, J );
      Spow = DoubleUtils.New( J, J );
      freq = DoubleUtils.New( J, J );
      logFreq = DoubleUtils.New( J, J );
      normFreq = DoubleUtils.New( J, J );
      logNormFreq = DoubleUtils.New( J, J );
      for ( int i = 0; i < J; i ++ ) freq[ i ][ i ] = normFreq[ i ][ i ] = 1.0;

      try {
	 int line = 0;
	 Enumeration lines = MyUtils.ReadFileLines( matrix ).elements();
	 while( lines.hasMoreElements() && line < J ) {
	    String str = (String) lines.nextElement();
	    if ( str == null || str.startsWith( "#" ) || str.startsWith( "//" ) || 
		   str.equals( "" ) ) continue;
	    if ( str.charAt( 0 ) == ' ' ) {
	       int i = 0, ind = 0;
	       char[] fptr = str.toCharArray();
	       char[] temp = new char[ J ];
	       while( i < J ) {
		  ind ++;
		  if ( fptr[ ind ] == ' ' ) continue;
		  else if ( J == 20 && ( fptr[ ind ] == 'B' || fptr[ ind ] == 'J' || 
					 fptr[ ind ] == 'O' || fptr[ ind ] == 'U' ||
					 fptr[ ind ] == 'X' || fptr[ ind ] >= 'Z' || 
					 fptr[ ind ] < 'A' ) ) continue;
		  temp[ i ++ ] = fptr[ ind ];
	       }
	       residues = new Sequence( new String( temp ) );
	    }
	    if ( residues.GetAlphabetSize() != J ) {
	       System.err.println( "ERROR: maxtrix size " + residues.GetAlphabetSize() +
				   " is incompatible with sequence alphabet size " + J + 
				   ". Quitting ungracefully." );
	       MyUtils.Exit( 1 );
	    }
	    char[] fptr = str.toCharArray();
	    if ( fptr[0] < 'A' || fptr[0] > 'Z' ) continue;
	    int ind = 0;
	    while( ( fptr[ ind ] < '0' || fptr[ ind ]  > '9' ) && fptr[ ind ]  != '+' && 
		   fptr[ ind ]  != '-' && fptr[ ind ]  != '.' ) ind ++;
	    double[] f = DoubleUtils.FTokenize( str.substring( ind ), " " );
	    short[] res = residues.GetResidues();
	    int ind2 = res[ line ];
	    for ( int i = 0; i < J; i ++ ) {
	       int ind1 = res[ i ];
	       S[ ind2 ][ ind1 ] = f[ i ];
	    }
	    line ++;
	 }
      } catch( Exception e ) {
         System.err.println( "Could not load matrix file " + matrix );
         e.printStackTrace();
      }

      // Matrices are defined as int(2*log2(P(i,j)/e(i,j))) where e(i,j)=P(i)*P(j).
      // Need to convert to acid frequencies (if protein). Use the SwissProt bg freqs
      // (assumes we're using BLOSUM matrices).
      double[] sums = DoubleUtils.New( J );
      for ( int i = 0; i < J; i ++ ) {
	 for ( int j = 0; j < J; j ++ ) {
	    Spow[ i ][ j ] = Math.pow( 2.0, S[ i ][ j ] / 2.0 );
	    if ( J == 20 ) {
	       if ( i == j ) freq[ i ][ j ] = Spow[ i ][ j ] * spFreqs[i] * spFreqs[j];
	       else freq[ i ][ j ] = Spow[ i ][ j ] * 2.0 * spFreqs[i] * spFreqs[j];
	    }
	    logFreq[ i ][ j ] = DoubleUtils.Log2( freq[ i ][ j ] );
	    sums[ i ] += freq[ i ][ j ];
	 }
      }
      for ( int i = 0; i < J; i ++ ) {
	 for ( int j = 0; j < J; j ++ ) {
	    normFreq[ i ][ j ] = freq[ i ][ j ] / sums[ i ];
	    logNormFreq[ i ][ j ] = DoubleUtils.Log2( normFreq[ i ][ j ] );
	 }
      }
   }
}
