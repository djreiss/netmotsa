package djr.motif.model;
import djr.util.bio.Sequence;
import djr.util.array.*;
import djr.util.MyUtils;

/**
 * <code>BackgroundModel</code> class.
 *
 * @author <a href="mailto:dreiss@systemsbiology.org">David Reiss</a>
 * @version 1.9978 (Fri Nov 07 05:56:26 PST 2003)
 */
public class BackgroundModel implements java.io.Serializable {
   double freq4[][][][], freq3[][][], freq2[][], freq1[], pseudo;
   int counts4[][][][], counts3[][][], counts2[][], counts1[];
   short J, order;
   long sumCounts1, sumCounts2[], sumCounts3[][], sumCounts4[][][];

   public BackgroundModel( String seqFile, short ord, double pseudo ) {
      Sequence bgS[] = null;
      if ( ! MyUtils.IsNullString( seqFile ) ) {
	 try {
	    bgS = Sequence.ReadSequences( seqFile );
	 } catch( Exception e ) {
	    bgS = null;
	    System.err.println( "File '" + seqFile + 
				"' not found. Using input sequences to generate bg model.\n" );
	 }
      }
      if ( bgS != null ) {
	 J = bgS[ 0 ].GetAlphabetSize();
	 int ssumLen = 0;
	 for ( int ii = 0; ii < bgS.length; ii ++ ) 
	    if ( bgS[ ii ] != null ) ssumLen += bgS[ ii ].GetLength();
	 if ( Math.pow( J, ord ) > ssumLen * 10 )
	    System.err.println( "WARNING: Background model order is probably too high for " +
				"amount of sequence data provided. Continuing anyway." );
      }
      Initialize( bgS, ord, pseudo );
   }

   public BackgroundModel( Sequence[] seqs, short ord ) {
      this( seqs, ord, 0.0 );
   }

   public BackgroundModel( Sequence[] seqs, short ord, double pseud ) {
      Initialize( seqs, ord, pseud );
   }

   public int GetOrder() {
      return order;
   }

   protected void Initialize( Sequence[] seqs, short ord, double pseud ) {
      pseudo = pseud;
      Sequence S[] = seqs;
      J = S[ 0 ].GetAlphabetSize();
      order = ord;
      order ++;

      counts1 = IntUtils.New( J );
      freq1 = DoubleUtils.New( J );
      if ( order > 1 ) {
	 counts2 = IntUtils.New( J, J );
	 freq2 = DoubleUtils.New( J, J );
	 sumCounts2 = LongUtils.New( J );
      }
      if ( order > 2 ) {
	 counts3 = IntUtils.New( J, J, J );
	 freq3 = DoubleUtils.New( J, J, J );
	 sumCounts3 = LongUtils.New( J, J );
      }
      if ( order > 3 ) {
	 counts4 = IntUtils.New( J, J, J, J );
	 freq4 = DoubleUtils.New( J, J, J, J );
	 sumCounts4 = LongUtils.New( J, J, J );
      }
      FillTheModel( S );
   }

   public short GetAlphabetSize() { return J; }

   protected void FillTheModel( Sequence S[] ) {
      for ( int i = 0; i < S.length; i ++ ) {
	 short[] res = S[ i ].GetResidues();
	 for ( int j = 0; j < S[ i ].GetLength(); j ++ ) {
	    if ( ! S[ i ].IsMasked( j ) ) {
	       counts1[ res[ j ] ] ++;
	       sumCounts1 ++;
	    }
	 }
      }
      for ( int i = 0; i < J; i ++ ) {
	 freq1[ i ] = ( (double) counts1[ i ] + pseudo * sumCounts1 ) / 
	    ( ( 1.0 + pseudo ) * sumCounts1 );
      }
   
      if ( order > 1 ) {
	 for ( int i = 0; i < S.length; i ++ ) {
	    short[] res = S[ i ].GetResidues(); int len = S[ i ].GetLength();
	    for ( int j = 1; j < len; j ++ ) {
	       if ( ! S[ i ].IsMaskedInRange( j - 1, j ) ) {
		  counts2[ res[ j-1 ] ][ res[ j ] ] ++;
		  sumCounts2[ res[ j-1 ] ] ++;
	       }
	    }
	 }
	 for ( int i = 0; i < J; i ++ ) {
	    double sum = 0.0;
	    for ( int j = 0; j < J; j ++ ) {
	       freq2[ i ][ j ] = ( (double) counts2[ i ][ j ] + pseudo * sumCounts2[ i ] ) / 
		  ( ( 1.0 + pseudo ) * sumCounts2[ i ] ); // * freq1[ i ];
	       sum += freq2[ i ][ j ];
	    }
	    DoubleUtils.Divide( freq2[ i ], sum );
	 }
      }

      if ( order > 2 ) {
	 for ( int i = 0; i < S.length; i ++ ) {
	    short[] res = S[ i ].GetResidues(); int len = S[ i ].GetLength();
	    for ( int j = 2; j < len; j ++ ) {
	       if ( ! S[ i ].IsMaskedInRange( j - 2, j ) ) {
		  counts3[ res[ j-2 ] ][ res[ j-1 ] ][ res[ j ] ] ++;
		  sumCounts3[ res[ j-2 ] ][ res[ j-1 ] ] ++;
	       }
	    }
	 }
	 for ( int i = 0; i < J; i ++ ) {
	    for ( int j = 0; j < J; j ++ ) {
	       double sum = 0.0;
	       for ( int k = 0; k < J; k ++ ) {
		  freq3[ i ][ j ][ k ] = ( (double) counts3[ i ][ j ][ k ] + 
					   pseudo * sumCounts3[ i ][ j ] ) / 
		     ( ( 1.0 + pseudo ) * sumCounts3[ i ][ j ] ); // * freq2[ i ][ j ];
		  sum += freq3[ i ][ j ][ k ];
	       }
	       DoubleUtils.Divide( freq3[ i ][ j ], sum );
	    }
	 }
      }

      if ( order > 3 ) {
	 for ( int i = 0; i < S.length; i ++ ) {
	    short[] res = S[ i ].GetResidues(); int len = S[ i ].GetLength();
	    for ( int j = 3; j < len; j ++ ) {
	       if ( ! S[ i ].IsMaskedInRange( j - 3, j ) ) {
		  counts4[ res[ j-3 ] ][ res[ j-2 ] ][ res[ j-1 ] ][ res[ j ] ] ++;
		  sumCounts4[ res[ j-3 ] ][ res[ j-2 ] ][ res[ j-1 ] ] ++;
	       }
	    }
	 }
	 for ( int i = 0; i < J; i ++ ) {
	    for ( int j = 0; j < J; j ++ ) {
	       for ( int k = 0; k < J; k ++ ) {
		  double sum = 0.0;
		  for ( int l = 0; l < J; l ++ ) {
		     freq4[ i ][ j ][ k ][ l ] = ( (double) counts4[ i ][ j ][ k ][ l ] + 
						   pseudo * sumCounts4[ i ][ j ][ k ] ) / 
			( ( 1.0 + pseudo ) * sumCounts4[ i ][ j ][ k ] ); //*freq3[ i ][ j ][ k ];
		     sum += freq4[ i ][ j ][ k ][ l ];
		  }
		  DoubleUtils.Divide( freq4[ i ][ j ][ k ], sum );
	       }
	    }
	 }
      }
   }

   public double GetProb( int res ) {
      return freq1[ res ];
   }

   public double GetProb( short res[], int ind ) {
      if ( res[ ind ] < 0 ) return 0.0;
      if ( order <= 1 || ind == 0 ) {
	 return freq1[ res[ ind ] ];
      } else if ( order == 2 || ind == 1 ) {
	 if ( res[ ind - 1 ] < 0 ) return 0.0;
	 return freq2[ res[ ind - 1 ] ][ res[ ind ] ];
      } else if ( order == 3 || ind == 2 ) {
	 if ( res[ ind - 2 ] < 0 || res[ ind - 1 ] < 0 ) return 0.0;
	 return freq3[ res[ ind - 2 ] ][ res[ ind - 1 ] ][ res[ ind ] ];
      } else if ( ind >= 3 ) {
	 if ( res[ ind - 3 ] < 0 || res[ ind - 2 ] < 0 || res[ ind - 1 ] < 0 ) return 0.0;
	 return freq4[ res[ ind - 3 ] ][ res[ ind - 2 ] ][ res[ ind - 1 ] ][ res[ ind ] ];
      }
      return 0.0;
   }

   public long GetTotalCounts() { 
      return sumCounts1; 
   }

   public int GetCounts( int res ) {
      return counts1[ res ];
   }

   public int GetCounts( short res[], int ind ) {
      if ( res[ ind ] < 0 ) return 0;
      if ( order <= 1 || ind == 0 ) {
	 return counts1[ res[ ind ] ];
      } else if ( order == 2 || ind == 1 ) {
	 if ( res[ ind - 1 ] < 0 ) return 0;
	 return counts2[ res[ ind - 1 ] ][ res[ ind ] ];
      } else if ( order == 3 || ind == 2 ) {
	 if ( res[ ind - 2 ] < 0 || res[ ind - 1 ] < 0 ) return 0;
	 return counts3[ res[ ind - 2 ] ][ res[ ind - 1 ] ][ res[ ind ] ];
      } else if ( ind >= 3 ) {
	 if ( res[ ind - 3 ] < 0 || res[ ind - 2 ] < 0 || res[ ind - 1 ] < 0 ) return 0;
	 return counts4[ res[ ind - 3 ] ][ res[ ind - 2 ] ][ res[ ind - 1 ] ][ res[ ind ] ];
      }
      return 0;
   }

   public short[] GenerateSequence( int length ) {
      short out[] = ShortUtils.New( length );
      double probs[] = DoubleUtils.New( J );
      for ( int i = 0; i < length; i ++ ) {
	 for ( int j = 0; j < J; j ++ ) {
	    out[ i ] = (short) j;
	    probs[ j ] = GetProb( out, i );
	 }
	 short res = (short) DoubleUtils.Sample( probs, J );
	 if ( res >= 0 ) out[ i ] = res;
	 else out[ i ] = ShortUtils.RandChoose( 0, J-1 );
      }
      return out;
   }
}
