package djr.motif.model.prior;
import djr.util.array.*;
import djr.motif.*;
import djr.motif.model.*;

/**
 * Class <code>AverageScorePrior</code>
 *
 * @author <a href="mailto:reiss@uw.edu">David Reiss</a>
 * @version 1.9978 (Fri Nov 07 05:56:26 PST 2003)
 */
public class AverageScorePrior extends PseudocountPrior {
   protected ScoringMatrix sm;
   protected boolean isDiag;

   public AverageScorePrior( ScoringMatrix mat, double pseudo, int alphSize ) {
      super( pseudo, alphSize );
      this.sm = mat;
      this.isDiag = sm.IsDiagonal();
   }

   public AverageScorePrior( ScoringMatrix mat, double pseudo, BackgroundModel bg ) {
      super( pseudo, bg );
      this.sm = mat;
      this.isDiag = sm.IsDiagonal();
   }

   public void GetProbs( double[] cts, double[] probs ) {
      if ( isDiag ) {
	 super.GetProbs( cts, probs );
	 return;
      }

      double sum = 0.0, sumcounts = DoubleUtils.Sum( cts );
      for ( short i = 0; i < alphSize; i ++ ) {
	 probs[ i ] = 0;
	 for ( short j = 0; j < alphSize; j ++ ) probs[ i ] += sm.GetRawScore( i, j ) * cts[ j ]; 
	 probs[ i ] = counts[ i ] * Math.pow( 2.0, probs[ i ] / sumcounts ) + 
	    pseudocounts[ i ] + sumcounts * cts[ i ];
	 sum += probs[ i ];
      }
      DoubleUtils.Divide( probs, sum );
   }   
}
