package djr.motif.model.prior;

import djr.util.array.*;
import djr.motif.*;
import djr.motif.model.*;

/**
 * Class <code>PseudocountPrior</code>
 *
 * @author <a href="mailto:dreiss@systemsbiology.org">David Reiss</a>
 * @version 1.9978 (Fri Nov 07 05:56:26 PST 2003)
 */
public class PseudocountPrior extends MotifModelPrior {
   protected double pseudo, pseudocounts[];
   protected BackgroundModel bg;
   protected int counts[], alphSize;

   public PseudocountPrior( double pseudo, int alphSize ) {
      this.pseudo = pseudo;
      this.bg = bg;
      this.alphSize = alphSize;
      counts = IntUtils.New( alphSize );
      pseudocounts = DoubleUtils.New( alphSize );
      for ( int i = 0; i < alphSize; i ++ ) {
	 counts[ i ] = 1;
	 pseudocounts[ i ] = counts[ i ] * pseudo;
      }
   }

   public PseudocountPrior( double pseudo, BackgroundModel bg ) {
      this.pseudo = pseudo;
      this.bg = bg;
      this.alphSize = bg.GetAlphabetSize();
      counts = IntUtils.New( alphSize );
      pseudocounts = DoubleUtils.New( alphSize );
      for ( int i = 0; i < alphSize; i ++ ) {
	 counts[ i ] = bg.GetCounts( i );
	 pseudocounts[ i ] = counts[ i ] * pseudo;
      }
   }

   public int GetAlphabetSize() { return alphSize; }

   public void GetProbs( double[] cts, double[] probs ) {
      double sum = 0.0;
      for ( int i = 0; i < alphSize; i ++ ) {
	 probs[ i ] = cts[ i ] + pseudocounts[ i ];
	 sum += probs[ i ];
      }
      DoubleUtils.Divide( probs, sum );
   }
}
