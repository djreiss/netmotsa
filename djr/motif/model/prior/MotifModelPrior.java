package djr.motif.model.prior;
import djr.util.array.*;
import djr.motif.*;
import djr.motif.model.*;

/**
 * Class <code>MotifModelPrior</code>
 *
 * @author <a href="mailto:reiss@uw.edu">David Reiss</a>
 * @version 1.9978 (Fri Nov 07 05:56:26 PST 2003)
 */
public abstract class MotifModelPrior implements java.io.Serializable {
   public abstract void GetProbs( double[] cts, double[] probs );

   public abstract int GetAlphabetSize();

   public double[] GetProbs( double[] cts ) {
      double[] out = DoubleUtils.New( cts.length );
      GetProbs( cts, out );
      return out;
   }

   public static final MotifModelPrior GetPrior( short J, String fgType, String dirichName,
						 String matrixName, double pseudo, 
						 BackgroundModel bgmodel ) {
      MotifModelPrior fgmodel = null;
      if ( J == 20 && ( "dirichlet".equals( fgType ) || "dirich".equals( fgType ) ) ) 
	 fgmodel = new DirichletMixturePrior( dirichName, pseudo, bgmodel );
      else if ( "avgscore".equals( fgType ) || "matrix".equals( fgType ) ) 
	 fgmodel = new AverageScorePrior( new ScoringMatrix( matrixName, J ), pseudo, bgmodel );
      else {
	 if ( ! "pseudo".equals( fgType ) ) System.err.println( "WARNING: " + fgType + 
		     " foreground model requested but pseudocount model being used." );
	 fgmodel = new PseudocountPrior( pseudo, bgmodel );
      }
      return fgmodel;
   }
}
