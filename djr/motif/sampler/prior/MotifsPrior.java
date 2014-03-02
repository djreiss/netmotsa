package djr.motif.sampler.prior;
import djr.motif.sampler.sites.*;
import djr.motif.sampler.Sampler;
import djr.util.array.*;

/**
 * Class <code>MotifsPrior</code>
 *
 * @author <a href="mailto:dreiss@systemsbiology.org">David Reiss</a>
 * @version 1.9978 (Fri Nov 07 05:56:26 PST 2003)
 */
public class MotifsPrior extends SitesPrior {
   double Ni, ai, Ai, WW;

   // Needs to be implemented correctly
   public MotifsPrior( Sampler samp ) {
      super( samp );
   }
   
   public MotifsPrior( int NS, int Ni, double WW ) {
      this.WW = WW;
      this.Ni = (double) Ni;
      this.ai = NS * WW; // / ( 1 - WW ); 
      this.Ai = Ni * WW; // / ( 1 - WW ); 
   }

   public double GetPriorValue( Sites sites, short seq[]/*ignored*/, 
				int site/*ignored*/, int mot ) {
      return GetPriorValue( sites, mot );
   }

   public double GetPriorValue( Sites sites, int mot ) {
      double pp = ( sites.GetMotifCounts( mot ) + ai ) / ( Ni + Ai );
      return pp; // / ( 1.0 - pp );
   }
}
