package djr.motif.sampler.prior;
import djr.motif.sampler.sites.Sites;
import djr.motif.sampler.Sampler;

/**
 * Interface <code>SitesPrior</code>
 *
 * @author <a href="mailto:dreiss@systemsbiology.org">David Reiss</a>
 * @version 1.9978 (Fri Nov 07 05:56:26 PST 2003)
 */
public abstract class SitesPrior implements java.io.Serializable {
   public SitesPrior() { };
   public SitesPrior( Sampler samp ) { };
   public abstract double GetPriorValue( Sites sites, short seq[], int site, int mot );
}
