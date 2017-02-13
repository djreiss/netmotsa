package djr.motif.sampler;

import djr.util.array.*;
import java.util.Vector;

/**
 * Abstract class <code>PostProcessor</code>
 *
 * @author <a href="mailto:reiss@uw.edu">David Reiss</a>
 * @version 1.9978 (Fri Nov 07 05:56:26 PST 2003)
 */
public abstract class PostProcessor {
   protected Sampler sampler;

   public PostProcessor( Sampler samp ) {
      this.sampler = samp;
   }

   public abstract void run();
}
