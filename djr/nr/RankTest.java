package djr.nr;

import djr.util.array.*;

public class RankTest {
   public double mean1, mean2, alldata[], t = 0, prob = 0;
   public int ranks[];

   public RankTest() { };

   public RankTest( double data1[], double data2[] ) throws NRException {
      recompute( data1, data2 );
   }

   public void recompute( double data1[], double data2[] ) throws NRException {
      int n1 = data1.length, n2 = data2.length, ntot = n1 + n2;
      if ( n1 == 0 || n2 == 0 ) { t = -Double.MAX_VALUE; prob = 0.0; return; }
      if ( alldata == null || alldata.length < ntot ) {
	 alldata = DoubleUtils.New( ntot );
	 ranks = IntUtils.New( ntot );
      }
      System.arraycopy( data1, 0, alldata, 0, n1 );
      System.arraycopy( data2, 0, alldata, n1, n2 );
      Indexx.irank( alldata, ranks );
      mean1 = IntUtils.Mean( ranks, 0, n1 - 1 );
      mean2 = IntUtils.Mean( ranks, n1, ntot - 1 );
      t = ( mean2 - mean1 ) / ntot / Math.sqrt( ( (double) ntot + 1 ) / ( ( (double) 12 * n1 * n2 ) ) );
      double df = ntot - 2; // Degrees of freedom. 
      prob = TTest.betai( 0.5 * df, 0.5, df / ( df + t * t ) ); // See equation (6.4.9) in NR
   }

   public double getTValue() { return t; }
   public double getProb() { return prob; }   
}
