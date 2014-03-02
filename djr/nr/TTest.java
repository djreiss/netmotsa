package djr.nr;

import djr.util.array.*;

/**
 * Class <code>TTest</code>
 * Use the ttest method from Numerical Recipes
 *
 * @author <a href="mailto:dreiss@systemsbiology.org">David Reiss</a>
 * @version 1.9978 (Fri Nov 07 05:56:26 PST 2003)
 */
public class TTest {
   protected double t = 0, prob = 0, temp[] = DoubleUtils.New( 2 );

   /* Given the arrays data1[1..n1] and data2[1..n2], this routine returns 
      Student's t as t, and its signi cance as prob, small values of prob 
      indicating that the arrays have significantly different means. 
      The data arrays are assumed to be drawn from populations with the 
      same true variance. */

   public TTest() { };

   public TTest( double data1[], double data2[] ) throws NRException {
      recompute( data1, data2 );
   }

   public void recompute( double data1[], double data2[] ) throws NRException {
      int n1 = data1.length, n2 = data2.length;
      if ( n1 <= 0 || n2 <= 0 || n1 + n2 <= 2 ) { 
	 t = -Double.MAX_VALUE; prob = 0.0; return; }
      double ave1 = 0, var1 = 0, ave2 = 0, var2 = 0;
      if ( n1 > 1 ) {
	 avevar( data1, temp );
	 ave1 = temp[ 0 ]; var1 = temp[ 1 ];
	 if ( n2 == 1 ) {
	    t = ( data2[ 0 ] - ave1 ) / var1 / Math.sqrt( 2 ); // use t-statistic for 1 data pt
	    double df = n1 - 1;
	    prob=betai(0.5*df,0.5,df/(df+t*t));
	    return;
	 }
      }
      if ( n2 > 1 ) {
	 avevar( data2, temp );
	 ave2 = temp[ 0 ]; var2 = temp[ 1 ];
	 if ( n1 == 1 ) {
	    t = ( data1[ 0 ] - ave2 ) / var2 / Math.sqrt( 2 ); // use t-statistic for 1 data pt
	    double df = n2 - 1;
	    prob=betai(0.5*df,0.5,df/(df+t*t));
	    return;
	 }
      }
      //System.err.println("MEAN1="+ave1+"  MEAN2="+ave2+"  VAR1="+var1+"  VAR2="+var2);
      double df=n1+n2-2; // Degrees of freedom. 
      double svar=((n1-1)*var1+(n2-1)*var2)/df; // Pooled variance. 
      t=(ave2-ave1)/Math.sqrt(svar*(1.0/n1+1.0/n2)); 
      prob=betai(0.5*df,0.5,df/(df+t*t)); // See equation (6.4.9). 
   }

   public double getTValue() { return t; }
   public double getProb() { return prob; }

   /* Given array data[1..n], returns its mean as ave and its variance as var. */
   static void avevar(double data[], double outavevar[] ) {
      int j, n = data.length; double ave, var, ep; 
      for (ave=0.0,j=0;j<n;j++) ave += data[j]; 
      ave /= n; 
      var=ep=0.0; 
      for (j=0;j<n;j++) { 
	 double s=data[j]-ave; 
	 ep += s; 
	 var += s*s; 
      } 
      var=(var-ep*ep/n)/(n-1); // Corrected two-pass formula (14.1.8). 
      outavevar[ 0 ] = ave; outavevar[ 1 ] = var;
   }

   /* Returns the incomplete beta function I x (a; b). */
   public static double betai(double a, double b, double x) throws NRException { 
      double bt = 0.0; 
      if (x < 0.0 || x > 1.0) throw new NRException("Bad x in routine betai"); 
      if (x == 0.0 || x == 1.0) bt=0.0; 
      else // Factors in front of the continued fraction. 
	 bt=Math.exp(DoubleUtils.logGamma(a+b)-DoubleUtils.logGamma(a)-
		     DoubleUtils.logGamma(b)+a*Math.log(x)+b*Math.log(1.0-x)); 
      if (x < (a+1.0)/(a+b+2.0)) // Use continued fraction directly. 
	 return bt*betacf(a,b,x)/a; 
      else // Use continued fraction after making the symmetry transformation. 
	 return 1.0-bt*betacf(b,a,1.0-x)/b; 
   } 

   static final int MAXIT = 100;
   static final double EPS = 3.0e-7;
   static final double FPMIN = 1.0e-30;

   /* Used by betai: Evaluates continued fraction for incomplete beta function by 
      modified Lentz's method ( x 5.2). */
   static double betacf(double a, double b, double x) throws NRException { 
      double qab=a+b; // These q's will be used in factors that occur in the coefficients (6.4.6).
      double qap=a+1.0; 
      double qam=a-1.0; 
      double c=1.0; // First step of Lentz's method. 
      double d=1.0-qab*x/qap; 
      if (Math.abs(d) < FPMIN) d=FPMIN; 
      d=1.0/d; 
      double h=d; 
      int m;
      for (m=1;m<=MAXIT;m++) { 
	 int m2=2*m; 
	 double aa=m*(b-m)*x/((qam+m2)*(a+m2));
	 d=1.0+aa*d; // One step (the even one) of the recurrence. 
	 if (Math.abs(d) < FPMIN) d=FPMIN; 
	 c=1.0+aa/c; 
	 if (Math.abs(c) < FPMIN) c=FPMIN; 
	 d=1.0/d; 
	 h *= d*c; 
	 aa = -(a+m)*(qab+m)*x/((a+m2)*(qap+m2)); 
	 d=1.0+aa*d; // Next step of the recurrence (the odd one). 
	 if (Math.abs(d) < FPMIN) d=FPMIN; 
	 c=1.0+aa/c; 
	 if (Math.abs(c) < FPMIN) c=FPMIN; 
	 d=1.0/d; 
	 double del=d*c; 
	 h *= del; 
	 if (Math.abs(del-1.0) < EPS) break; // Are we done? 
      } 
      if (m > MAXIT) throw new NRException("a or b too big, or MAXIT too small in betacf"); 
      return h;
   }

   /*public static void main( String args[] ) {
      double d1[] = DoubleUtils.New( 100 );
      DoubleUtils.Normal( d1, 0.0, 3.0 );
      DoubleUtils.PrintStats( d1 );
      System.out.println();
      djr.util.gui.MyPlot.PlotHistogram( "d1", "x", "count", d1 );

      double d2[] = DoubleUtils.New( 100 );
      DoubleUtils.Normal( d2, 10.0, 30.0 );
      DoubleUtils.PrintStats( d2 );
      System.out.println();
      djr.util.gui.MyPlot.PlotHistogram( "d2", "x", "count", d2 );

      try { 
	 TTest ttest = new TTest( d1, d2 );
	 System.err.println( ttest.getTValue()+" "+ttest.getProb());
      } catch( Exception e ) { e.printStackTrace(); }
   }
   */   
}
