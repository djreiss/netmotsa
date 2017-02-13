package djr.nr;

import djr.util.array.*;

/**
 * Class <code>Chisq2</code>
 * Use the chstwo method from Numerical Recipes
 *
 * @author <a href="mailto:reiss@uw.edu">David Reiss</a>
 * @version 1.9978 (Fri Nov 07 05:56:26 PST 2003)
 */
public class Chisq2 {
   double df, chsq, prob;

   public Chisq2( double bins1[], double bins2[], int knstrn ) throws NRException {
      int j, nbins = bins1.length;
      double temp;
      
      df=nbins-knstrn;
      chsq=0.0;
      for (j=0;j<nbins;j++) {
	 if (bins1[j] == 0.0 && bins2[j] == 0.0)
	    --df;
	 else {
	    temp=bins1[j]-bins2[j];
	    chsq += temp*temp/(bins1[j]+bins2[j]);
	 }
      }
      prob=gammq(0.5*(df),0.5*(chsq));
   }

   public double getProb() { return prob; }

   double gammq(double a, double x) throws NRException {
      double gamser = 0,gammcf = 0,gln = 0, args[] = new double[ 2 ];
      if (x < 0.0 || a <= 0.0) throw new NRException("Invalid arguments in routine gammq");
      if (x < (a+1.0)) {
	 args[0] = gamser; args[1] = gln;
	 gser(args,a,x);
	 return 1.0-args[0];
      } else {
	 args[0] = gammcf; args[1] = gln;
	 gcf(args,a,x);
	 return args[0];
      }
   }

   final int ITMAX = 100;
   final double EPS = 3.0e-7;
   final double FPMIN = 1.0e-30;

   void gcf(double args[],/*double *gammcf,*/ double a, double x/*, double *gln*/) throws NRException {
	int i;
	double an,b,c,d,del,h;

	args[1]=gammln(a);
	b=x+1.0-a;
	c=1.0/FPMIN;
	d=1.0/b;
	h=d;
	for (i=1;i<=ITMAX;i++) {
	   an = -i*(i-a);
	   b += 2.0;
	   d=an*d+b;
	   if (Math.abs(d) < FPMIN) d=FPMIN;
	   c=b+an/c;
	   if (Math.abs(c) < FPMIN) c=FPMIN;
	   d=1.0/d;
	   del=d*c;
	   h *= del;
	   if (Math.abs(del-1.0) < EPS) break;
	}
	if (i > ITMAX) throw new NRException("a too large, ITMAX too small in gcf");
	args[0]=Math.exp(-x+a*Math.log(x)-(args[1]))*h;
   }

   void gser(double args[],/*double *gamser,*/ double a, double x/*, double *gln*/) throws NRException {
      int n;
      double sum,del,ap;

      args[1]=gammln(a);
      if (x <= 0.0) {
	 if (x < 0.0) throw new NRException("x less than 0 in routine gser");
	 args[0]=0.0;
	 return;
      } else {
	 ap=a;
	 del=sum=1.0/a;
	 for (n=1;n<=ITMAX;n++) {
	    ++ap;
	    del *= x/ap;
	    sum += del;
	    if (Math.abs(del) < Math.abs(sum)*EPS) {
	       args[0]=sum*Math.exp(-x+a*Math.log(x)-(args[1]));
	       return;
	    }
	 }
	 throw new NRException("a too large, ITMAX too small in routine gser");
      } 
   }

   static double cof[]={76.18009172947146,-86.50532032941677,
			24.01409824083091,-1.231739572450155,
			0.1208650973866179e-2,-0.5395239384953e-5};

   double gammln(double xx) {
      double x,y,tmp,ser;
      int j;

      y=x=xx;
      tmp=x+5.5;
      tmp -= (x+0.5)*Math.log(tmp);
      ser=1.000000000190015;
      for (j=0;j<=5;j++) ser += cof[j]/++y;
      return -tmp+Math.log(2.5066282746310005*ser/x);
   }
}
