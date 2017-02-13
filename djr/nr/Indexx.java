package djr.nr;

import djr.util.array.*;

/**
 * Class <code>Indexx</code>
 *
 * @author <a href="mailto:reiss@uw.edu">David Reiss</a>
 * @version 1.9978 (Fri Nov 07 05:56:26 PST 2003)
 */
public class Indexx {
   static final int M = 7;
   static final int NSTACK = 50;

   static final void SWAP( int a[], int ai, int b[], int bi ) {
      int itemp = a[ai]; a[ai] = b[bi]; b[bi] = itemp;
   }

   public static final int[] irank( double arr[], int rnk[] ) throws NRException {
      int indx[] = null; indx = indexx( arr, indx );
      if ( rnk == null || rnk.length < arr.length ) rnk = IntUtils.New( arr.length );
      int j, n = arr.length; for (j=0;j<n;j++) rnk[ indx[ j ] ] = j; 
      return rnk;
   }

   public static final int[] indexx( double arr[], int indx[] ) throws NRException {
      int n=arr.length,i,indxt,ir=n-1,itemp,j,k,l=0;
      int jstack=0,istack[];
      double a;

      if ( indx == null || indx.length < arr.length ) indx = IntUtils.New( arr.length );
      istack=IntUtils.New(NSTACK);
      for (j=0;j<n;j++) indx[j]=j;
      for (;;) {
	 if (ir-l < M) {
	    for (j=l+1;j<=ir;j++) {
	       indxt=indx[j];
	       a=arr[indxt];
	       for (i=j-1;i>=0;i--) {
		  if (arr[indx[i]] <= a) break;
		  indx[i+1]=indx[i];
	       }
	       indx[i+1]=indxt;
	    }
	    if (jstack == 0) break;
	    ir=istack[jstack--];
	    l=istack[jstack--];
	 } else {
	    k=(l+ir) >> 1;
	    SWAP(indx,k,indx,l+1);
	    if (arr[indx[l+1]] > arr[indx[ir]]) SWAP(indx,l+1,indx,ir);
	    if (arr[indx[l]] > arr[indx[ir]]) SWAP(indx,l,indx,ir);
	    if (arr[indx[l+1]] > arr[indx[l]]) SWAP(indx,l+1,indx,l);
	    i=l+1;
	    j=ir;
	    indxt=indx[l];
	    a=arr[indxt];
	    for (;;) {
	       do i++; while (arr[indx[i]] < a);
	       do j--; while (arr[indx[j]] > a);
	       if (j < i) break;
	       SWAP(indx,i,indx,j);
	    }
	    indx[l]=indx[j];
	    indx[j]=indxt;
	    jstack += 2;
	    if (jstack > NSTACK) throw new NRException("NSTACK too small in indexx.");
	    if (ir-i+1 >= j-l) {
	       istack[jstack]=ir;
	       istack[jstack-1]=i;
	       ir=j-1;
	    } else {
	       istack[jstack]=j-1;
	       istack[jstack-1]=l;
	       l=i;
	    }
	 }
      }
      return indx;
   }

   /*public static void main( String args[] ) {
      double arr[] = DoubleUtils.New( 10 );
      DoubleUtils.Random( arr );
      DoubleUtils.Mult( arr, 10.0 );
      int ind[] = DoubleUtils.Indexx( arr, null );
      DoubleUtils.Printf( "%.3f ", arr );
      IntUtils.Print( ind );
      DoubleUtils.Printf( "%.3f ", DoubleUtils.Index( arr, ind ) );
      }*/
}
