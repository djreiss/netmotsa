package djr.util.array;

import java.util.*;
import corejava.*;
import cern.jet.stat.*;
import cern.jet.random.Normal;
//import cern.jet.random.engine.MersenneTwister; //DRand;
import edu.cornell.lassp.houle.RngPack.Ranecu;
import djr.util.MyUtils;

/**
 * Template class <code>BoolUtils</code>.
 *
 * This is a "template" class and should not be directly compiled.
 * Commands to create real file before compiling (e.g. for ints) are:
 * cp Utils.java IntUtils.java
 * rpl Boolean Integer IntUtils.java ; rpl Bool Int IntUtils.java ; rpl boolean int IntUtils.java
 * javac IntUtils.java
 *
 * @author <a href="mailto:reiss@uw.edu">David Reiss</a>
 * @version 1.9978 (Fri Nov 07 05:56:26 PST 2003)
 */
public class BoolUtils {
   public static double version = 1.9978;
   //public static MersenneTwister rand = null;
   public static Ranecu rand = null;
   public static Normal nrand = null;
   public static final double LOG2 = Math.log( 2.0 );
   public static final double LOG10 = Math.log( 10.0 );

   static { SetSeed(); }

   // All the typeUtils classes share the rand in DoubleUtils.
   public static final void SetSeed() { 
      if ( DoubleUtils.rand == null ) SetSeed( System.currentTimeMillis() );
      else rand = DoubleUtils.rand; }
   public static final void SetSeed( long seed ) { 
      DoubleUtils.rand = rand = new Ranecu( seed );
      Random(); } // The first number is always 0.9997something, for some reason
   /* NO BOOL
   public static final boolean Random() {
      return (boolean) rand.raw(); }
   public static final long GetSeed() {
      return DoubleUtils.rand.getSeed(); }
   public static final void Random( double arr[] ) { 
      DoubleUtils.rand.raw( arr, arr.length ); }
   public static final boolean Random( int max ) { 
      return (boolean) ( DoubleUtils.Random() * max ); }
   public static final boolean Random( double min, double max ) { 
      return (boolean) ( min + ( max - min ) * DoubleUtils.Random() ); }
   public static final boolean Random( int min, int max ) { 
      return Random( (double) min, (double) max ); }
   public static final boolean RandChoose( int lo, int hi ) {
      return (boolean) ( (long) lo + (long) ( ( 1L + (long) hi - (long) lo ) * 
					   DoubleUtils.Random() ) ); }
   public static final boolean RandChoose( int hi ) {
      return RandChoose( 0, hi ); }
   END NO BOOL */
   public static final boolean RandChoose( boolean[] in ) {
      return in[ IntUtils.RandChoose( 0, in.length - 1 ) ]; }

   // BOOL ONLY
   public static final boolean Random() {
      return rand.raw() > 0.5 ? false : true; }
   // END BOOL ONLY */

   /* NO BOOL
   public static final boolean Normal( double mean, double stddev ) { 
      return (boolean) Normal.staticNextDouble( mean, stddev ); }
   public static final boolean Normal() { 
      return Normal( 0.0, 1.0 ); }

   public static final double Log( boolean in ) { return in <= 0 ? -999 : Math.log( in ); }
   public static final double Log2( boolean in ) { return DoubleUtils.Log( in ) / LOG2; }
   public static final double Log10( boolean in ) { return DoubleUtils.Log( in ) / LOG10; }
   public static final double Sign( boolean in ) { return in < 0 ? -1.0 : 1.0; }
   END NO BOOL */

   public static final boolean[] New( int nx ) { 
      return new boolean[ nx ]; }
   public static final boolean[][] New( int nx, int ny ) { 
      return new boolean[ nx ][ ny ]; }
   public static final boolean[][][] New( int nx, int ny, int nz ) { 
      return new boolean[ nx ][ ny ][ nz ]; }
   public static final boolean[][][][] New( int nx, int ny, int nz, int na ) { 
      return new boolean[ nx ][ ny ][ nz ][ na ]; }

   public static final boolean[] Resize( boolean[] in, int newx ) {
      boolean[] out = New( newx ); Copy( out, in ); return out; }
   public static final boolean[][] Resize( boolean[][] in, int newx, int newy ) {
      boolean[][] out = New( newx, newy ); Copy( out, in ); return out; }
   public static final boolean[][][] Resize( boolean[][][] in, int newx, int newy, int newz ) {
      boolean[][][] out = New( newx, newy, newz ); Copy( out, in ); return out; }
   public static final boolean[][][][] Resize( boolean[][][][] in, int newx, int newy, int newz, int newa ) {
      boolean[][][][] out = New( newx, newy, newz, newa ); Copy( out, in ); return out; }

   public static final boolean[] New( boolean[] in ) {
      if ( in == null ) return null;
      boolean[] out = New( in.length ); Copy( out, in ); return out; }
   public static final boolean[][] New( boolean[][] in ) {
      boolean[][] out = new boolean[ in.length ][];
      for ( int i = 0, s = in.length; i < s; i ++ ) out[ i ] = New( in[ i ] ); return out; }
   public static final boolean[][][] New( boolean[][][] in ) {
      boolean[][][] out = new boolean[ in.length ][][];
      for ( int i = 0, s = in.length; i < s; i ++ ) out[ i ] = New( in[ i ] ); return out; }

   public static final Object NewObj( Object obj ) {
      if ( obj instanceof boolean[] ) return New( (boolean[]) obj );
      else if ( obj instanceof boolean[][] ) return New( (boolean[][]) obj );
      else if ( obj instanceof boolean[][][] ) return New( (boolean[][][]) obj );
      return null; }

   public static final boolean[] SubArr( boolean[] in, int start, int end ) {
      boolean[] out = New( end - start + 1 );
      System.arraycopy( in, start, out, 0, end - start + 1 );
      return out; }

   /* NO BOOL
   public static final boolean[] Sequence( int length ) {
      boolean[] out = New( length ); for ( int i = 0; i < length; i ++ ) out[ i ] = true;
      return out; }
   public static final double[] Convert( Object in ) {
      if ( in instanceof double[] ) return (double[]) in;
      boolean iin[] = (boolean[]) in;
      double out[] = new double[ iin.length ];
      for ( int i = 0, size = iin.length; i < size; i ++ ) out[ i ] = (double) iin[ i ];
      return out; }
   END NO BOOL */

   public static final void Print( Object in ) {
      if ( in instanceof boolean[] ) Print( (boolean[]) in, 0, ( (boolean[]) in ).length ); 
      else for ( int i = 0, size = ( (Object[]) in ).length; i < size; i ++ ) 
	 Print( ( (Object[]) in )[ i ] ); }
   public static final void Print( String start, Object in ) {
      System.out.print( start ); Print( in ); }

   public static final void Print( boolean[] in, int start, int end ) {
      if ( in == null ) System.out.println( "NULL" );
      for ( int i = start; i < end; i ++ ) System.out.print( in[ i ] + " " );
      System.out.println(); }

   /* NO BOOL
   public static final void Printf( String fmt, boolean in ) {
      if ( Double.isInfinite( (double) in ) ) in = - Boolean.MAX_VALUE;
      Format.print( System.out, fmt, in ); };
   public static final String SPrintf( String fmt, boolean in ) {
      if ( Double.isInfinite( (double) in ) ) in = - Boolean.MAX_VALUE;
      return ( new Format( fmt ) ).form( in ); }

   public static final void PrintCols( String fmt, boolean in1[], boolean in2[] ) {
      PrintCols( fmt, new boolean[][]{ in1, in2 } ); }
   public static final void PrintCols( String fmt, boolean in[][], boolean lines ) {
      for ( int i = 0, size = in[ 0 ].length; i < size; i ++ ) {
	 if ( lines ) System.out.print( i + " " );
	 for ( int j = 0, ss = in.length; j < ss; j ++ ) 
	    System.out.print( SPrintf( fmt, in[ j ][ i ] ) );
	 System.out.println(); } }
   public static final void PrintCols( String fmt, boolean in[][] ) {
      PrintCols( fmt, in, false ); }

   public static final void Printf( String fmt, Object in ) {
      if ( in == null ) System.out.println( "NULL" );
      if ( in instanceof boolean[] ) { for ( int i = 0, size = ( (boolean[]) in ).length; i < size; i ++ )
	 Printf( fmt, ( (boolean[]) in )[ i ] ); }
      else { for ( int i = 0, size = ( (Object[]) in ).length; i < size; i ++ )
	 Printf( fmt, ( (Object[]) in )[ i ] ); }
      System.out.println(); }
   public static final void Printf( String start, String fmt, Object in ) {
      System.out.print( start ); Printf( fmt, in ); }

   public static final String SPrintf( String fmt, Object in ) {
      if ( in == null ) return "NULL";
      String out = "";
      if ( in instanceof boolean[] ) { for ( int i = 0, size = ( (boolean[]) in ).length; i < size; i ++ )
	 out += SPrintf( fmt, ( (boolean[]) in )[ i ] ); }
      else { for ( int i = 0, size = ( (Object[]) in ).length; i < size; i ++ )
	 out += SPrintf( fmt, ( (Object[]) in )[ i ] ); }
      return out; }

   public static final void PrintT( String fmt, boolean[][] in ) {
      for ( int y = 0, size = in[0].length; y < size; y ++ ) {
	 for ( int x = 0, s = in.length; x < s; x ++ ) Format.print( System.out, fmt, in[ x ][ y ] );
	 System.out.println(); } }
   public static final void PrintT( boolean[][] in ) {
      for ( int y = 0, size = in[0].length; y < size; y ++ ) {
	 for ( int x = 0, s = in.length; x < s; x ++ ) System.out.print( in[ x ][ y ] + " " );
	 System.out.println(); } }
   END NO BOOL */

   public static final boolean[][] Transpose( boolean[][] in ) {
      boolean out[][] = New( in[ 0 ].length, in.length );
      for ( int i = 0; i < in[ 0 ].length; i ++ ) for ( int j = 0; j < in.length; j ++ ) 
	 out[ i ][ j ] = in[ j ][ i ]; return out; }

   public static final void Copy( Object dest, Object src ) {
      if ( src == null || dest == null ) return;
      if ( src instanceof boolean[] ) System.arraycopy( (boolean[]) src, 0, (boolean[]) dest, 0, 
			IntUtils.Min( ( (boolean[]) src ).length, ( (boolean[]) dest ).length ) );
      else for ( int i = 0, s = IntUtils.Min( ( (Object[]) src ).length, 
	     ( (Object[]) dest ).length ); i < s; i ++ ) 
	 Copy( ( (Object[]) dest )[ i ], ( (Object[]) src )[ i ] );
   }

   public static final boolean[] Flatten( boolean[][] in, boolean[] out ) {
      if ( out == null || out.length < NElem( in ) ) out = new boolean[ NElem( in ) ]; 
      int ind = 0; for ( int i = 0, s = in.length; i < s; i ++ ) {
	 System.arraycopy( in[ i ], 0, out, ind, in[ i ].length ); ind += in[ i ].length; }
      return out; }
   public static final boolean[] Flatten( boolean[][] in ) {
      return Flatten( in, null ); }

   public static final boolean[] Flatten( boolean[][][] in, boolean[] out ) {
      if ( out == null || out.length < NElem( in ) ) out = new boolean[ NElem( in ) ]; 
      int ind = 0; for ( int i = 0, s = in.length; i < s; i ++ ) {
	 boolean temp[] = Flatten( in[ i ] );
	 System.arraycopy( temp, 0, out, ind, temp.length ); ind += temp.length; }
      return out; }
   public static final boolean[] Flatten( boolean[][][] in ) {
      return Flatten( in, null ); }

   private static boolean[][] UnflattenPriv( boolean in[], int start, int nx, int ny, boolean out[][] ) {
      if ( out == null || out.length < nx ) out = new boolean[ nx ][]; int ind = start;
      for ( int i = 0; i < nx; i ++ ) {
	 if ( out[ i ] == null || out[ i ].length < ny ) out[ i ] = new boolean[ ny ];
	 System.arraycopy( in, ind, out[ i ], 0, ny ); ind += ny;
      } return out; }
   public static boolean[][] Unflatten( boolean in[], int nx, int ny, boolean out[][] ) {
      return UnflattenPriv( in, 0, nx, ny, out ); }
   public static boolean[][] Unflatten( boolean in[], int nx, int ny ) {
      return UnflattenPriv( in, 0, nx, ny, null ); }

   public static boolean[][][] Unflatten( boolean in[], int nx, int ny, int nz, boolean out[][][] ) {
      if ( out == null || out.length < nx ) out = new boolean[ nx ][][]; int ind = 0, step = ny *nz;
      for ( int i = 0; i < nx; i ++ ) { out[ i ] = UnflattenPriv( in, ind, ny, nz, out[ i ] );
      ind += step; } return out; }
   public static boolean[][][] Unflatten( boolean in[], int nx, int ny, int nz ) {
      return Unflatten( in, nx, ny, nz, null ); }

   public static final void Zero( Object in, int start, int end ) {
      if ( in == null ) return;
      if ( in instanceof boolean[] ) Arrays.fill( (boolean[]) in, start, end, (boolean) false );
      else for ( int i = 0, s = ( (Object[]) in ).length; i < s; i ++ ) 
	 Zero( ( (Object[]) in )[ i ], start, end ); }
   public static final void Zero( Object in ) {
      if ( in == null ) return;
      if ( in instanceof boolean[] ) Arrays.fill( (boolean[]) in, (boolean) false );
      else for ( int i = 0, s = ( (Object[]) in ).length; i < s; i ++ ) 
	 Zero( ( (Object[]) in )[ i ] ); }
   public static final void Zero( boolean in[], int inds[] ) {
      if ( in == null ) return;
      else for ( int i = 0, s = inds.length; i < s; i ++ ) in[ inds[ i ] ] = false;
   }

   public static final void Set( Object in, int from, int to, boolean val ) {
      if ( in instanceof boolean[] ) Arrays.fill( (boolean[]) in, from, to, val );
      else for ( int i = 0, s = ( (Object[]) in ).length; i < s; i ++ ) 
	 Set( ( (Object[]) in )[ i ], from, to, val ); }
   public static final void Set( Object in, boolean val ) {
      if ( in == null ) return;
      if ( in instanceof boolean[] ) Arrays.fill( (boolean[]) in, val );
      else for ( int i = 0, s = ( (Object[]) in ).length; i < s; i ++ ) 
	 Set( ( (Object[]) in )[ i ], val ); }
   public static final void Set( boolean in[], int inds[], boolean val ) {
      if ( in == null ) return;
      else for ( int i = 0, s = inds.length; i < s; i ++ ) in[ inds[ i ] ] = val;
   }

   /* NO BOOL
   public static final boolean Sum( Object in, int min, int max ) {
      if ( in == null ) return false;
   END NO BOOL */
   // BOOL ONLY
   public static final int Sum( Object in, int min, int max ) {
   // END BOOL ONLY */
      if ( in instanceof boolean[] ) { boolean[] inn = (boolean[]) in; 
	 if ( max >= inn.length ) max = inn.length - 1;
   /* NO BOOL
	 boolean out = 0; for ( int i = min; i <= max; i ++ ) out += inn[ i ]; return out; 
   END NO BOOL */
   // BOOL ONLY
         int out = 0; for ( int i = min; i <= max; i ++ ) out += inn[ i ] ? 1 : 0; 
	 return out; 
   // END BOOL ONLY */
      } else {
   /* NO BOOL
	 boolean out = 0; 
   END NO BOOL */
   // BOOL ONLY
         int out = 0; 
   // END BOOL ONLY */
	 for ( int i = 0, s = ( (Object[]) in ).length; i < s; i ++ ) 
	    out += Sum( ( (Object[]) in )[ i ], min, max ); return out;
      } }

   /* NO BOOL
   public static final boolean Sum( Object in ) {
   END NO BOOL */
   // BOOL ONLY
   public static final int Sum( Object in ) {
   // END BOOL ONLY */
      if ( in instanceof boolean[] ) return Sum( (boolean[]) in, 0, ( (boolean[]) in ).length ); 
      else {
   /* NO BOOL
	 boolean out = 0; 
   END NO BOOL */
   // BOOL ONLY
         int out = 0; 
   // END BOOL ONLY */
	 for ( int i = 0, s = ( (Object[]) in).length; i < s; i ++ ) 
	    out += Sum( ( (Object[]) in )[ i ] ); return out; }
   }

   /* NO BOOL
   public static final boolean SumSquared( Object in ) {
      if ( in == null ) return false;
      if ( in instanceof boolean[] ) { boolean[] inn = (boolean[]) in; boolean out = 0; 
      for ( int i = 0, s = inn.length; i < s; i ++ ) out += inn[ i ] * inn[ i ]; return out;
      } else {
	 boolean out = 0; for ( int i = 0, s = ( (Object[]) in).length; i < s; i ++ ) 
	    out += SumSquared( ( (Object[]) in )[ i ] ); return out; }
   }

   public static final boolean SumSquared2( Object in ) {
      if ( in == null ) return false;
      if ( in instanceof boolean[] ) { boolean[] inn = (boolean[]) in; boolean out = 0; 
      for ( int i = 0, s = inn.length; i < s; i ++ ) out += Sign( inn[ i ] ) * inn[ i ] * inn[ i ]; return out;
      } else {
	 boolean out = 0; for ( int i = 0, s = ( (Object[]) in).length; i < s; i ++ ) 
	    out += SumSquared2( ( (Object[]) in )[ i ] ); return out; }
   }

   public static final boolean Prod( Object in ) {
      if ( in == null ) return false;
      if ( in instanceof boolean[] ) { boolean[] inn = (boolean[]) in; boolean out = inn[ 0 ]; for ( int i = 1, s = inn.length; i < s; i ++ ) out *= inn[ i ]; return out;
      } else {
	 boolean out = Prod( ( (Object[]) in )[ 0 ] ); for ( int i = 1, s = ( (Object[]) in).length; i < s; i ++ ) 
	    out *= Prod( ( (Object[]) in )[ i ] ); return out; }
   }

   public static final boolean Max( Object in ) { boolean out = - Boolean.MAX_VALUE;
      if ( in == null ) return false;
      if ( in instanceof boolean[] ) { boolean[] inn = (boolean[]) in; 
      for ( int i = 0, s = inn.length; i < s; i ++ ) out = Max( out, inn[ i ] ); return out; 
      } else { for ( int i = 0, s = ( (Object[]) in ).length; i < s; i ++ ) 
	 out = Max( out, Max( ( (Object[]) in )[ i ] ) ); return out; } }
   
   public static final boolean Min( Object in ) { boolean out = Boolean.MAX_VALUE;
      if ( in == null ) return false;
      if ( in instanceof boolean[] ) { boolean[] inn = (boolean[]) in; 
      for ( int i = 0, s = inn.length; i < s; i ++ ) out = Min( out, inn[ i ] ); return out; 
      } else { for ( int i = 0, s = ( (Object[]) in ).length; i < s; i ++ ) 
	 out = Min( out, Min( ( (Object[]) in )[ i ] ) ); return out; } }

   public static final boolean[] MinMax( Object in, boolean inout[] ) { 
      if ( in == null ) return null;
      if ( inout == null ) { inout = New( 2 ); 
      inout[ 0 ] = Boolean.MAX_VALUE; inout[ 1 ] = - Boolean.MAX_VALUE; }
      if ( in instanceof boolean[] ) { boolean[] inn = (boolean[]) in; 
      for ( int i = 0, s = inn.length; i < s; i ++ ) {
	 inout[ 0 ] = Min( inout[ 0 ], inn[ i ] ); 
	 inout[ 1 ] = Max( inout[ 1 ], inn[ i ] ); } return inout;
      } else { for ( int i = 0, s = ( (Object[]) in ).length; i < s; i ++ ) 
	 inout = MinMax( ( (Object[]) in )[ i ], inout ); return inout; } }
   public static final boolean[] MinMax( Object in ) { 
      return MinMax( in, null ); }

   public static final void Round( Object in ) { 
      if ( in == null ) return;
      if ( in instanceof boolean[] ) { boolean[] inn = (boolean[]) in; 
      for ( int i = 0, s = inn.length; i < s; i ++ ) inn[ i ] = (boolean) (int) inn[ i ]; }
      else { for ( int i = 0, s = ( (Object[]) in ).length; i < s; i ++ ) 
	 Round( ( (Object[]) in )[ i ] ); }
   }

   public static final void Floor( Object in, boolean val ) { 
      if ( in == null ) return;
      if ( in instanceof boolean[] ) { boolean[] inn = (boolean[]) in; 
      for ( int i = 0, s = inn.length; i < s; i ++ ) if ( inn[ i ] < val ) inn[ i ] = val; }
      else { for ( int i = 0, s = ( (Object[]) in ).length; i < s; i ++ ) 
	 Floor( ( (Object[]) in )[ i ], val ); }
   }

   public static final void Ceil( Object in, boolean val ) { 
      if ( in == null ) return;
      if ( in instanceof boolean[] ) { boolean[] inn = (boolean[]) in; 
      for ( int i = 0, s = inn.length; i < s; i ++ ) if ( inn[ i ] > val ) inn[ i ] = val; }
      else { for ( int i = 0, s = ( (Object[]) in ).length; i < s; i ++ ) 
	 Ceil( ( (Object[]) in )[ i ], val ); }
   }
   
   public static final void Mult( Object in, boolean val ) { 
      if ( in == null ) return;
      if ( in instanceof boolean[] ) { boolean[] inn = (boolean[]) in; 
      for ( int i = 0, s = inn.length; i < s; i ++ ) inn[ i ] *= val; }
      else { for ( int i = 0, s = ( (Object[]) in ).length; i < s; i ++ ) 
	 Mult( ( (Object[]) in )[ i ], val ); }
   }

   public static final void Mult( Object in, Object in2 ) { 
      if ( in == null || in2 == null ) return;
      if ( in instanceof boolean[] ) { boolean[] inn = (boolean[]) in; 
      for ( int i = 0, s = inn.length; i < s; i ++ ) inn[ i ] *= ( (boolean[]) in2 )[ i ]; }
      else { for ( int i = 0, s = ( (Object[]) in ).length; i < s; i ++ ) 
	 Mult( ( (Object[]) in )[ i ], ( (Object[]) in2 )[ i ] ); }
   }

   public static final void Divide( Object in, boolean val ) { 
      if ( in == null ) return;
      if ( val == 0 ) return;
      if ( in instanceof boolean[] ) { boolean[] inn = (boolean[]) in; 
      for ( int i = 0, s = inn.length; i < s; i ++ ) inn[ i ] /= val; }
      else { for ( int i = 0, s = ( (Object[]) in ).length; i < s; i ++ ) 
	 Divide( ( (Object[]) in )[ i ], val ); }
   }

   public static final void Divide( boolean val, Object in ) { 
      if ( in == null ) return;
      if ( val == 0 ) return;
      if ( in instanceof boolean[] ) { boolean[] inn = (boolean[]) in; 
      for ( int i = 0, s = inn.length; i < s; i ++ ) inn[ i ] = (boolean) ( val / inn[ i ] ); }
      else { for ( int i = 0, s = ( (Object[]) in ).length; i < s; i ++ ) 
	 Divide( ( (Object[]) in )[ i ], val ); }
   }

   public static final void Divide( Object in, Object in2 ) { 
      if ( in == null || in2 == null ) return;
      if ( in instanceof boolean[] ) { boolean[] inn = (boolean[]) in; 
      for ( int i = 0, s = inn.length; i < s; i ++ ) inn[ i ] /= ( (boolean[]) in2 )[ i ]; }
      else { for ( int i = 0, s = ( (Object[]) in ).length; i < s; i ++ ) 
	 Divide( ( (Object[]) in )[ i ], ( (Object[]) in2 )[ i ] ); }
   }

   public static final void Add( Object in, boolean val ) { 
      if ( in == null ) return;
      if ( in instanceof boolean[] ) { boolean[] inn = (boolean[]) in; 
      for ( int i = 0, s = inn.length; i < s; i ++ ) inn[ i ] += val; }
      else { for ( int i = 0, s = ( (Object[]) in ).length; i < s; i ++ ) 
	 Add( ( (Object[]) in )[ i ], val ); }
   }

   public static final void Add( Object in, Object in2 ) { 
      if ( in == null || in2 == null ) return;
      if ( in instanceof boolean[] ) { boolean[] inn = (boolean[]) in; 
      for ( int i = 0, s = inn.length; i < s; i ++ ) inn[ i ] += ( (boolean[]) in2 )[ i ]; }
      else { for ( int i = 0, s = IntUtils.Min( ( (Object[]) in ).length, 
						( (Object[]) in2 ).length ); i < s; i ++ ) 
	 Add( ( (Object[]) in )[ i ], ( (Object[]) in2 )[ i ] ); }
   }

   public static final void Sub( Object in, boolean val ) { 
      Add( in, (boolean) -val ); }

   public static final void Sub( Object in, Object in2 ) { 
      if ( in == null || in2 == null ) return;
      if ( in instanceof boolean[] ) { boolean[] inn = (boolean[]) in; 
      for ( int i = 0, s = inn.length; i < s; i ++ ) inn[ i ] -= ( (boolean[]) in2 )[ i ]; }
      else { for ( int i = 0, s = ( (Object[]) in ).length; i < s; i ++ ) 
	 Sub( ( (Object[]) in )[ i ], ( (Object[]) in2 )[ i ] ); }
   }

   public static final void Invert( Object in ) {
      if ( in == null ) return;
      if ( in instanceof boolean[] ) { boolean[] inn = (boolean[]) in; 
      for ( int i = 0, s = inn.length; i < s; i ++ ) 
	 inn[ i ] = (boolean) ( 1.0 / (double) inn[ i ] ); }
      else { for ( int i = 0, s = ( (Object[]) in ).length; i < s; i ++ ) 
	 Invert( ( (Object[]) in )[ i ] ); }
   }

   END NO BOOL */

   public static final void Random( Object in ) {
      if ( in == null ) return;
      if ( in instanceof boolean[] ) { boolean[] inn = (boolean[]) in; 
      for ( int i = 0, s = inn.length; i < s; i ++ ) inn[ i ] = Random(); }
      else { for ( int i = 0, s = ( (Object[]) in ).length; i < s; i ++ ) 
	 Random( ( (Object[]) in )[ i ] ); }
   }

   /* NO BOOL
   public static final void Normal( Object in, double mean, double stddev ) {
      if ( in == null ) return;
      if ( nrand == null ) nrand = new Normal( mean, stddev, DoubleUtils.rand );
      else nrand.setState( mean, stddev );
      if ( in instanceof boolean[] ) { boolean[] inn = (boolean[]) in; 
      for ( int i = 0, s = inn.length; i < s; i ++ ) inn[ i ] = (boolean) nrand.nextDouble(); }
      else { for ( int i = 0, s = ( (Object[]) in ).length; i < s; i ++ ) 
	 Normal( ( (Object[]) in )[ i ], mean, stddev ); }
   }

   public static final void NormalAdd( Object in, double mean, double stddev ) {
      if ( in == null ) return;
      if ( nrand == null ) nrand = new Normal( mean, stddev, DoubleUtils.rand );
      else nrand.setState( mean, stddev );
      if ( in instanceof boolean[] ) { boolean[] inn = (boolean[]) in; 
      for ( int i = 0, s = inn.length; i < s; i ++ ) inn[ i ] += (boolean) nrand.nextDouble(); }
      else { for ( int i = 0, s = ( (Object[]) in ).length; i < s; i ++ ) 
	 NormalAdd( ( (Object[]) in )[ i ], mean, stddev ); }
   }

   public static final void Log10( Object in ) { 
      if ( in == null ) return;
      if ( in instanceof boolean[] ) { boolean[] inn = (boolean[]) in; 
      for ( int i = 0, s = inn.length; i < s; i ++ ) inn[ i ] = (boolean) Log10( inn[ i ] ); }
      else { for ( int i = 0, s = ( (Object[]) in ).length; i < s; i ++ ) 
	 Log10( ( (Object[]) in )[ i ] ); }
   }

   public static final void Log2( Object in ) { 
      if ( in == null ) return;
      if ( in instanceof boolean[] ) { boolean[] inn = (boolean[]) in; 
      for ( int i = 0, s = inn.length; i < s; i ++ ) inn[ i ] = (boolean) Log2( inn[ i ] ); }
      else { for ( int i = 0, s = ( (Object[]) in ).length; i < s; i ++ ) 
	 Log2( ( (Object[]) in )[ i ] ); }
   }

   public static final void Log( Object in ) { 
      if ( in == null ) return;
      if ( in instanceof boolean[] ) { boolean[] inn = (boolean[]) in; 
      for ( int i = 0, s = inn.length; i < s; i ++ ) inn[ i ] = (boolean) Log( inn[ i ] ); }
      else { for ( int i = 0, s = ( (Object[]) in ).length; i < s; i ++ ) 
	 Log( ( (Object[]) in )[ i ] ); }
   }

   public static final void Pow( boolean val, Object in ) { 
      if ( in == null ) return;
      if ( in instanceof boolean[] ) { boolean[] inn = (boolean[]) in; 
      for ( int i = 0, s = inn.length; i < s; i ++ ) 
	 inn[ i ] = (boolean) Math.pow( val, (double) inn[ i ] ); }
      else { for ( int i = 0, s = ( (Object[]) in ).length; i < s; i ++ ) 
	 Pow( ( (Object[]) in )[ i ], val ); }
   }

   public static final void Pow( Object in, boolean val ) { 
      if ( in == null ) return;
      if ( in instanceof boolean[] ) { boolean[] inn = (boolean[]) in; 
      for ( int i = 0, s = inn.length; i < s; i ++ ) 
	 inn[ i ] = (boolean) Math.pow( (double) inn[ i ], val ); }
      else { for ( int i = 0, s = ( (Object[]) in ).length; i < s; i ++ ) 
	 Pow( ( (Object[]) in )[ i ], val ); }
   }

   public static final void Pow2( Object in, boolean val ) { 
      if ( in == null ) return;
      if ( in instanceof boolean[] ) { boolean[] inn = (boolean[]) in; 
      for ( int i = 0, s = inn.length; i < s; i ++ ) 
	 inn[ i ] = (boolean) ( Sign( inn[ i ] ) * Math.pow( (double) inn[ i ], val ) ); }
      else { for ( int i = 0, s = ( (Object[]) in ).length; i < s; i ++ ) 
	 Pow2( ( (Object[]) in )[ i ], val ); }
   }

   public static final void Abs( Object in ) { 
      if ( in == null ) return;
      if ( in instanceof boolean[] ) { boolean[] inn = (boolean[]) in; 
      for ( int i = 0, s = inn.length; i < s; i ++ ) {
	 boolean x = inn[ i ]; inn[ i ] = (boolean) ( x > 0 ? x : -x ); } }
      else { for ( int i = 0, s = ( (Object[]) in ).length; i < s; i ++ ) 
	 Abs( ( (Object[]) in )[ i ] ); }
   }

   public static final void Replace( Object in, boolean from, boolean to ) {
      if ( in == null ) return;
      if ( in instanceof boolean[] ) { boolean[] inn = (boolean[]) in; 
      for ( int i = 0, s = inn.length; i < s; i ++ ) 
	 if ( inn[ i ] == from ) inn[ i ] = to; }
      else { for ( int i = 0, s = ( (Object[]) in ).length; i < s; i ++ ) 
	 Replace( ( (Object[]) in )[ i ], from, to ); }      
   }

   public static final int WhereMax( boolean[] in, int maxind ) { 
      if ( in == null ) return -1;
      boolean out = - Boolean.MAX_VALUE ; int ind = -1;
      maxind = Math.min( in.length, maxind+1 );
      for ( int i = 0; i < maxind; i ++ ) 
	 if ( out < in[ i ] ) { out = in[ i ]; ind = i; } return ind; }
   public static final int WhereMax( boolean[] in ) {
      return WhereMax( in, in.length ); }

   public static final int WhereMin( boolean[] in, int maxind ) { 
      boolean out = Boolean.MAX_VALUE ; int ind = -1;
      maxind = Math.min( in.length, maxind+1 );
      for ( int i = 0; i < maxind; i ++ ) if ( out > in[ i ] ) { out = in[ i ]; ind = i; } 
      return ind; }
   public static final int WhereMin( boolean[] in ) { 
      return WhereMin( in, in.length ); }

   public static final void MaxNorm( Object in ) { 
      boolean max = Max( in ); Divide( in, max ); }
   public static final void UnitNorm( Object in ) {
      boolean sum = Sum( in ); Divide( in, sum ); }

   public static final double Integral( boolean[] inx, boolean[] iny ) { // trapezoid rule
      double sum = 0.0;
      for ( int i = 1, s = inx.length; i < s; i ++ )
	 sum += ( inx[ i ] + inx[ i - 1 ] ) / 2.0 * ( iny[ i ] - iny[ i - 1 ] ) ;
      return sum;
   }

   public static final boolean[] Integrate( boolean[] in, boolean[] out ) {
      if ( out == null || out.length < in.length ) out = New( in.length );
      boolean sum = 0;
      for ( int i = 0, s = in.length; i < s; i ++ ) { sum += in[ i ]; out[ i ] = sum; }
      return out;
   }
   END NO BOOL */

   public static final boolean Equals( Object in1, Object in2 ) {
      if ( in1 instanceof boolean[] ) return Arrays.equals( (boolean[]) in1, (boolean[]) in2 );
      else {
	 if ( ( (Object[]) in1 ).length != ( (Object[]) in2 ).length ) return false;
	 for ( int i = 0, s = ( (Object[]) in1 ).length; i < s; i ++ ) 
	    if ( ! Equals( ( (Object[]) in1 )[ i ], ( (Object[]) in2 )[ i ] ) ) return false; 
	 return true; } }

   public static int NEqualTo( Object in, boolean val ) { 
      int out = 0; 
      if ( in instanceof boolean[] ) { boolean[] inn = (boolean[]) in; 
      for ( int i = 0, s = inn.length; i < s; i ++ ) if ( inn[ i ] == val ) out ++; }
      else { for ( int i = 0, s = ( (Object[]) in ).length; i < s; i ++ ) 
	 out += NEqualTo( ( (Object[]) in )[ i ], val ); }
      return out; }

   public static int NNotEqualTo( Object in, boolean val ) { 
      return NElem( in ) - NEqualTo( in, val ); }

   public static boolean AnyEqualTo( Object in, boolean val ) { 
      if ( in instanceof boolean[] ) { boolean[] inn = (boolean[]) in; 
      for ( int i = 0, s = inn.length; i < s; i ++ ) if ( inn[ i ] == val ) return true; }
      else { for ( int i = 0, s = ( (Object[]) in ).length; i < s; i ++ ) 
	 if ( AnyEqualTo( ( (Object[]) in )[ i ], val ) ) return true; }
      return false; }   

   public static boolean AnyNotEqualTo( Object in, boolean val ) { 
      if ( in instanceof boolean[] ) { boolean[] inn = (boolean[]) in; 
      for ( int i = 0, s = inn.length; i < s; i ++ ) if ( inn[ i ] != val ) return true; }
      else { for ( int i = 0, s = ( (Object[]) in ).length; i < s; i ++ ) 
	 if ( AnyNotEqualTo( ( (Object[]) in )[ i ], val ) ) return true; }
      return false; }   

   public static final void Shift( boolean[] in, int amt, boolean spare ) {
      if ( amt > 0 ) { System.arraycopy( in, 0, in, amt, in.length - amt ); Arrays.fill( in, 0, amt, spare ); }
      else if ( amt < 0 ) { System.arraycopy( in, -amt, in, 0, in.length + amt ); Arrays.fill( in, in.length+amt, in.length, spare ); } }
   public static final void Shift( boolean[][] in, int amt, boolean spare ) {
      for ( int i = 0, s = in.length; i < s; i ++ ) Shift( in[ i ], amt, spare ); }

   public static final void Mask( boolean[] in, boolean[] flags, boolean val ) {
      for ( int i = 0, s = in.length; i < s; i ++ ) if ( flags[ i ] == val ) in[ i ] = false; }
   public static final int IndexOf( boolean[] in, boolean val ) {
      for ( int i = 0, s = in.length; i < s; i ++ ) if ( in[ i ] == val ) return i; return -1; }
   public static final boolean Contains( boolean[] in, boolean val ) {
      return ( IndexOf( in, val ) != -1 ); }
   public static final int LastIndexOf( boolean[] in, boolean val ) {
      for ( int i = in.length - 1; i >= 0; i -- ) if ( in[ i ] == val ) return i; return -1; }
   public static final boolean[] Index( boolean in[], int index[] ) {
      boolean out[] = New( in.length ); 
      for ( int i = 0, s = in.length; i < s; i ++ ) out[ i ] = in[ index[ i ] ]; return out; }

   public static final boolean[] Reverse( boolean[] in ) {
      for ( int i = 0, s = in.length; i < s/2; i ++ ) {
	 boolean temp = in[ i ]; in[ i ] = in[ s-i-1 ]; in[ s-i-1 ] = temp; }
      return in; }

   /* NO BOOL
   public static final void Sort( boolean[] in ) { Arrays.sort( in ); }
   public static final int[] Indexx( boolean[] in, int[] out ) {
      double darr[] = Convert( in );
      try { return djr.nr.Indexx.indexx( darr, out ); } catch( djr.nr.NRException e ) { }; 
      return null; }
   public static final int[] IRank( boolean[] in, int[] out ) {
      double darr[] = Convert( in );
      try { return djr.nr.Indexx.irank( darr, out ); } catch( djr.nr.NRException e ) { }; 
      return null; }

   public static final int Search( boolean[] in, boolean val, boolean sort ) {
      if ( sort ) Sort( in ); return Arrays.binarySearch( in, val ); }
   public static final int Search( boolean[] in, boolean val ) {
      return Search( in, val, true ); }

   public static final boolean[] Uniq( boolean[] in ) {
      BoolVector v = new BoolVector();
      for ( int i = 0, s = in.length; i < s; i ++ ) 
	 if ( ! v.contains( in[ i ] ) ) v.addElement( true );
      return v.data();
   }
   END NO BOOL */

   public static final int NElem( Object in ) { 
      if ( in instanceof boolean[] ) return ( (boolean[]) in ).length; 
      else { int sum = 0; for ( int i = 0, s = ( (Object[]) in ).length; i < s; i ++ ) 
	 sum += NElem( ( (Object[]) in )[ i ] ); return sum; } }
      
   /* NO BOOL
   public static final double Mean( Object in ) {
      return ( (double) Sum( in ) / (double) NElem( in ) ); }

   public static final double Mean( Object in, int i1, int i2 ) {
      return ( (double) Sum( in, i1, i2 ) / (double) ( i2 - i1 + 1 ) ); }

   public static final double Median( boolean[] in ) {
      final int nx = in.length; boolean temp[] = New( nx );
      Copy( temp, in ); Sort( temp ); double out = (double) temp[ (nx+1)/2 ];
      if ( nx % 2 != 1 ) out = ( (double) temp[ nx/2 ] + (double) temp[ nx/2+1 ] ) / 2.0;
      return out; }

   public static final double Variance( boolean[] in, double mean ) {
      double sum = 0.0; int nelem = 0;
      for ( int i = 0, s = in.length; i < s; i ++ ) 
	 sum += ( in[ i ] - mean ) * ( in[ i ] - mean );
      return sum / ( (double) in.length - 1.0 ); }

   public static final double Stddev( boolean[] in, double mean ) {
      return Math.sqrt( Variance( in, mean ) ); }

   public static final double Stddev( boolean[] in ) {
      return Stddev( in, Mean( in ) ); }

   public static final void Stats( boolean[] in, double[] out ) {
      out[ 0 ] = Mean( in ); out[ 1 ] = Median( in ); out[ 2 ] = Stddev( in ); }
   public static final double[] Stats( boolean[] in ) {
      double out[] = new double[ 3 ]; Stats( in, out ); return out; }

   public static void PrintStats( boolean[] data ) {
      System.out.println      ( "NELEM   = " +      NElem(data) );
      Format.print( System.out, "MEAN    = %.3f\n", Mean(data) );
      Format.print( System.out, "MEDIAN  = %.3f\n", Median(data) );
      Format.print( System.out, "MODE    = %.3f\n", Mode(data) );
      Format.print( System.out, "STDDEV  = %.3f\n", Stddev(data) );
      double stddevs = (double) (Max(data) - Mode(data)) / (double) Stddev(data);
      Format.print( System.out, "STDDEVS = %.3f\n", stddevs );
      Format.print( System.out, "ERF     = %.5f\n", Probability.errorFunction( stddevs ) );
      Format.print( System.out, "NORMAL  = %.5f\n", Probability.normal( stddevs ) );
   }

   public static final double[][] Histogram( boolean[] in, int nbins ) {
      double min = (double) Min( in ), max = (double) Max( in );
      double binsize = ( max - min ) / ( (double) nbins ), bs2 = binsize / 2.0;
      double[][] out = DoubleUtils.New( 2, nbins + 1 ); int ind = 0;
      for ( double i = min; i <= max; i += binsize ) {
	 out[ 0 ][ ind ] = i;
	 for ( int j = 0, size = in.length; j < size; j ++ )
	    if ( in[ j ] != 0 && in[ j ] >= i - bs2 && in[ j ] < i + bs2 ) out[ 1 ][ ind ] ++;
	 if ( ind ++ >= nbins ) break; }
      return out; }

   public static final double[] Histogram( boolean[] in, boolean[] bins ) {
      double out[] = DoubleUtils.New( bins.length );
      for ( int i = 0, s = bins.length; i < s; i ++ ) {
	 double bin = (double) bins[ i ];
	 double lower = (double) ( i == 0 ? -Double.MAX_VALUE : ( bins[ i ] + bins[ i-1 ] ) / 2 );
	 double upper = (double) ( i == s-1 ? Double.MAX_VALUE : ( bins[ i ] + bins[ i+1 ] ) / 2 );
	 for ( int j = 0, size = in.length; j < size; j ++ )
	    if ( in[ j ] != 0 && in[ j ] >= lower && in[ j ] < upper ) out[ i ] ++;
      }
      return out;
   }

   public static final double Mode( boolean[] in ) {
      int nbins = Math.min( (int) ( in.length / 20.0 ), 20 );
      double[][] hist = Histogram( in, nbins );
      return hist[ 0 ][ DoubleUtils.WhereMax( hist[ 1 ] ) ]; }

   public static final boolean Max( boolean x1, boolean x2 ) {
      return ( x1 > x2 ? x1 : x2 ); }
   public static final boolean Max( boolean x1, boolean x2, boolean x3 ) {
      return Max( x1, Max( x2, x3 ) ); }
   public static final boolean Max( boolean x1, boolean x2, boolean x3, boolean x4 ) {
      return Max( Max( x1, x2 ), Max( x3, x4 ) ); }
   public static final boolean Min( boolean x1, boolean x2 ) {
      return ( x1 < x2 ? x1 : x2 ); }
   public static final boolean Approx( boolean x1, boolean x2 ) {
      return ( Math.abs( x1 - x2 ) < 1e-8 ); }

   public static final int Sample( boolean[] x ) {
      return Sample( x, x.length, 500 ); }
   public static final int Sample( boolean[] x, int len ) {
      return Sample( x, len, 500 ); }

   // Assumes x[] is normalized so its max. value is 1.
   public static final int Sample( boolean[] x, int len, int maxCount ) {
      int r1 = -1, count = 0; double r2, xx; len --;
      while( count ++ < maxCount ) {
	 r1 = IntUtils.RandChoose( 0, len );
	 xx = x[ r1 ];
	 if ( xx <= 0.0 ) { r1 = -1; continue; }
	 r2 = DoubleUtils.Random(); if ( r2 <= xx ) break;
	 r1 = -1;
      }
      return r1; }

   public static final double[] Convolve( boolean x[], boolean kern[], double out[] ) {
      if ( out == null || out.length < x.length ) out = DoubleUtils.New( x.length );
      else DoubleUtils.Zero( out );
      UnitNorm( kern );
      int ksize = kern.length, k2 = ksize / 2;
      for ( int i = 0, s = x.length, min = i-k2, max = i+k2; 
	    i < s; i ++, min ++, max ++ ) {
	 for ( int j = min, k = 0; j <= max; j ++, k ++ ) {
	    if ( j < 0 ) continue; else if ( j >= s ) break;
	    out[ i ] += (double) ( x[ j ] * kern[ k ] );
	 }
      }
      return out;
   }

   public static final double logGamma( boolean xx ) {
      return cern.jet.stat.Gamma.logGamma( (double) xx );
   }
   END NO BOOL */

   public static final boolean[] FTokenize( String str, final String tok ) { 
      return FTokenize( MyUtils.Tokenize( str, tok ) ); }

   public static final boolean[] FTokenize( String str[] ) { 
      return FTokenize( str, 0, str.length ); }

   public static final boolean[] FTokenize( String str[], int start, int end ) { 
      boolean[] out = New( end - start );
      for ( int i = start; i < end; i ++ ) {
   /* NO BOOL
	 try { out[ i - start ] = (boolean) Boolean.parseBool( str[ i ] ) ; } 
   END NO BOOL */
   // BOOL ONLY
	 try { out[ i - start ] = Boolean.valueOf( str[ i ] ).booleanValue() ; } 
   // END BOOL ONLY */
	 catch( Exception e ) { out[ i ] = false; } } return out; }

   public static final boolean[] FromVector( Vector v ) {
      if ( v == null || v.size() <= 0 ) return null;
      boolean out[] = New( v.size() );
      try {
	 for ( int i = 0, s = v.size(); i < s; i ++ )
	    out[ i ] = ( (Boolean) v.elementAt( i ) ).booleanValue();
      } catch( Exception e ) { };
      return out;
   }

   public static final Vector ToVector( boolean[] in ) {
      if ( in == null ) return null;
      Vector v = new Vector();
      for ( int i = 0, s = in.length; i < s; i ++ ) v.addElement( new Boolean ( in[ i ] ) );
      return v;
   }
}
