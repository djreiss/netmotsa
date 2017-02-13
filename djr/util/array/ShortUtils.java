package djr.util.array;

import java.util.*;
import corejava.*;
import cern.jet.stat.*;
import cern.jet.random.Normal;
//import cern.jet.random.engine.MersenneTwister; //DRand;
import edu.cornell.lassp.houle.RngPack.Ranecu;
import djr.util.MyUtils;

/**
 * Template class <code>ShortUtils</code>.
 *
 * This is a "template" class and should not be directly compiled.
 * Commands to create real file before compiling (e.g. for ints) are:
 * cp Utils.java IntUtils.java
 * rpl Short Integer IntUtils.java ; rpl Short Int IntUtils.java ; rpl short int IntUtils.java
 * javac IntUtils.java
 *
 * @author <a href="mailto:reiss@uw.edu">David Reiss</a>
 * @version 1.9978 (Fri Nov 07 05:56:26 PST 2003)
 */
public class ShortUtils {
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
   // NO BOOL
   public static final short Random() {
      return (short) rand.raw(); }
   public static final long GetSeed() {
      return DoubleUtils.rand.getSeed(); }
   public static final void Random( double arr[] ) { 
      DoubleUtils.rand.raw( arr, arr.length ); }
   public static final short Random( int max ) { 
      return (short) ( DoubleUtils.Random() * max ); }
   public static final short Random( double min, double max ) { 
      return (short) ( min + ( max - min ) * DoubleUtils.Random() ); }
   public static final short Random( int min, int max ) { 
      return Random( (double) min, (double) max ); }
   public static final short RandChoose( int lo, int hi ) {
      return (short) ( (long) lo + (long) ( ( 1L + (long) hi - (long) lo ) * 
					   DoubleUtils.Random() ) ); }
   public static final short RandChoose( int hi ) {
      return RandChoose( 0, hi ); }
   // END NO BOOL
   public static final short RandChoose( short[] in ) {
      return in[ IntUtils.RandChoose( 0, in.length - 1 ) ]; }

   /* BOOL ONLY
   public static final boolean Random() {
      return rand.raw() > 0.5 ? false : true; }
   // END BOOL ONLY */

   // NO BOOL
   public static final short Normal( double mean, double stddev ) { 
      return (short) Normal.staticNextDouble( mean, stddev ); }
   public static final short Normal() { 
      return Normal( 0.0, 1.0 ); }

   public static final double Log( short in ) { return in <= 0 ? -999 : Math.log( in ); }
   public static final double Log2( short in ) { return DoubleUtils.Log( in ) / LOG2; }
   public static final double Log10( short in ) { return DoubleUtils.Log( in ) / LOG10; }
   public static final double Sign( short in ) { return in < 0 ? -1.0 : 1.0; }
   // END NO BOOL

   public static final short[] New( int nx ) { 
      return new short[ nx ]; }
   public static final short[][] New( int nx, int ny ) { 
      return new short[ nx ][ ny ]; }
   public static final short[][][] New( int nx, int ny, int nz ) { 
      return new short[ nx ][ ny ][ nz ]; }
   public static final short[][][][] New( int nx, int ny, int nz, int na ) { 
      return new short[ nx ][ ny ][ nz ][ na ]; }

   public static final short[] Resize( short[] in, int newx ) {
      short[] out = New( newx ); Copy( out, in ); return out; }
   public static final short[][] Resize( short[][] in, int newx, int newy ) {
      short[][] out = New( newx, newy ); Copy( out, in ); return out; }
   public static final short[][][] Resize( short[][][] in, int newx, int newy, int newz ) {
      short[][][] out = New( newx, newy, newz ); Copy( out, in ); return out; }
   public static final short[][][][] Resize( short[][][][] in, int newx, int newy, int newz, int newa ) {
      short[][][][] out = New( newx, newy, newz, newa ); Copy( out, in ); return out; }

   public static final short[] New( short[] in ) {
      if ( in == null ) return null;
      short[] out = New( in.length ); Copy( out, in ); return out; }
   public static final short[][] New( short[][] in ) {
      short[][] out = new short[ in.length ][];
      for ( int i = 0, s = in.length; i < s; i ++ ) out[ i ] = New( in[ i ] ); return out; }
   public static final short[][][] New( short[][][] in ) {
      short[][][] out = new short[ in.length ][][];
      for ( int i = 0, s = in.length; i < s; i ++ ) out[ i ] = New( in[ i ] ); return out; }

   public static final Object NewObj( Object obj ) {
      if ( obj instanceof short[] ) return New( (short[]) obj );
      else if ( obj instanceof short[][] ) return New( (short[][]) obj );
      else if ( obj instanceof short[][][] ) return New( (short[][][]) obj );
      return null; }

   public static final short[] SubArr( short[] in, int start, int end ) {
      short[] out = New( end - start + 1 );
      System.arraycopy( in, start, out, 0, end - start + 1 );
      return out; }

   // NO BOOL
   public static final short[] Sequence( int length ) {
      short[] out = New( length ); for ( int i = 0; i < length; i ++ ) out[ i ] = (short) i;
      return out; }
   public static final double[] Convert( Object in ) {
      if ( in instanceof double[] ) return (double[]) in;
      short iin[] = (short[]) in;
      double out[] = new double[ iin.length ];
      for ( int i = 0, size = iin.length; i < size; i ++ ) out[ i ] = (double) iin[ i ];
      return out; }
   // END NO BOOL

   public static final void Print( Object in ) {
      if ( in instanceof short[] ) Print( (short[]) in, 0, ( (short[]) in ).length ); 
      else for ( int i = 0, size = ( (Object[]) in ).length; i < size; i ++ ) 
	 Print( ( (Object[]) in )[ i ] ); }
   public static final void Print( String start, Object in ) {
      System.out.print( start ); Print( in ); }

   public static final void Print( short[] in, int start, int end ) {
      if ( in == null ) System.out.println( "NULL" );
      for ( int i = start; i < end; i ++ ) System.out.print( in[ i ] + " " );
      System.out.println(); }

   // NO BOOL
   public static final void Printf( String fmt, short in ) {
      if ( Double.isInfinite( (double) in ) ) in = - Short.MAX_VALUE;
      Format.print( System.out, fmt, in ); };
   public static final String SPrintf( String fmt, short in ) {
      if ( Double.isInfinite( (double) in ) ) in = - Short.MAX_VALUE;
      return ( new Format( fmt ) ).form( in ); }

   public static final void PrintCols( String fmt, short in1[], short in2[] ) {
      PrintCols( fmt, new short[][]{ in1, in2 } ); }
   public static final void PrintCols( String fmt, short in[][], boolean lines ) {
      for ( int i = 0, size = in[ 0 ].length; i < size; i ++ ) {
	 if ( lines ) System.out.print( i + " " );
	 for ( int j = 0, ss = in.length; j < ss; j ++ ) 
	    System.out.print( SPrintf( fmt, in[ j ][ i ] ) );
	 System.out.println(); } }
   public static final void PrintCols( String fmt, short in[][] ) {
      PrintCols( fmt, in, false ); }

   public static final void Printf( String fmt, Object in ) {
      if ( in == null ) System.out.println( "NULL" );
      if ( in instanceof short[] ) { for ( int i = 0, size = ( (short[]) in ).length; i < size; i ++ )
	 Printf( fmt, ( (short[]) in )[ i ] ); }
      else { for ( int i = 0, size = ( (Object[]) in ).length; i < size; i ++ )
	 Printf( fmt, ( (Object[]) in )[ i ] ); }
      System.out.println(); }
   public static final void Printf( String start, String fmt, Object in ) {
      System.out.print( start ); Printf( fmt, in ); }

   public static final String SPrintf( String fmt, Object in ) {
      if ( in == null ) return "NULL";
      String out = "";
      if ( in instanceof short[] ) { for ( int i = 0, size = ( (short[]) in ).length; i < size; i ++ )
	 out += SPrintf( fmt, ( (short[]) in )[ i ] ); }
      else { for ( int i = 0, size = ( (Object[]) in ).length; i < size; i ++ )
	 out += SPrintf( fmt, ( (Object[]) in )[ i ] ); }
      return out; }

   public static final void PrintT( String fmt, short[][] in ) {
      for ( int y = 0, size = in[0].length; y < size; y ++ ) {
	 for ( int x = 0, s = in.length; x < s; x ++ ) Format.print( System.out, fmt, in[ x ][ y ] );
	 System.out.println(); } }
   public static final void PrintT( short[][] in ) {
      for ( int y = 0, size = in[0].length; y < size; y ++ ) {
	 for ( int x = 0, s = in.length; x < s; x ++ ) System.out.print( in[ x ][ y ] + " " );
	 System.out.println(); } }
   // END NO BOOL

   public static final short[][] Transpose( short[][] in ) {
      short out[][] = New( in[ 0 ].length, in.length );
      for ( int i = 0; i < in[ 0 ].length; i ++ ) for ( int j = 0; j < in.length; j ++ ) 
	 out[ i ][ j ] = in[ j ][ i ]; return out; }

   public static final void Copy( Object dest, Object src ) {
      if ( src == null || dest == null ) return;
      if ( src instanceof short[] ) System.arraycopy( (short[]) src, 0, (short[]) dest, 0, 
			IntUtils.Min( ( (short[]) src ).length, ( (short[]) dest ).length ) );
      else for ( int i = 0, s = IntUtils.Min( ( (Object[]) src ).length, 
	     ( (Object[]) dest ).length ); i < s; i ++ ) 
	 Copy( ( (Object[]) dest )[ i ], ( (Object[]) src )[ i ] );
   }

   public static final short[] Flatten( short[][] in, short[] out ) {
      if ( out == null || out.length < NElem( in ) ) out = new short[ NElem( in ) ]; 
      int ind = 0; for ( int i = 0, s = in.length; i < s; i ++ ) {
	 System.arraycopy( in[ i ], 0, out, ind, in[ i ].length ); ind += in[ i ].length; }
      return out; }
   public static final short[] Flatten( short[][] in ) {
      return Flatten( in, null ); }

   public static final short[] Flatten( short[][][] in, short[] out ) {
      if ( out == null || out.length < NElem( in ) ) out = new short[ NElem( in ) ]; 
      int ind = 0; for ( int i = 0, s = in.length; i < s; i ++ ) {
	 short temp[] = Flatten( in[ i ] );
	 System.arraycopy( temp, 0, out, ind, temp.length ); ind += temp.length; }
      return out; }
   public static final short[] Flatten( short[][][] in ) {
      return Flatten( in, null ); }

   private static short[][] UnflattenPriv( short in[], int start, int nx, int ny, short out[][] ) {
      if ( out == null || out.length < nx ) out = new short[ nx ][]; int ind = start;
      for ( int i = 0; i < nx; i ++ ) {
	 if ( out[ i ] == null || out[ i ].length < ny ) out[ i ] = new short[ ny ];
	 System.arraycopy( in, ind, out[ i ], 0, ny ); ind += ny;
      } return out; }
   public static short[][] Unflatten( short in[], int nx, int ny, short out[][] ) {
      return UnflattenPriv( in, 0, nx, ny, out ); }
   public static short[][] Unflatten( short in[], int nx, int ny ) {
      return UnflattenPriv( in, 0, nx, ny, null ); }

   public static short[][][] Unflatten( short in[], int nx, int ny, int nz, short out[][][] ) {
      if ( out == null || out.length < nx ) out = new short[ nx ][][]; int ind = 0, step = ny *nz;
      for ( int i = 0; i < nx; i ++ ) { out[ i ] = UnflattenPriv( in, ind, ny, nz, out[ i ] );
      ind += step; } return out; }
   public static short[][][] Unflatten( short in[], int nx, int ny, int nz ) {
      return Unflatten( in, nx, ny, nz, null ); }

   public static final void Zero( Object in, int start, int end ) {
      if ( in == null ) return;
      if ( in instanceof short[] ) Arrays.fill( (short[]) in, start, end, (short) 0 );
      else for ( int i = 0, s = ( (Object[]) in ).length; i < s; i ++ ) 
	 Zero( ( (Object[]) in )[ i ], start, end ); }
   public static final void Zero( Object in ) {
      if ( in == null ) return;
      if ( in instanceof short[] ) Arrays.fill( (short[]) in, (short) 0 );
      else for ( int i = 0, s = ( (Object[]) in ).length; i < s; i ++ ) 
	 Zero( ( (Object[]) in )[ i ] ); }
   public static final void Zero( short in[], int inds[] ) {
      if ( in == null ) return;
      else for ( int i = 0, s = inds.length; i < s; i ++ ) in[ inds[ i ] ] = 0;
   }

   public static final void Set( Object in, int from, int to, short val ) {
      if ( in instanceof short[] ) Arrays.fill( (short[]) in, from, to, val );
      else for ( int i = 0, s = ( (Object[]) in ).length; i < s; i ++ ) 
	 Set( ( (Object[]) in )[ i ], from, to, val ); }
   public static final void Set( Object in, short val ) {
      if ( in == null ) return;
      if ( in instanceof short[] ) Arrays.fill( (short[]) in, val );
      else for ( int i = 0, s = ( (Object[]) in ).length; i < s; i ++ ) 
	 Set( ( (Object[]) in )[ i ], val ); }
   public static final void Set( short in[], int inds[], short val ) {
      if ( in == null ) return;
      else for ( int i = 0, s = inds.length; i < s; i ++ ) in[ inds[ i ] ] = val;
   }

   // NO BOOL
   public static final short Sum( Object in, int min, int max ) {
      if ( in == null ) return 0;
   // END NO BOOL
   /* BOOL ONLY
   public static final int Sum( Object in, int min, int max ) {
   // END BOOL ONLY */
      if ( in instanceof short[] ) { short[] inn = (short[]) in; 
	 if ( max >= inn.length ) max = inn.length - 1;
   // NO BOOL
	 short out = 0; for ( int i = min; i <= max; i ++ ) out += inn[ i ]; return out; 
   // END NO BOOL
   /* BOOL ONLY
         int out = 0; for ( int i = min; i <= max; i ++ ) out += inn[ i ] ? 1 : 0; 
	 return out; 
   // END BOOL ONLY */
      } else {
   // NO BOOL
	 short out = 0; 
   // END NO BOOL
   /* BOOL ONLY
         int out = 0; 
   // END BOOL ONLY */
	 for ( int i = 0, s = ( (Object[]) in ).length; i < s; i ++ ) 
	    out += Sum( ( (Object[]) in )[ i ], min, max ); return out;
      } }

   // NO BOOL
   public static final short Sum( Object in ) {
   // END NO BOOL
   /* BOOL ONLY
   public static final int Sum( Object in ) {
   // END BOOL ONLY */
      if ( in instanceof short[] ) return Sum( (short[]) in, 0, ( (short[]) in ).length ); 
      else {
   // NO BOOL
	 short out = 0; 
   // END NO BOOL
   /* BOOL ONLY
         int out = 0; 
   // END BOOL ONLY */
	 for ( int i = 0, s = ( (Object[]) in).length; i < s; i ++ ) 
	    out += Sum( ( (Object[]) in )[ i ] ); return out; }
   }

   // NO BOOL
   public static final short SumSquared( Object in ) {
      if ( in == null ) return 0;
      if ( in instanceof short[] ) { short[] inn = (short[]) in; short out = 0; 
      for ( int i = 0, s = inn.length; i < s; i ++ ) out += inn[ i ] * inn[ i ]; return out;
      } else {
	 short out = 0; for ( int i = 0, s = ( (Object[]) in).length; i < s; i ++ ) 
	    out += SumSquared( ( (Object[]) in )[ i ] ); return out; }
   }

   public static final short SumSquared2( Object in ) {
      if ( in == null ) return 0;
      if ( in instanceof short[] ) { short[] inn = (short[]) in; short out = 0; 
      for ( int i = 0, s = inn.length; i < s; i ++ ) out += Sign( inn[ i ] ) * inn[ i ] * inn[ i ]; return out;
      } else {
	 short out = 0; for ( int i = 0, s = ( (Object[]) in).length; i < s; i ++ ) 
	    out += SumSquared2( ( (Object[]) in )[ i ] ); return out; }
   }

   public static final short Prod( Object in ) {
      if ( in == null ) return 0;
      if ( in instanceof short[] ) { short[] inn = (short[]) in; short out = inn[ 0 ]; for ( int i = 1, s = inn.length; i < s; i ++ ) out *= inn[ i ]; return out;
      } else {
	 short out = Prod( ( (Object[]) in )[ 0 ] ); for ( int i = 1, s = ( (Object[]) in).length; i < s; i ++ ) 
	    out *= Prod( ( (Object[]) in )[ i ] ); return out; }
   }

   public static final short Max( Object in ) { short out = - Short.MAX_VALUE;
      if ( in == null ) return 0;
      if ( in instanceof short[] ) { short[] inn = (short[]) in; 
      for ( int i = 0, s = inn.length; i < s; i ++ ) out = Max( out, inn[ i ] ); return out; 
      } else { for ( int i = 0, s = ( (Object[]) in ).length; i < s; i ++ ) 
	 out = Max( out, Max( ( (Object[]) in )[ i ] ) ); return out; } }
   
   public static final short Min( Object in ) { short out = Short.MAX_VALUE;
      if ( in == null ) return 0;
      if ( in instanceof short[] ) { short[] inn = (short[]) in; 
      for ( int i = 0, s = inn.length; i < s; i ++ ) out = Min( out, inn[ i ] ); return out; 
      } else { for ( int i = 0, s = ( (Object[]) in ).length; i < s; i ++ ) 
	 out = Min( out, Min( ( (Object[]) in )[ i ] ) ); return out; } }

   public static final short[] MinMax( Object in, short inout[] ) { 
      if ( in == null ) return null;
      if ( inout == null ) { inout = New( 2 ); 
      inout[ 0 ] = Short.MAX_VALUE; inout[ 1 ] = - Short.MAX_VALUE; }
      if ( in instanceof short[] ) { short[] inn = (short[]) in; 
      for ( int i = 0, s = inn.length; i < s; i ++ ) {
	 inout[ 0 ] = Min( inout[ 0 ], inn[ i ] ); 
	 inout[ 1 ] = Max( inout[ 1 ], inn[ i ] ); } return inout;
      } else { for ( int i = 0, s = ( (Object[]) in ).length; i < s; i ++ ) 
	 inout = MinMax( ( (Object[]) in )[ i ], inout ); return inout; } }
   public static final short[] MinMax( Object in ) { 
      return MinMax( in, null ); }

   public static final void Round( Object in ) { 
      if ( in == null ) return;
      if ( in instanceof short[] ) { short[] inn = (short[]) in; 
      for ( int i = 0, s = inn.length; i < s; i ++ ) inn[ i ] = (short) (int) inn[ i ]; }
      else { for ( int i = 0, s = ( (Object[]) in ).length; i < s; i ++ ) 
	 Round( ( (Object[]) in )[ i ] ); }
   }

   public static final void Floor( Object in, short val ) { 
      if ( in == null ) return;
      if ( in instanceof short[] ) { short[] inn = (short[]) in; 
      for ( int i = 0, s = inn.length; i < s; i ++ ) if ( inn[ i ] < val ) inn[ i ] = val; }
      else { for ( int i = 0, s = ( (Object[]) in ).length; i < s; i ++ ) 
	 Floor( ( (Object[]) in )[ i ], val ); }
   }

   public static final void Ceil( Object in, short val ) { 
      if ( in == null ) return;
      if ( in instanceof short[] ) { short[] inn = (short[]) in; 
      for ( int i = 0, s = inn.length; i < s; i ++ ) if ( inn[ i ] > val ) inn[ i ] = val; }
      else { for ( int i = 0, s = ( (Object[]) in ).length; i < s; i ++ ) 
	 Ceil( ( (Object[]) in )[ i ], val ); }
   }
   
   public static final void Mult( Object in, short val ) { 
      if ( in == null ) return;
      if ( in instanceof short[] ) { short[] inn = (short[]) in; 
      for ( int i = 0, s = inn.length; i < s; i ++ ) inn[ i ] *= val; }
      else { for ( int i = 0, s = ( (Object[]) in ).length; i < s; i ++ ) 
	 Mult( ( (Object[]) in )[ i ], val ); }
   }

   public static final void Mult( Object in, Object in2 ) { 
      if ( in == null || in2 == null ) return;
      if ( in instanceof short[] ) { short[] inn = (short[]) in; 
      for ( int i = 0, s = inn.length; i < s; i ++ ) inn[ i ] *= ( (short[]) in2 )[ i ]; }
      else { for ( int i = 0, s = ( (Object[]) in ).length; i < s; i ++ ) 
	 Mult( ( (Object[]) in )[ i ], ( (Object[]) in2 )[ i ] ); }
   }

   public static final void Divide( Object in, short val ) { 
      if ( in == null ) return;
      if ( val == 0 ) return;
      if ( in instanceof short[] ) { short[] inn = (short[]) in; 
      for ( int i = 0, s = inn.length; i < s; i ++ ) inn[ i ] /= val; }
      else { for ( int i = 0, s = ( (Object[]) in ).length; i < s; i ++ ) 
	 Divide( ( (Object[]) in )[ i ], val ); }
   }

   public static final void Divide( short val, Object in ) { 
      if ( in == null ) return;
      if ( val == 0 ) return;
      if ( in instanceof short[] ) { short[] inn = (short[]) in; 
      for ( int i = 0, s = inn.length; i < s; i ++ ) inn[ i ] = (short) ( val / inn[ i ] ); }
      else { for ( int i = 0, s = ( (Object[]) in ).length; i < s; i ++ ) 
	 Divide( ( (Object[]) in )[ i ], val ); }
   }

   public static final void Divide( Object in, Object in2 ) { 
      if ( in == null || in2 == null ) return;
      if ( in instanceof short[] ) { short[] inn = (short[]) in; 
      for ( int i = 0, s = inn.length; i < s; i ++ ) inn[ i ] /= ( (short[]) in2 )[ i ]; }
      else { for ( int i = 0, s = ( (Object[]) in ).length; i < s; i ++ ) 
	 Divide( ( (Object[]) in )[ i ], ( (Object[]) in2 )[ i ] ); }
   }

   public static final void Add( Object in, short val ) { 
      if ( in == null ) return;
      if ( in instanceof short[] ) { short[] inn = (short[]) in; 
      for ( int i = 0, s = inn.length; i < s; i ++ ) inn[ i ] += val; }
      else { for ( int i = 0, s = ( (Object[]) in ).length; i < s; i ++ ) 
	 Add( ( (Object[]) in )[ i ], val ); }
   }

   public static final void Add( Object in, Object in2 ) { 
      if ( in == null || in2 == null ) return;
      if ( in instanceof short[] ) { short[] inn = (short[]) in; 
      for ( int i = 0, s = inn.length; i < s; i ++ ) inn[ i ] += ( (short[]) in2 )[ i ]; }
      else { for ( int i = 0, s = IntUtils.Min( ( (Object[]) in ).length, 
						( (Object[]) in2 ).length ); i < s; i ++ ) 
	 Add( ( (Object[]) in )[ i ], ( (Object[]) in2 )[ i ] ); }
   }

   public static final void Sub( Object in, short val ) { 
      Add( in, (short) -val ); }

   public static final void Sub( Object in, Object in2 ) { 
      if ( in == null || in2 == null ) return;
      if ( in instanceof short[] ) { short[] inn = (short[]) in; 
      for ( int i = 0, s = inn.length; i < s; i ++ ) inn[ i ] -= ( (short[]) in2 )[ i ]; }
      else { for ( int i = 0, s = ( (Object[]) in ).length; i < s; i ++ ) 
	 Sub( ( (Object[]) in )[ i ], ( (Object[]) in2 )[ i ] ); }
   }

   public static final void Invert( Object in ) {
      if ( in == null ) return;
      if ( in instanceof short[] ) { short[] inn = (short[]) in; 
      for ( int i = 0, s = inn.length; i < s; i ++ ) 
	 inn[ i ] = (short) ( 1.0 / (double) inn[ i ] ); }
      else { for ( int i = 0, s = ( (Object[]) in ).length; i < s; i ++ ) 
	 Invert( ( (Object[]) in )[ i ] ); }
   }

   // END NO BOOL

   public static final void Random( Object in ) {
      if ( in == null ) return;
      if ( in instanceof short[] ) { short[] inn = (short[]) in; 
      for ( int i = 0, s = inn.length; i < s; i ++ ) inn[ i ] = Random(); }
      else { for ( int i = 0, s = ( (Object[]) in ).length; i < s; i ++ ) 
	 Random( ( (Object[]) in )[ i ] ); }
   }

   // NO BOOL
   public static final void Normal( Object in, double mean, double stddev ) {
      if ( in == null ) return;
      if ( nrand == null ) nrand = new Normal( mean, stddev, DoubleUtils.rand );
      else nrand.setState( mean, stddev );
      if ( in instanceof short[] ) { short[] inn = (short[]) in; 
      for ( int i = 0, s = inn.length; i < s; i ++ ) inn[ i ] = (short) nrand.nextDouble(); }
      else { for ( int i = 0, s = ( (Object[]) in ).length; i < s; i ++ ) 
	 Normal( ( (Object[]) in )[ i ], mean, stddev ); }
   }

   public static final void NormalAdd( Object in, double mean, double stddev ) {
      if ( in == null ) return;
      if ( nrand == null ) nrand = new Normal( mean, stddev, DoubleUtils.rand );
      else nrand.setState( mean, stddev );
      if ( in instanceof short[] ) { short[] inn = (short[]) in; 
      for ( int i = 0, s = inn.length; i < s; i ++ ) inn[ i ] += (short) nrand.nextDouble(); }
      else { for ( int i = 0, s = ( (Object[]) in ).length; i < s; i ++ ) 
	 NormalAdd( ( (Object[]) in )[ i ], mean, stddev ); }
   }

   public static final void Log10( Object in ) { 
      if ( in == null ) return;
      if ( in instanceof short[] ) { short[] inn = (short[]) in; 
      for ( int i = 0, s = inn.length; i < s; i ++ ) inn[ i ] = (short) Log10( inn[ i ] ); }
      else { for ( int i = 0, s = ( (Object[]) in ).length; i < s; i ++ ) 
	 Log10( ( (Object[]) in )[ i ] ); }
   }

   public static final void Log2( Object in ) { 
      if ( in == null ) return;
      if ( in instanceof short[] ) { short[] inn = (short[]) in; 
      for ( int i = 0, s = inn.length; i < s; i ++ ) inn[ i ] = (short) Log2( inn[ i ] ); }
      else { for ( int i = 0, s = ( (Object[]) in ).length; i < s; i ++ ) 
	 Log2( ( (Object[]) in )[ i ] ); }
   }

   public static final void Log( Object in ) { 
      if ( in == null ) return;
      if ( in instanceof short[] ) { short[] inn = (short[]) in; 
      for ( int i = 0, s = inn.length; i < s; i ++ ) inn[ i ] = (short) Log( inn[ i ] ); }
      else { for ( int i = 0, s = ( (Object[]) in ).length; i < s; i ++ ) 
	 Log( ( (Object[]) in )[ i ] ); }
   }

   public static final void Pow( short val, Object in ) { 
      if ( in == null ) return;
      if ( in instanceof short[] ) { short[] inn = (short[]) in; 
      for ( int i = 0, s = inn.length; i < s; i ++ ) 
	 inn[ i ] = (short) Math.pow( val, (double) inn[ i ] ); }
      else { for ( int i = 0, s = ( (Object[]) in ).length; i < s; i ++ ) 
	 Pow( ( (Object[]) in )[ i ], val ); }
   }

   public static final void Pow( Object in, short val ) { 
      if ( in == null ) return;
      if ( in instanceof short[] ) { short[] inn = (short[]) in; 
      for ( int i = 0, s = inn.length; i < s; i ++ ) 
	 inn[ i ] = (short) Math.pow( (double) inn[ i ], val ); }
      else { for ( int i = 0, s = ( (Object[]) in ).length; i < s; i ++ ) 
	 Pow( ( (Object[]) in )[ i ], val ); }
   }

   public static final void Pow2( Object in, short val ) { 
      if ( in == null ) return;
      if ( in instanceof short[] ) { short[] inn = (short[]) in; 
      for ( int i = 0, s = inn.length; i < s; i ++ ) 
	 inn[ i ] = (short) ( Sign( inn[ i ] ) * Math.pow( (double) inn[ i ], val ) ); }
      else { for ( int i = 0, s = ( (Object[]) in ).length; i < s; i ++ ) 
	 Pow2( ( (Object[]) in )[ i ], val ); }
   }

   public static final void Abs( Object in ) { 
      if ( in == null ) return;
      if ( in instanceof short[] ) { short[] inn = (short[]) in; 
      for ( int i = 0, s = inn.length; i < s; i ++ ) {
	 short x = inn[ i ]; inn[ i ] = (short) ( x > 0 ? x : -x ); } }
      else { for ( int i = 0, s = ( (Object[]) in ).length; i < s; i ++ ) 
	 Abs( ( (Object[]) in )[ i ] ); }
   }

   public static final void Replace( Object in, short from, short to ) {
      if ( in == null ) return;
      if ( in instanceof short[] ) { short[] inn = (short[]) in; 
      for ( int i = 0, s = inn.length; i < s; i ++ ) 
	 if ( inn[ i ] == from ) inn[ i ] = to; }
      else { for ( int i = 0, s = ( (Object[]) in ).length; i < s; i ++ ) 
	 Replace( ( (Object[]) in )[ i ], from, to ); }      
   }

   public static final int WhereMax( short[] in, int maxind ) { 
      if ( in == null ) return -1;
      short out = - Short.MAX_VALUE ; int ind = -1;
      maxind = Math.min( in.length, maxind+1 );
      for ( int i = 0; i < maxind; i ++ ) 
	 if ( out < in[ i ] ) { out = in[ i ]; ind = i; } return ind; }
   public static final int WhereMax( short[] in ) {
      return WhereMax( in, in.length ); }

   public static final int WhereMin( short[] in, int maxind ) { 
      short out = Short.MAX_VALUE ; int ind = -1;
      maxind = Math.min( in.length, maxind+1 );
      for ( int i = 0; i < maxind; i ++ ) if ( out > in[ i ] ) { out = in[ i ]; ind = i; } 
      return ind; }
   public static final int WhereMin( short[] in ) { 
      return WhereMin( in, in.length ); }

   public static final void MaxNorm( Object in ) { 
      short max = Max( in ); Divide( in, max ); }
   public static final void UnitNorm( Object in ) {
      short sum = Sum( in ); Divide( in, sum ); }

   public static final double Integral( short[] inx, short[] iny ) { // trapezoid rule
      double sum = 0.0;
      for ( int i = 1, s = inx.length; i < s; i ++ )
	 sum += ( inx[ i ] + inx[ i - 1 ] ) / 2.0 * ( iny[ i ] - iny[ i - 1 ] ) ;
      return sum;
   }

   public static final short[] Integrate( short[] in, short[] out ) {
      if ( out == null || out.length < in.length ) out = New( in.length );
      short sum = 0;
      for ( int i = 0, s = in.length; i < s; i ++ ) { sum += in[ i ]; out[ i ] = sum; }
      return out;
   }
   // END NO BOOL

   public static final boolean Equals( Object in1, Object in2 ) {
      if ( in1 instanceof short[] ) return Arrays.equals( (short[]) in1, (short[]) in2 );
      else {
	 if ( ( (Object[]) in1 ).length != ( (Object[]) in2 ).length ) return false;
	 for ( int i = 0, s = ( (Object[]) in1 ).length; i < s; i ++ ) 
	    if ( ! Equals( ( (Object[]) in1 )[ i ], ( (Object[]) in2 )[ i ] ) ) return false; 
	 return true; } }

   public static int NEqualTo( Object in, short val ) { 
      int out = 0; 
      if ( in instanceof short[] ) { short[] inn = (short[]) in; 
      for ( int i = 0, s = inn.length; i < s; i ++ ) if ( inn[ i ] == val ) out ++; }
      else { for ( int i = 0, s = ( (Object[]) in ).length; i < s; i ++ ) 
	 out += NEqualTo( ( (Object[]) in )[ i ], val ); }
      return out; }

   public static int NNotEqualTo( Object in, short val ) { 
      return NElem( in ) - NEqualTo( in, val ); }

   public static boolean AnyEqualTo( Object in, short val ) { 
      if ( in instanceof short[] ) { short[] inn = (short[]) in; 
      for ( int i = 0, s = inn.length; i < s; i ++ ) if ( inn[ i ] == val ) return true; }
      else { for ( int i = 0, s = ( (Object[]) in ).length; i < s; i ++ ) 
	 if ( AnyEqualTo( ( (Object[]) in )[ i ], val ) ) return true; }
      return false; }   

   public static boolean AnyNotEqualTo( Object in, short val ) { 
      if ( in instanceof short[] ) { short[] inn = (short[]) in; 
      for ( int i = 0, s = inn.length; i < s; i ++ ) if ( inn[ i ] != val ) return true; }
      else { for ( int i = 0, s = ( (Object[]) in ).length; i < s; i ++ ) 
	 if ( AnyNotEqualTo( ( (Object[]) in )[ i ], val ) ) return true; }
      return false; }   

   public static final void Shift( short[] in, int amt, short spare ) {
      if ( amt > 0 ) { System.arraycopy( in, 0, in, amt, in.length - amt ); Arrays.fill( in, 0, amt, spare ); }
      else if ( amt < 0 ) { System.arraycopy( in, -amt, in, 0, in.length + amt ); Arrays.fill( in, in.length+amt, in.length, spare ); } }
   public static final void Shift( short[][] in, int amt, short spare ) {
      for ( int i = 0, s = in.length; i < s; i ++ ) Shift( in[ i ], amt, spare ); }

   public static final void Mask( short[] in, boolean[] flags, boolean val ) {
      for ( int i = 0, s = in.length; i < s; i ++ ) if ( flags[ i ] == val ) in[ i ] = 0; }
   public static final int IndexOf( short[] in, short val ) {
      for ( int i = 0, s = in.length; i < s; i ++ ) if ( in[ i ] == val ) return i; return -1; }
   public static final boolean Contains( short[] in, short val ) {
      return ( IndexOf( in, val ) != -1 ); }
   public static final int LastIndexOf( short[] in, short val ) {
      for ( int i = in.length - 1; i >= 0; i -- ) if ( in[ i ] == val ) return i; return -1; }
   public static final short[] Index( short in[], int index[] ) {
      short out[] = New( in.length ); 
      for ( int i = 0, s = in.length; i < s; i ++ ) out[ i ] = in[ index[ i ] ]; return out; }

   public static final short[] Reverse( short[] in ) {
      for ( int i = 0, s = in.length; i < s/2; i ++ ) {
	 short temp = in[ i ]; in[ i ] = in[ s-i-1 ]; in[ s-i-1 ] = temp; }
      return in; }

   // NO BOOL
   public static final void Sort( short[] in ) { Arrays.sort( in ); }
   public static final int[] Indexx( short[] in, int[] out ) {
      double darr[] = Convert( in );
      try { return djr.nr.Indexx.indexx( darr, out ); } catch( djr.nr.NRException e ) { }; 
      return null; }
   public static final int[] IRank( short[] in, int[] out ) {
      double darr[] = Convert( in );
      try { return djr.nr.Indexx.irank( darr, out ); } catch( djr.nr.NRException e ) { }; 
      return null; }

   public static final int Search( short[] in, short val, boolean sort ) {
      if ( sort ) Sort( in ); return Arrays.binarySearch( in, val ); }
   public static final int Search( short[] in, short val ) {
      return Search( in, val, true ); }

   public static final short[] Uniq( short[] in ) {
      ShortVector v = new ShortVector();
      for ( int i = 0, s = in.length; i < s; i ++ ) 
	 if ( ! v.contains( in[ i ] ) ) v.addElement( (short) i );
      return v.data();
   }
   // END NO BOOL

   public static final int NElem( Object in ) { 
      if ( in instanceof short[] ) return ( (short[]) in ).length; 
      else { int sum = 0; for ( int i = 0, s = ( (Object[]) in ).length; i < s; i ++ ) 
	 sum += NElem( ( (Object[]) in )[ i ] ); return sum; } }
      
   // NO BOOL
   public static final double Mean( Object in ) {
      return ( (double) Sum( in ) / (double) NElem( in ) ); }

   public static final double Mean( Object in, int i1, int i2 ) {
      return ( (double) Sum( in, i1, i2 ) / (double) ( i2 - i1 + 1 ) ); }

   public static final double Median( short[] in ) {
      final int nx = in.length; short temp[] = New( nx );
      Copy( temp, in ); Sort( temp ); double out = (double) temp[ (nx+1)/2 ];
      if ( nx % 2 != 1 ) out = ( (double) temp[ nx/2 ] + (double) temp[ nx/2+1 ] ) / 2.0;
      return out; }

   public static final double Variance( short[] in, double mean ) {
      double sum = 0.0; int nelem = 0;
      for ( int i = 0, s = in.length; i < s; i ++ ) 
	 sum += ( in[ i ] - mean ) * ( in[ i ] - mean );
      return sum / ( (double) in.length - 1.0 ); }

   public static final double Stddev( short[] in, double mean ) {
      return Math.sqrt( Variance( in, mean ) ); }

   public static final double Stddev( short[] in ) {
      return Stddev( in, Mean( in ) ); }

   public static final void Stats( short[] in, double[] out ) {
      out[ 0 ] = Mean( in ); out[ 1 ] = Median( in ); out[ 2 ] = Stddev( in ); }
   public static final double[] Stats( short[] in ) {
      double out[] = new double[ 3 ]; Stats( in, out ); return out; }

   public static void PrintStats( short[] data ) {
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

   public static final double[][] Histogram( short[] in, int nbins ) {
      double min = (double) Min( in ), max = (double) Max( in );
      double binsize = ( max - min ) / ( (double) nbins ), bs2 = binsize / 2.0;
      double[][] out = DoubleUtils.New( 2, nbins + 1 ); int ind = 0;
      for ( double i = min; i <= max; i += binsize ) {
	 out[ 0 ][ ind ] = i;
	 for ( int j = 0, size = in.length; j < size; j ++ )
	    if ( in[ j ] != 0 && in[ j ] >= i - bs2 && in[ j ] < i + bs2 ) out[ 1 ][ ind ] ++;
	 if ( ind ++ >= nbins ) break; }
      return out; }

   public static final double[] Histogram( short[] in, short[] bins ) {
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

   public static final double Mode( short[] in ) {
      int nbins = Math.min( (int) ( in.length / 20.0 ), 20 );
      double[][] hist = Histogram( in, nbins );
      return hist[ 0 ][ DoubleUtils.WhereMax( hist[ 1 ] ) ]; }

   public static final short Max( short x1, short x2 ) {
      return ( x1 > x2 ? x1 : x2 ); }
   public static final short Max( short x1, short x2, short x3 ) {
      return Max( x1, Max( x2, x3 ) ); }
   public static final short Max( short x1, short x2, short x3, short x4 ) {
      return Max( Max( x1, x2 ), Max( x3, x4 ) ); }
   public static final short Min( short x1, short x2 ) {
      return ( x1 < x2 ? x1 : x2 ); }
   public static final boolean Approx( short x1, short x2 ) {
      return ( Math.abs( x1 - x2 ) < 1e-8 ); }

   public static final int Sample( short[] x ) {
      return Sample( x, x.length, 500 ); }
   public static final int Sample( short[] x, int len ) {
      return Sample( x, len, 500 ); }

   // Assumes x[] is normalized so its max. value is 1.
   public static final int Sample( short[] x, int len, int maxCount ) {
      int r1 = -1, count = 0; double r2, xx; len --;
      while( count ++ < maxCount ) {
	 r1 = IntUtils.RandChoose( 0, len );
	 xx = x[ r1 ];
	 if ( xx <= 0.0 ) { r1 = -1; continue; }
	 r2 = DoubleUtils.Random(); if ( r2 <= xx ) break;
	 r1 = -1;
      }
      return r1; }

   public static final double[] Convolve( short x[], short kern[], double out[] ) {
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

   public static final double logGamma( short xx ) {
      return cern.jet.stat.Gamma.logGamma( (double) xx );
   }
   // END NO BOOL

   public static final short[] FTokenize( String str, final String tok ) { 
      return FTokenize( MyUtils.Tokenize( str, tok ) ); }

   public static final short[] FTokenize( String str[] ) { 
      return FTokenize( str, 0, str.length ); }

   public static final short[] FTokenize( String str[], int start, int end ) { 
      short[] out = New( end - start );
      for ( int i = start; i < end; i ++ ) {
   // NO BOOL
	 try { out[ i - start ] = (short) Short.parseShort( str[ i ] ) ; } 
   // END NO BOOL
   /* BOOL ONLY
	 try { out[ i - start ] = Short.valueOf( str[ i ] ).ShortValue() ; } 
   // END BOOL ONLY */
	 catch( Exception e ) { out[ i ] = 0; } } return out; }

   public static final short[] FromVector( Vector v ) {
      if ( v == null || v.size() <= 0 ) return null;
      short out[] = New( v.size() );
      try {
	 for ( int i = 0, s = v.size(); i < s; i ++ )
	    out[ i ] = ( (Short) v.elementAt( i ) ).shortValue();
      } catch( Exception e ) { };
      return out;
   }

   public static final Vector ToVector( short[] in ) {
      if ( in == null ) return null;
      Vector v = new Vector();
      for ( int i = 0, s = in.length; i < s; i ++ ) v.addElement( new Short ( in[ i ] ) );
      return v;
   }
}
