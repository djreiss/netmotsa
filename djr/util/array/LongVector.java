package djr.util.array;

/**
 * Template class <code>LongVector</code>
 * See djr.util.Utils.java for a description of what goes on here.
 *
 * @author <a href="mailto:dreiss@systemsbiology.org">David Reiss</a>
 * @version 1.9978 (Fri Nov 07 05:56:26 PST 2003)
 */
public class LongVector extends MyVector {
   protected long[] data;

   public LongVector() {
      this( 10 );
   }

   public LongVector( int initialSize ) {
      this( initialSize, -1 );
   }

   public LongVector( int initialSize, int growth ) {
      data = new long[ initialSize ];
      this.inuse = 0;
      this.growth = growth;
      this.initialSize = initialSize;
   }

   public LongVector( long arr[] ) {
      this( arr.length );
      System.arraycopy( arr, 0, data, 0, arr.length );
      this.inuse = arr.length;
      this.initialSize = arr.length;
   }

   public void add( long value ) {
      if ( size() >= data.length ) {
	 int leng = data.length;
	 long temp[] = new long[ growth <= 0 && leng > 0 ? leng * 2 : leng + growth ];
	 System.arraycopy( data, 0, temp, 0, size() );
	 data = temp;
      }
      data[ inuse ++ ] = value;
   }

   public void addElement( long value ) {
      add( value );
   }

   public long get( int which ) {
      if ( which < 0 || which > data.length ) return 0;
      return data[ which ];
   }

   public long elementAt( int which ) /*throws ArrayIndexOutOfBoundsException*/ {
      //if ( which < 0 || which >= size() ) 
      //	 throw new ArrayIndexOutOfBoundsException( which );
      return get( which );
   }

   public long set( int which, long value ) {
      if ( which < 0 || which > data.length ) return 0;
      data[ which ] = value;
      return value;
   }

   public long setElementAt( int which, long value ) /*throws ArrayIndexOutOfBoundsException*/ {
      //if ( which < 0 || which >= size() ) 
      // throw new ArrayIndexOutOfBoundsException( which );
      return set( which, value );
   }

   public long remove( int which ) {
      if ( which < 0 || which > data.length ) return 0;
      int j = size() - which - 1;
      long out = data[ which ];
      if ( j > 0 ) System.arraycopy( data, which + 1, data, which, j);
      data[ -- inuse ] = 0;
      return out;
   }

   public long removeFast( int which ) {
      if ( which < 0 || which > data.length ) return 0;
      long out = data[ which ];
      data[ which ] = data[ inuse - 1 ];
      data[ -- inuse ] = 0;
      return out;
   }

   public long removeElementAt( int which ) /*throws ArrayIndexOutOfBoundsException*/ {
      //if ( which < 0 || which >= size() ) 
      //   throw new ArrayIndexOutOfBoundsException( which );
      return remove( which );
   }

   public boolean removeAll() {
      if ( size() <= 0 ) return false;
      inuse = 0;
      //java.util.Arrays.fill( data, (long) 0 );
      data = new long[ initialSize ];
      return true;
   }

   public boolean clear() { return removeAll(); }
   public boolean erase() { return removeAll(); }

   public int indexOf( long val ) {
      for ( int i = 0, s = size(); i < s; i ++ ) if ( val == data[ i ] ) return i;
      //int out = java.util.Arrays.binarySearch( (long[]) data, val );
      //if ( out >= 0 && out < size() ) return out;
      return -1;
   }

   public boolean contains( long val ) {
      return indexOf( val ) >= 0;
   }

   public boolean delete( long val ) {
      return removeElement( val );
   }

   public boolean deleteFast( long val ) {
      int ind = indexOf( val );
      if ( ind >= 0 ) { removeFast( ind ); return true; }
      return false;
   }

   public boolean removeElement( long val ) {
      int ind = indexOf( val );
      if ( ind >= 0 ) { remove( ind ); return true; }
      return false;
   }

   public int capacity() { return data.length; }

   /*public int search( long o ) {
      int i = LongUtils.LastIndexOf( data, o );
      if ( i >= 0 ) return size() - i;
      return -1;
      }*/

   public long[] toArray() { return data(); }

   public long[] rawData() { return data; }

   public long[] data() {
      if ( capacity() == size() ) return data; 
      long out[] = new long[ size() ];
      System.arraycopy( data, 0, out, 0, size() );
      return out;
   }

   public long[] data( int start, int end ) {
      int len = end - start + 1;
      if ( len > size() ) return data();
      if ( end >= size() ) { end = size() - 1; len = end - start + 1; }
      long out[] = new long[ len ];
      System.arraycopy( data, start, out, 0, len );
      return out;
   }

   public long[] data( long out[] ) { 
      if ( out == null || out.length < size() ) return data();
      System.arraycopy( data, 0, out, 0, size() );
      return out;
   }

   public void removeEveryOther( int steps ) {
      for ( int i = size()-1; i > 0; i -= steps ) removeElementAt( i );
   }

   public String toString( int ind ) {
      return Long.toString( get( ind ) );
   }

   public java.util.Enumeration elements() {
      return new java.util.Enumeration() { int cur = 0;
	    public boolean hasMoreElements() { return cur < size(); }
	    public Object nextElement() { return new Long( get( cur ++ ) ); }
	    public long nextElement2() { return get( cur ++ ); }
	 };
   }
}
