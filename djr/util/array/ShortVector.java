package djr.util.array;

/**
 * Template class <code>ShortVector</code>
 * See djr.util.Utils.java for a description of what goes on here.
 *
 * @author <a href="mailto:reiss@uw.edu">David Reiss</a>
 * @version 1.9978 (Fri Nov 07 05:56:26 PST 2003)
 */
public class ShortVector extends MyVector {
   protected short[] data;

   public ShortVector() {
      this( 10 );
   }

   public ShortVector( int initialSize ) {
      this( initialSize, -1 );
   }

   public ShortVector( int initialSize, int growth ) {
      data = new short[ initialSize ];
      this.inuse = 0;
      this.growth = growth;
      this.initialSize = initialSize;
   }

   public ShortVector( short arr[] ) {
      this( arr.length );
      System.arraycopy( arr, 0, data, 0, arr.length );
      this.inuse = arr.length;
      this.initialSize = arr.length;
   }

   public void add( short value ) {
      if ( size() >= data.length ) {
	 int leng = data.length;
	 short temp[] = new short[ growth <= 0 && leng > 0 ? leng * 2 : leng + growth ];
	 System.arraycopy( data, 0, temp, 0, size() );
	 data = temp;
      }
      data[ inuse ++ ] = value;
   }

   public void addElement( short value ) {
      add( value );
   }

   public short get( int which ) {
      if ( which < 0 || which > data.length ) return 0;
      return data[ which ];
   }

   public short elementAt( int which ) /*throws ArrayIndexOutOfBoundsException*/ {
      //if ( which < 0 || which >= size() ) 
      //	 throw new ArrayIndexOutOfBoundsException( which );
      return get( which );
   }

   public short set( int which, short value ) {
      if ( which < 0 || which > data.length ) return 0;
      data[ which ] = value;
      return value;
   }

   public short setElementAt( int which, short value ) /*throws ArrayIndexOutOfBoundsException*/ {
      //if ( which < 0 || which >= size() ) 
      // throw new ArrayIndexOutOfBoundsException( which );
      return set( which, value );
   }

   public short remove( int which ) {
      if ( which < 0 || which > data.length ) return 0;
      int j = size() - which - 1;
      short out = data[ which ];
      if ( j > 0 ) System.arraycopy( data, which + 1, data, which, j);
      data[ -- inuse ] = 0;
      return out;
   }

   public short removeFast( int which ) {
      if ( which < 0 || which > data.length ) return 0;
      short out = data[ which ];
      data[ which ] = data[ inuse - 1 ];
      data[ -- inuse ] = 0;
      return out;
   }

   public short removeElementAt( int which ) /*throws ArrayIndexOutOfBoundsException*/ {
      //if ( which < 0 || which >= size() ) 
      //   throw new ArrayIndexOutOfBoundsException( which );
      return remove( which );
   }

   public boolean removeAll() {
      if ( size() <= 0 ) return false;
      inuse = 0;
      //java.util.Arrays.fill( data, (short) 0 );
      data = new short[ initialSize ];
      return true;
   }

   public boolean clear() { return removeAll(); }
   public boolean erase() { return removeAll(); }

   public int indexOf( short val ) {
      for ( int i = 0, s = size(); i < s; i ++ ) if ( val == data[ i ] ) return i;
      //int out = java.util.Arrays.binarySearch( (short[]) data, val );
      //if ( out >= 0 && out < size() ) return out;
      return -1;
   }

   public boolean contains( short val ) {
      return indexOf( val ) >= 0;
   }

   public boolean delete( short val ) {
      return removeElement( val );
   }

   public boolean deleteFast( short val ) {
      int ind = indexOf( val );
      if ( ind >= 0 ) { removeFast( ind ); return true; }
      return false;
   }

   public boolean removeElement( short val ) {
      int ind = indexOf( val );
      if ( ind >= 0 ) { remove( ind ); return true; }
      return false;
   }

   public int capacity() { return data.length; }

   /*public int search( short o ) {
      int i = ShortUtils.LastIndexOf( data, o );
      if ( i >= 0 ) return size() - i;
      return -1;
      }*/

   public short[] toArray() { return data(); }

   public short[] rawData() { return data; }

   public short[] data() {
      if ( capacity() == size() ) return data; 
      short out[] = new short[ size() ];
      System.arraycopy( data, 0, out, 0, size() );
      return out;
   }

   public short[] data( int start, int end ) {
      int len = end - start + 1;
      if ( len > size() ) return data();
      if ( end >= size() ) { end = size() - 1; len = end - start + 1; }
      short out[] = new short[ len ];
      System.arraycopy( data, start, out, 0, len );
      return out;
   }

   public short[] data( short out[] ) { 
      if ( out == null || out.length < size() ) return data();
      System.arraycopy( data, 0, out, 0, size() );
      return out;
   }

   public void removeEveryOther( int steps ) {
      for ( int i = size()-1; i > 0; i -= steps ) removeElementAt( i );
   }

   public String toString( int ind ) {
      return Short.toString( get( ind ) );
   }

   public java.util.Enumeration elements() {
      return new java.util.Enumeration() { int cur = 0;
	    public boolean hasMoreElements() { return cur < size(); }
	    public Object nextElement() { return new Short( get( cur ++ ) ); }
	    public short nextElement2() { return get( cur ++ ); }
	 };
   }
}
