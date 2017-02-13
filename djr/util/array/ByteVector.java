package djr.util.array;

/**
 * Template class <code>ByteVector</code>
 * See djr.util.Utils.java for a description of what goes on here.
 *
 * @author <a href="mailto:reiss@uw.edu">David Reiss</a>
 * @version 1.9978 (Fri Nov 07 05:56:26 PST 2003)
 */
public class ByteVector extends MyVector {
   protected byte[] data;

   public ByteVector() {
      this( 10 );
   }

   public ByteVector( int initialSize ) {
      this( initialSize, -1 );
   }

   public ByteVector( int initialSize, int growth ) {
      data = new byte[ initialSize ];
      this.inuse = 0;
      this.growth = growth;
      this.initialSize = initialSize;
   }

   public ByteVector( byte arr[] ) {
      this( arr.length );
      System.arraycopy( arr, 0, data, 0, arr.length );
      this.inuse = arr.length;
      this.initialSize = arr.length;
   }

   public void add( byte value ) {
      if ( size() >= data.length ) {
	 int leng = data.length;
	 byte temp[] = new byte[ growth <= 0 && leng > 0 ? leng * 2 : leng + growth ];
	 System.arraycopy( data, 0, temp, 0, size() );
	 data = temp;
      }
      data[ inuse ++ ] = value;
   }

   public void addElement( byte value ) {
      add( value );
   }

   public byte get( int which ) {
      if ( which < 0 || which > data.length ) return 0;
      return data[ which ];
   }

   public byte elementAt( int which ) /*throws ArrayIndexOutOfBoundsException*/ {
      //if ( which < 0 || which >= size() ) 
      //	 throw new ArrayIndexOutOfBoundsException( which );
      return get( which );
   }

   public byte set( int which, byte value ) {
      if ( which < 0 || which > data.length ) return 0;
      data[ which ] = value;
      return value;
   }

   public byte setElementAt( int which, byte value ) /*throws ArrayIndexOutOfBoundsException*/ {
      //if ( which < 0 || which >= size() ) 
      // throw new ArrayIndexOutOfBoundsException( which );
      return set( which, value );
   }

   public byte remove( int which ) {
      if ( which < 0 || which > data.length ) return 0;
      int j = size() - which - 1;
      byte out = data[ which ];
      if ( j > 0 ) System.arraycopy( data, which + 1, data, which, j);
      data[ -- inuse ] = 0;
      return out;
   }

   public byte removeFast( int which ) {
      if ( which < 0 || which > data.length ) return 0;
      byte out = data[ which ];
      data[ which ] = data[ inuse - 1 ];
      data[ -- inuse ] = 0;
      return out;
   }

   public byte removeElementAt( int which ) /*throws ArrayIndexOutOfBoundsException*/ {
      //if ( which < 0 || which >= size() ) 
      //   throw new ArrayIndexOutOfBoundsException( which );
      return remove( which );
   }

   public boolean removeAll() {
      if ( size() <= 0 ) return false;
      inuse = 0;
      //java.util.Arrays.fill( data, (byte) 0 );
      data = new byte[ initialSize ];
      return true;
   }

   public boolean clear() { return removeAll(); }
   public boolean erase() { return removeAll(); }

   public int indexOf( byte val ) {
      for ( int i = 0, s = size(); i < s; i ++ ) if ( val == data[ i ] ) return i;
      //int out = java.util.Arrays.binarySearch( (byte[]) data, val );
      //if ( out >= 0 && out < size() ) return out;
      return -1;
   }

   public boolean contains( byte val ) {
      return indexOf( val ) >= 0;
   }

   public boolean delete( byte val ) {
      return removeElement( val );
   }

   public boolean deleteFast( byte val ) {
      int ind = indexOf( val );
      if ( ind >= 0 ) { removeFast( ind ); return true; }
      return false;
   }

   public boolean removeElement( byte val ) {
      int ind = indexOf( val );
      if ( ind >= 0 ) { remove( ind ); return true; }
      return false;
   }

   public int capacity() { return data.length; }

   /*public int search( byte o ) {
      int i = ByteUtils.LastIndexOf( data, o );
      if ( i >= 0 ) return size() - i;
      return -1;
      }*/

   public byte[] toArray() { return data(); }

   public byte[] rawData() { return data; }

   public byte[] data() {
      if ( capacity() == size() ) return data; 
      byte out[] = new byte[ size() ];
      System.arraycopy( data, 0, out, 0, size() );
      return out;
   }

   public byte[] data( int start, int end ) {
      int len = end - start + 1;
      if ( len > size() ) return data();
      if ( end >= size() ) { end = size() - 1; len = end - start + 1; }
      byte out[] = new byte[ len ];
      System.arraycopy( data, start, out, 0, len );
      return out;
   }

   public byte[] data( byte out[] ) { 
      if ( out == null || out.length < size() ) return data();
      System.arraycopy( data, 0, out, 0, size() );
      return out;
   }

   public void removeEveryOther( int steps ) {
      for ( int i = size()-1; i > 0; i -= steps ) removeElementAt( i );
   }

   public String toString( int ind ) {
      return Byte.toString( get( ind ) );
   }

   public java.util.Enumeration elements() {
      return new java.util.Enumeration() { int cur = 0;
	    public boolean hasMoreElements() { return cur < size(); }
	    public Object nextElement() { return new Byte( get( cur ++ ) ); }
	    public byte nextElement2() { return get( cur ++ ); }
	 };
   }
}
