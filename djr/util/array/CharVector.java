package djr.util.array;

/**
 * Template class <code>CharVector</code>
 * See djr.util.Utils.java for a description of what goes on here.
 *
 * @author <a href="mailto:reiss@uw.edu">David Reiss</a>
 * @version 1.9978 (Fri Nov 07 05:56:26 PST 2003)
 */
public class CharVector extends MyVector {
   protected char[] data;

   public CharVector() {
      this( 10 );
   }

   public CharVector( int initialSize ) {
      this( initialSize, -1 );
   }

   public CharVector( int initialSize, int growth ) {
      data = new char[ initialSize ];
      this.inuse = 0;
      this.growth = growth;
      this.initialSize = initialSize;
   }

   public CharVector( char arr[] ) {
      this( arr.length );
      System.arraycopy( arr, 0, data, 0, arr.length );
      this.inuse = arr.length;
      this.initialSize = arr.length;
   }

   public void add( char value ) {
      if ( size() >= data.length ) {
	 int leng = data.length;
	 char temp[] = new char[ growth <= 0 && leng > 0 ? leng * 2 : leng + growth ];
	 System.arraycopy( data, 0, temp, 0, size() );
	 data = temp;
      }
      data[ inuse ++ ] = value;
   }

   public void addElement( char value ) {
      add( value );
   }

   public char get( int which ) {
      if ( which < 0 || which > data.length ) return 0;
      return data[ which ];
   }

   public char elementAt( int which ) /*throws ArrayIndexOutOfBoundsException*/ {
      //if ( which < 0 || which >= size() ) 
      //	 throw new ArrayIndexOutOfBoundsException( which );
      return get( which );
   }

   public char set( int which, char value ) {
      if ( which < 0 || which > data.length ) return 0;
      data[ which ] = value;
      return value;
   }

   public char setElementAt( int which, char value ) /*throws ArrayIndexOutOfBoundsException*/ {
      //if ( which < 0 || which >= size() ) 
      // throw new ArrayIndexOutOfBoundsException( which );
      return set( which, value );
   }

   public char remove( int which ) {
      if ( which < 0 || which > data.length ) return 0;
      int j = size() - which - 1;
      char out = data[ which ];
      if ( j > 0 ) System.arraycopy( data, which + 1, data, which, j);
      data[ -- inuse ] = 0;
      return out;
   }

   public char removeFast( int which ) {
      if ( which < 0 || which > data.length ) return 0;
      char out = data[ which ];
      data[ which ] = data[ inuse - 1 ];
      data[ -- inuse ] = 0;
      return out;
   }

   public char removeElementAt( int which ) /*throws ArrayIndexOutOfBoundsException*/ {
      //if ( which < 0 || which >= size() ) 
      //   throw new ArrayIndexOutOfBoundsException( which );
      return remove( which );
   }

   public boolean removeAll() {
      if ( size() <= 0 ) return false;
      inuse = 0;
      //java.util.Arrays.fill( data, (char) 0 );
      data = new char[ initialSize ];
      return true;
   }

   public boolean clear() { return removeAll(); }
   public boolean erase() { return removeAll(); }

   public int indexOf( char val ) {
      for ( int i = 0, s = size(); i < s; i ++ ) if ( val == data[ i ] ) return i;
      //int out = java.util.Arrays.binarySearch( (char[]) data, val );
      //if ( out >= 0 && out < size() ) return out;
      return -1;
   }

   public boolean contains( char val ) {
      return indexOf( val ) >= 0;
   }

   public boolean delete( char val ) {
      return removeElement( val );
   }

   public boolean deleteFast( char val ) {
      int ind = indexOf( val );
      if ( ind >= 0 ) { removeFast( ind ); return true; }
      return false;
   }

   public boolean removeElement( char val ) {
      int ind = indexOf( val );
      if ( ind >= 0 ) { remove( ind ); return true; }
      return false;
   }

   public int capacity() { return data.length; }

   /*public int search( char o ) {
      int i = CharUtils.LastIndexOf( data, o );
      if ( i >= 0 ) return size() - i;
      return -1;
      }*/

   public char[] toArray() { return data(); }

   public char[] rawData() { return data; }

   public char[] data() {
      if ( capacity() == size() ) return data; 
      char out[] = new char[ size() ];
      System.arraycopy( data, 0, out, 0, size() );
      return out;
   }

   public char[] data( int start, int end ) {
      int len = end - start + 1;
      if ( len > size() ) return data();
      if ( end >= size() ) { end = size() - 1; len = end - start + 1; }
      char out[] = new char[ len ];
      System.arraycopy( data, start, out, 0, len );
      return out;
   }

   public char[] data( char out[] ) { 
      if ( out == null || out.length < size() ) return data();
      System.arraycopy( data, 0, out, 0, size() );
      return out;
   }

   public void removeEveryOther( int steps ) {
      for ( int i = size()-1; i > 0; i -= steps ) removeElementAt( i );
   }

   public String toString( int ind ) {
      return Character.toString( get( ind ) );
   }

   public java.util.Enumeration elements() {
      return new java.util.Enumeration() { int cur = 0;
	    public boolean hasMoreElements() { return cur < size(); }
	    public Object nextElement() { return new Character( get( cur ++ ) ); }
	    public char nextElement2() { return get( cur ++ ); }
	 };
   }
}
