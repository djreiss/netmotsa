package djr.util.array;

/**
 * Template class <code>BoolVector</code>
 * See djr.util.Utils.java for a description of what goes on here.
 *
 * @author <a href="mailto:dreiss@systemsbiology.org">David Reiss</a>
 * @version 1.9978 (Fri Nov 07 05:56:26 PST 2003)
 */
public class BoolVector extends MyVector {
   protected boolean[] data;

   public BoolVector() {
      this( 10 );
   }

   public BoolVector( int initialSize ) {
      this( initialSize, -1 );
   }

   public BoolVector( int initialSize, int growth ) {
      data = new boolean[ initialSize ];
      this.inuse = 0;
      this.growth = growth;
      this.initialSize = initialSize;
   }

   public BoolVector( boolean arr[] ) {
      this( arr.length );
      System.arraycopy( arr, 0, data, 0, arr.length );
      this.inuse = arr.length;
      this.initialSize = arr.length;
   }

   public void add( boolean value ) {
      if ( size() >= data.length ) {
	 int leng = data.length;
	 boolean temp[] = new boolean[ growth <= 0 && leng > 0 ? leng * 2 : leng + growth ];
	 System.arraycopy( data, 0, temp, 0, size() );
	 data = temp;
      }
      data[ inuse ++ ] = value;
   }

   public void addElement( boolean value ) {
      add( value );
   }

   public boolean get( int which ) {
      if ( which < 0 || which > data.length ) return false;
      return data[ which ];
   }

   public boolean elementAt( int which ) /*throws ArrayIndexOutOfBoundsException*/ {
      //if ( which < 0 || which >= size() ) 
      //	 throw new ArrayIndexOutOfBoundsException( which );
      return get( which );
   }

   public boolean set( int which, boolean value ) {
      if ( which < 0 || which > data.length ) return false;
      data[ which ] = value;
      return value;
   }

   public boolean setElementAt( int which, boolean value ) /*throws ArrayIndexOutOfBoundsException*/ {
      //if ( which < 0 || which >= size() ) 
      // throw new ArrayIndexOutOfBoundsException( which );
      return set( which, value );
   }

   public boolean remove( int which ) {
      if ( which < 0 || which > data.length ) return false;
      int j = size() - which - 1;
      boolean out = data[ which ];
      if ( j > 0 ) System.arraycopy( data, which + 1, data, which, j);
      data[ -- inuse ] = false;
      return out;
   }

   public boolean removeFast( int which ) {
      if ( which < 0 || which > data.length ) return false;
      boolean out = data[ which ];
      data[ which ] = data[ inuse - 1 ];
      data[ -- inuse ] = false;
      return out;
   }

   public boolean removeElementAt( int which ) /*throws ArrayIndexOutOfBoundsException*/ {
      //if ( which < 0 || which >= size() ) 
      //   throw new ArrayIndexOutOfBoundsException( which );
      return remove( which );
   }

   public boolean removeAll() {
      if ( size() <= 0 ) return false;
      inuse = 0;
      //java.util.Arrays.fill( data, (boolean) false );
      data = new boolean[ initialSize ];
      return true;
   }

   public boolean clear() { return removeAll(); }
   public boolean erase() { return removeAll(); }

   public int indexOf( boolean val ) {
      for ( int i = 0, s = size(); i < s; i ++ ) if ( val == data[ i ] ) return i;
      //int out = java.util.Arrays.binarySearch( (boolean[]) data, val );
      //if ( out >= 0 && out < size() ) return out;
      return -1;
   }

   public boolean contains( boolean val ) {
      return indexOf( val ) >= 0;
   }

   public boolean delete( boolean val ) {
      return removeElement( val );
   }

   public boolean deleteFast( boolean val ) {
      int ind = indexOf( val );
      if ( ind >= 0 ) { removeFast( ind ); return true; }
      return false;
   }

   public boolean removeElement( boolean val ) {
      int ind = indexOf( val );
      if ( ind >= 0 ) { remove( ind ); return true; }
      return false;
   }

   public int capacity() { return data.length; }

   /*public int search( boolean o ) {
      int i = BoolUtils.LastIndexOf( data, o );
      if ( i >= 0 ) return size() - i;
      return -1;
      }*/

   public boolean[] toArray() { return data(); }

   public boolean[] rawData() { return data; }

   public boolean[] data() {
      if ( capacity() == size() ) return data; 
      boolean out[] = new boolean[ size() ];
      System.arraycopy( data, 0, out, 0, size() );
      return out;
   }

   public boolean[] data( int start, int end ) {
      int len = end - start + 1;
      if ( len > size() ) return data();
      if ( end >= size() ) { end = size() - 1; len = end - start + 1; }
      boolean out[] = new boolean[ len ];
      System.arraycopy( data, start, out, 0, len );
      return out;
   }

   public boolean[] data( boolean out[] ) { 
      if ( out == null || out.length < size() ) return data();
      System.arraycopy( data, 0, out, 0, size() );
      return out;
   }

   public void removeEveryOther( int steps ) {
      for ( int i = size()-1; i > 0; i -= steps ) removeElementAt( i );
   }

   public String toString( int ind ) {
      return Boolean.toString( get( ind ) );
   }

   public java.util.Enumeration elements() {
      return new java.util.Enumeration() { int cur = 0;
	    public boolean hasMoreElements() { return cur < size(); }
	    public Object nextElement() { return new Boolean( get( cur ++ ) ); }
	    public boolean nextElement2() { return get( cur ++ ); }
	 };
   }
}
