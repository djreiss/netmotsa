package djr.util.array;

/**
 * Template class <code>ObjVector</code>
 * See djr.util.Utils.java for a description of what goes on here.
 *
 * @author <a href="mailto:dreiss@systemsbiology.org">David Reiss</a>
 * @version 1.9978 (Fri Nov 07 05:56:26 PST 2003)
 */
public class ObjVector extends MyVector {
   protected Object[] data;

   public ObjVector() {
      this( 10 );
   }

   public ObjVector( int initialSize ) {
      this( initialSize, -1 );
   }

   public ObjVector( int initialSize, int growth ) {
      data = new Object[ initialSize ];
      this.inuse = 0;
      this.growth = growth;
      this.initialSize = initialSize;
   }

   public ObjVector( Object arr[] ) {
      this( arr.length );
      System.arraycopy( arr, 0, data, 0, arr.length );
      this.inuse = arr.length;
      this.initialSize = arr.length;
   }

   public void add( Object value ) {
      if ( size() >= data.length ) {
	 int leng = data.length;
	 Object temp[] = new Object[ growth <= 0 && leng > 0 ? leng * 2 : leng + growth ];
	 System.arraycopy( data, 0, temp, 0, size() );
	 data = temp;
      }
      data[ inuse ++ ] = value;
   }

   public void addElement( Object value ) {
      add( value );
   }

   public Object get( int which ) {
      if ( which < 0 || which > data.length ) return null;
      return data[ which ];
   }

   public Object elementAt( int which ) /*throws ArrayIndexOutOfBoundsException*/ {
      //if ( which < 0 || which >= size() ) 
      //	 throw new ArrayIndexOutOfBoundsException( which );
      return get( which );
   }

   public Object set( int which, Object value ) {
      if ( which < 0 || which > data.length ) return null;
      data[ which ] = value;
      return value;
   }

   public Object setElementAt( int which, Object value ) /*throws ArrayIndexOutOfBoundsException*/ {
      //if ( which < 0 || which >= size() ) 
      // throw new ArrayIndexOutOfBoundsException( which );
      return set( which, value );
   }

   public Object remove( int which ) {
      if ( which < 0 || which > data.length ) return null;
      int j = size() - which - 1;
      Object out = data[ which ];
      if ( j > 0 ) System.arraycopy( data, which + 1, data, which, j);
      data[ -- inuse ] = null;
      return out;
   }

   public Object removeFast( int which ) {
      if ( which < 0 || which > data.length ) return null;
      Object out = data[ which ];
      data[ which ] = data[ inuse - 1 ];
      data[ -- inuse ] = null;
      return out;
   }

   public Object removeElementAt( int which ) /*throws ArrayIndexOutOfBoundsException*/ {
      //if ( which < 0 || which >= size() ) 
      //   throw new ArrayIndexOutOfBoundsException( which );
      return remove( which );
   }

   public boolean removeAll() {
      if ( size() <= 0 ) return false;
      inuse = 0;
      //java.util.Arrays.fill( data, (Object) null );
      data = new Object[ initialSize ];
      return true;
   }

   public boolean clear() { return removeAll(); }
   public boolean erase() { return removeAll(); }

   public int indexOf( Object val ) {
      for ( int i = 0, s = size(); i < s; i ++ ) if ( val.equals( data[ i ] ) ) return i;
      //int out = java.util.Arrays.binarySearch( (Object[]) data, val );
      //if ( out >= 0 && out < size() ) return out;
      return -1;
   }

   public boolean contains( Object val ) {
      return indexOf( val ) >= 0;
   }

   public boolean delete( Object val ) {
      return removeElement( val );
   }

   public boolean deleteFast( Object val ) {
      int ind = indexOf( val );
      if ( ind >= 0 ) { removeFast( ind ); return true; }
      return false;
   }

   public boolean removeElement( Object val ) {
      int ind = indexOf( val );
      if ( ind >= 0 ) { remove( ind ); return true; }
      return false;
   }

   public int capacity() { return data.length; }

   /*public int search( Object o ) {
      int i = ObjUtils.LastIndexOf( data, o );
      if ( i >= 0 ) return size() - i;
      return -1;
      }*/

   public Object[] toArray() { return data(); }

   public Object[] rawData() { return data; }

   public Object[] data() {
      if ( capacity() == size() ) return data; 
      Object out[] = new Object[ size() ];
      System.arraycopy( data, 0, out, 0, size() );
      return out;
   }

   public Object[] data( int start, int end ) {
      int len = end - start + 1;
      if ( len > size() ) return data();
      if ( end >= size() ) { end = size() - 1; len = end - start + 1; }
      Object out[] = new Object[ len ];
      System.arraycopy( data, start, out, 0, len );
      return out;
   }

   public Object[] data( Object out[] ) { 
      if ( out == null || out.length < size() ) return data();
      System.arraycopy( data, 0, out, 0, size() );
      return out;
   }

   public void removeEveryOther( int steps ) {
      for ( int i = size()-1; i > 0; i -= steps ) removeElementAt( i );
   }

   public String toString( int ind ) {
      return get( ind ).toString();
   }

   public java.util.Enumeration elements() {
      return new java.util.Enumeration() { int cur = 0;
	    public boolean hasMoreElements() { return cur < size(); }
	    public Object nextElement() { return get( cur ++ ); }
	    public Object nextElement2() { return get( cur ++ ); }
	 };
   }
}
