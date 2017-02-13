package djr.util.array;

/**
 * Abstract class <code>MyVector</code>
 * Abstract superclass for the TYPEVector classes created from TTYPEVector.jtempl
 *
 * @author <a href="mailto:reiss@uw.edu">David Reiss</a>
 * @version 1.9978 (Fri Nov 07 05:56:26 PST 2003)
 */
public abstract class MyVector implements java.io.Serializable {
   protected int inuse, growth, initialSize;
   
   public abstract boolean removeAll();

   public boolean empty() {
      return removeAll();
   }

   public int size() { return inuse; }

   public int length() { return size(); }

   public boolean isEmpty() {
      return size() == 0;
   }

   public abstract String toString( int ind );

   public String toString() {
      String out = ""; //this.getClass().getName();
      out = out + "["; //.substring( out.lastIndexOf( '.' ) + 1 ) + "{";
      if ( size() <= 0 ) out += "]"; //"}";
      else {
	 for ( int i = 0, s = size(); i < s - 1; i ++ ) out += toString( i ) + "; ";
	 out += toString( size() - 1 ) + "]"; //"}";
      }
      return out;
   }
}
