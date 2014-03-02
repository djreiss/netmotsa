package djr.util.bio;

import djr.util.array.ShortUtils;
import djr.util.ANSI;
import djr.util.HTML;

/**
 * <code>MaskedSequence</code> class.
 *
 * @author <a href="mailto:dreiss@systemsbiology.org">David Reiss</a>
 * @version 1.9978 (Fri Nov 07 05:56:26 PST 2003)
 */
public class MaskedSequence extends Sequence {
   short mask[];
   
   public MaskedSequence( String s, String h ) {
      super( s, h );
      mask = ShortUtils.New( GetLength() );
   }

   public MaskedSequence( String s ) {
      super( s );
      mask = ShortUtils.New( GetLength() );
   }

   public MaskedSequence( short[] res ) {
      this( res, null ); 
   }

   public MaskedSequence( short[] res, String head ) {
      super( res, head );
      mask = ShortUtils.New( GetLength() );
   }

   public MaskedSequence( Sequence in ) {
      super( in );
      if ( in instanceof MaskedSequence ) 
	 mask = ShortUtils.New( ( (MaskedSequence) in ).mask );
      else mask = ShortUtils.New( GetLength() );
   }

   public short[] GetMask() { return mask; }
   public short GetMask( int ind ) { return mask[ ind ]; }
   public void SetMask( int ind, int val ) { mask[ ind ] = (short) val; }
   public void SetMask( int ind ) { SetMask( ind, 1 ); }
   public void ClearMask( int ind ) { SetMask( ind, 0 ); }
   public void MaskRange( int start, int end, int val ) { 
      ShortUtils.Set( mask, start, end + 1, (short) val ); }
   public void MaskRange( int start, int end ) { MaskRange( start, end, 1 ); }
   public void ClearRange( int start, int end ) { MaskRange( start, end, 0 ); }
   public boolean IsMasked( int ind ) { return mask[ ind ] != 0; }
   public void ClearMask() { ShortUtils.Zero( mask ); }
   public boolean HasMask() { return ShortUtils.Sum( mask ) > 0; }
   public void SetMask( int inds[], int val ) { 
      ShortUtils.Set( mask, inds, (short) val ); }
   public void SetMask( int inds[] ) { SetMask( inds, 1 ); }
   public void ClearMask( int inds[], int val ) { SetMask( inds, 0 ); }
   public void ApplyMask() { ApplyMask( GetResidues() ); } // THIS IS DESTRUCTIVE!!!

   public void ApplyMask( short r[] ) { 
      for ( int i = 0, sz = GetLength(); i < sz; i ++ ) 
	 if ( GetMask( i ) != 0 ) r[ i ] = -1;
   }

   public boolean IsMaskedInRange( int start, int end ) {
      if ( start < 0 ) start = 0;
      if ( end >= GetLength() ) end = GetLength() - 1;
      for ( int i = start; i <= end; i ++ ) if ( mask[ i ] > 0 ) return true; 
      return false; }

   public String ToANSIString() {
      char c[] = GetSequence().toCharArray();
      short J = GetAlphabetSize();
      StringBuffer out = new StringBuffer();
      for ( int i = 0, s = c.length; i < s; i ++ ) {
	 if ( ! IsMasked( i ) ) out.append( setANSIColorFG( c[ i ], J ) );
	 else out.append( setANSIColorBG( c[ i ], J ) );
      }
      return out.toString();
   }

   public String ToHTMLString( boolean tableWrapper ) {
      char c[] = GetSequence().toCharArray();
      short J = GetAlphabetSize();
      StringBuffer out = new StringBuffer();
      if ( tableWrapper ) out.append( HTML.tableStart() );
      for ( int i = 0, s = c.length; i < s; i ++ ) {
	 boolean masked = IsMasked( i );
	 if ( masked ) {
	    out.append( setHTMLColorBG( c[ i ], J ) );
	 } else {
	    out.append( HTML.cellStartA( "center" ) );
	    out.append( setHTMLColorFG( c[ i ], J ) );
	    out.append( HTML.cellEnd() );
	 }
      }
      if ( tableWrapper ) out.append( HTML.reset() );
      return out.toString();
   }

   public String toString() {
      if ( "ANSI".equals( PRINT_FORMAT ) ) return ToANSIString();
      else if ( "HTML".equals( PRINT_FORMAT ) ) return ToHTMLString( true );
      else if ( "PLAIN".equals( PRINT_FORMAT ) ) {
	 char c[] = GetSequence().toCharArray();
	 short r[] = GetResidues(), J = GetAlphabetSize();
	 StringBuffer out = new StringBuffer();
	 for ( int i = 0, s = c.length; i < s; i ++ ) {
	    if ( IsMasked( i ) ) out.append( Character.toLowerCase( c[ i ] ) );
	    else out.append( c[ i ] );
	 }
	 return out.toString();
      }
      return super.toString();
      /*String out = super.toString();
      if ( HasMask() ) {
	 out += "\nMASKED: ";
	 for ( int i = 0; i < length; i ++ ) {
	    if ( IsMasked( i ) ) {
	       out += i;
	       int start = i;
	       while( IsMasked( ++ i ) );
	       if ( i > start + 1 ) out += " => " + ( i - 1 );
	       out += "; ";
	    }
	 }
	 return out.substring( 0, out.length() - 2 );
      }
      return out;*/
   }

   /*public static void main( String args[] ) {
      PRINT_FORMAT = "ANSI";
      Sequence S[] = RunTest( args );
      MaskedSequence ss = new MaskedSequence( S[ 0 ] );
      ss.MaskRange( 3, 7 );
      ss.SetMask( 11 );
      System.out.println( ss );
      ss.ClearMask( 5 );
      System.out.println( ss );
      ss.ClearMask();
      System.out.println( ss );
      ss.SetMask( new int[] { 3, 7, 9, 12, 13, 18 }, 1 );
      System.out.println( ss );
      }*/
}
