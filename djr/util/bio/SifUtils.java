package djr.util.bio;

import java.io.*;
import java.util.*;

import djr.util.*;
import djr.util.array.*;

/**
 * <code>SifUtils</code> class.
 *
 * @author <a href="mailto:dreiss@systemsbiology.org">David Reiss</a>
 * @version 1.9978 (Fri Nov 07 05:56:26 PST 2003)
 */
public class SifUtils {
   public static Map ReadSif( String fname ) {
      return ReadSif( fname, true ); }

   // Read a sif interaction file in as a hashmap of vectors.
   public static Map ReadSif( String fname, boolean bothWays ) {
      Map interaxns = new HashMap();
      try {
	 ObjVector vlines = MyUtils.ReadFileLines( fname );
         Enumeration lines = vlines.elements();
	 int lnum = -1;
         while( lines.hasMoreElements() ) {
	    lnum ++;
            String line = (String) lines.nextElement();
            String toks[] = MyUtils.Tokenize( line, " " );
            String p1 = toks[ 0 ].toUpperCase();
	    //String p2 = "pp".equals( toks[ 1 ] ) ? toks[ 2 ].toUpperCase() : toks[ 1 ].toUpperCase();
	    String p2 = toks[ 2 ].toUpperCase();
	    p1 = replaceSynonym( p1 ).toUpperCase();
	    p2 = replaceSynonym( p2 ).toUpperCase();
	    AddInteraction( interaxns, p1, p2, bothWays );
	 }
      } catch( Exception e ) { e.printStackTrace(); }
      return interaxns;
   }

   public static ObjVector ReadProteinList( String listName ) {
      ObjVector out = new ObjVector();
      try {
	 Enumeration lines = MyUtils.ReadFileLines( listName ).elements();
         while( lines.hasMoreElements() ) {
            String line = (String) lines.nextElement();
            String toks[] = MyUtils.Tokenize( line, " " );
	    toks[ 0 ] = replaceSynonym( toks[ 0 ] );
	    if ( "true".equals( toks[ 1 ] ) ) out.addElement( toks[ 0 ] );
	 }
      } catch( Exception e ) { };
      return out;
   }

   public static void WriteSif( Map map, String outName, boolean bothWays ) {
      try {
	 PrintStream out = new PrintStream( MyUtils.OpenOutputFile( outName ) );
	 Map temp = (Map) MyUtils.DeepCopy( map );
	 Iterator it = temp.keySet().iterator();
	 while( it.hasNext() ) {
	    Object key = it.next();
	    ObjVector vec = (ObjVector) temp.get( key );
	    for ( int i = 0; i < vec.size(); i ++ ) {
	       out.println( key + " pp " + vec.elementAt( i ) );
	       RemoveInteraction( temp, (String) key, (String) vec.elementAt( i ), bothWays );
	       it = temp.keySet().iterator();
	    }
	 }
	 out.flush(); out.close();
      } catch( Exception e ) { e.printStackTrace(); }
   }

   public static boolean HasInteraction( Map map, String p1, String p2, boolean bothWays ) {
      p1 = replaceSynonym( p1 ).toUpperCase();
      p2 = replaceSynonym( p2 ).toUpperCase();
      boolean out = false;
      if ( map.get( p1 ) != null && ( (ObjVector) map.get( p1 ) ).indexOf( p2 ) >= 0 ) out = true;
      if ( bothWays && map.get( p2 ) != null && ( (ObjVector) map.get( p2 ) ).indexOf( p1 ) >= 0 ) out = true;
      return out;
   }

   public static boolean RemoveInteraction( Map map, String p1, String p2, boolean bothWays ) {
      p1 = replaceSynonym( p1 ).toUpperCase();
      p2 = replaceSynonym( p2 ).toUpperCase();
      boolean out = false;
      if ( map.get( p1 ) != null ) { ( (ObjVector) map.get( p1 ) ).removeElement( p2 ); out = true; }
      if ( bothWays && map.get( p2 ) != null ) { ( (ObjVector) map.get( p2 ) ).removeElement( p1 ); out = true; }
      return out;
   }

   public static void AddInteraction( Map map, String p1, String p2, boolean bothWays ) {
      p1 = replaceSynonym( p1 ).toUpperCase();
      p2 = replaceSynonym( p2 ).toUpperCase();
      if ( map.get( p1 ) != null ) ( (ObjVector) map.get( p1 ) ).addElement( p2 );
      else { ObjVector vec = new ObjVector(); vec.addElement( p2 ); map.put( p1, vec ); }
      if ( bothWays ) {
	 if ( map.get( p2 ) != null ) ( (ObjVector) map.get( p2 ) ).addElement( p1 );
	 else { ObjVector vec = new ObjVector(); vec.addElement( p1 ); map.put( p2, vec ); }
      }
   }

   public static Map GetIntersection( Map map1, Map map2, boolean bothWays ) {
      Map out = new HashMap();
      Iterator it1 = map1.keySet().iterator();
      while( it1.hasNext() ) {
	 String p1 = (String) it1.next();
	 String P1 = p1.toUpperCase();
	 ObjVector vec1 = (ObjVector) map1.get( p1 );
	 for ( int i = 0; i < vec1.size(); i ++ ) {
	    String p2 = (String) vec1.elementAt( i );
	    String P2 = p2.toUpperCase();

	    if ( ! HasInteraction( out, P1, P2, false ) ) {
	       if ( HasInteraction( map2, P1, P2, false ) ) 
		  AddInteraction( out, P1, P2, false );
	    }
	    if ( bothWays && ! HasInteraction( out, P2, P1, false ) ) {
	       if ( HasInteraction( map2, P2, P1, false ) )
		  AddInteraction( out, P2, P1, false );
	    }
	 }
      }
      System.out.println(out);
      return out;
   }

   protected static String replaceSynonym( String p1 ) {  // A hack to replace synonyms...
      if ( "YHR114W".equals( p1 ) ) return "BZZ1";
      //if ( "YHR016C".equals( p1 ) ) return "YSC84";
      if ( "YSC84".equals( p1 ) ) return "YHR016C";
      if ( "YBR098W".equals( p1 ) ) return "MMS4";
      //if ( "YHL027W".equals( p1 ) ) return "RIM101";
      return p1;
   }

   public static void main( String args[] ) {
      Map sif1 = ReadSif( args[ 0 ], false );
      Map sif2 = ReadSif( args[ 1 ], false );
      Map out = GetIntersection( sif1, sif2, false );
      WriteSif( out, "out.sif", false );
   }
}
