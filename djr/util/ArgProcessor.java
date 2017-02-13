package djr.util;

import java.util.*;
import gnu.getopt.*;
import corejava.*;

import djr.util.array.*;

/**
 * Abstract class <code>ArgProcessor</code>.
 *
 * @author <a href="mailto:reiss@uw.edu">David Reiss</a>
 * @version 1.9978 (Fri Nov 07 05:56:26 PST 2003)
 */
public class ArgProcessor implements java.io.Serializable {
   protected String finalArgs = "", finalParams = "";
   protected Map argsHash;
   protected ArgProcessorInterface obj;
   protected ObjVector args = null;

   public ArgProcessor() { };

   public ArgProcessor( String args[] ) {
      ProcessArgs( null, args );
   }

   public Map getArgsMap() { return argsHash; }
   public void putArg( String arg, Object obj ) { argsHash.put( arg, obj ); }
   public void setArg( String arg, Object obj ) { putArg( arg, obj ); }

   public String getArg( String par ) { 
      String out = (String) getObjectArg( par ); 
      if ( out == null ) 
	 System.err.println( "Warning: no argument for requested parameter " + par );
      return out; }

   public boolean hasArg( String par ) {
      StringTokenizer t = new StringTokenizer( par, "|" );
      while( t.hasMoreTokens() ) { String s = t.nextToken(); if ( argsHash.containsKey( s ) ) 
	 return true; } return false; }

   public Object getObjectArg( String par ) { 
      StringTokenizer t = new StringTokenizer( par, "|" ); 
      while( t.hasMoreTokens() ) { String s = t.nextToken(); if ( argsHash.containsKey( s ) ) 
	 return argsHash.get( s ); } return null; }

   public int getIntArg( String par ) { 
      StringTokenizer t = new StringTokenizer( par, "|" ); 
      while( t.hasMoreTokens() ) { String s = t.nextToken(); if ( argsHash.containsKey( s ) ) 
	 return (Integer.valueOf( getArg( s ) )).intValue(); } return 0; }

   public double getFloatArg( String par ) { 
      StringTokenizer t = new StringTokenizer( par, "|" ); 
      while( t.hasMoreTokens() ) { String s = t.nextToken(); if ( argsHash.containsKey( s ) ) 
      return (Double.valueOf( getArg( s ) )).doubleValue(); } return 0.0; }

   public double getDoubleArg( String par ) { 
      StringTokenizer t = new StringTokenizer( par, "|" ); 
      while( t.hasMoreTokens() ) { String s = t.nextToken(); if ( argsHash.containsKey( s ) ) 
      return getFloatArg( s ); } return 0.0; }

   public boolean getBooleanArg( String par ) { 
      StringTokenizer t = new StringTokenizer( par, "|" ); 
      while( t.hasMoreTokens() ) { String s = t.nextToken(); if ( argsHash.containsKey( s ) ) 
      return (Boolean.valueOf( getArg( s ) )).booleanValue(); } return false; }

   public short getShortArg( String par ) { 
      StringTokenizer t = new StringTokenizer( par, "|" ); 
      while( t.hasMoreTokens() ) { String s = t.nextToken(); if ( argsHash.containsKey( s ) ) 
      return (Short.valueOf( getArg( s ) )).shortValue(); } return 0; }

   public long getLongArg( String par ) { 
      StringTokenizer t = new StringTokenizer( par, "|" ); 
      while( t.hasMoreTokens() ) { String s = t.nextToken(); if ( argsHash.containsKey( s ) ) 
      return (Long.valueOf( getArg( s ) )).longValue(); } return 0L; }

   public String getArg( String par, String dflt ) {
      if ( hasArg( par ) ) return getArg( par ); else putArg( par, dflt ); return dflt; }
   public Object getObjectArg( String par, Object dflt ) {
      if ( hasArg( par ) ) return getObjectArg( par ); else putArg( par, "" + dflt ); return dflt; }
   public int getIntArg( String par, int dflt ) {
      if ( hasArg( par ) ) return getIntArg( par ); else putArg( par, "" + dflt ); return dflt; }
   public double getFloatArg( String par, float dflt ) {
      if ( hasArg( par ) ) return getFloatArg( par ); else putArg( par, "" + dflt ); return (double) dflt; }
   public double getDoubleArg( String par, double dflt ) {
      if ( hasArg( par ) ) return getDoubleArg( par ); else putArg( par, "" + dflt ); return dflt; }
   public boolean getBooleanArg( String par, boolean dflt ) {
      if ( hasArg( par ) ) return getBooleanArg( par ); else putArg( par, "" + dflt ); return dflt; }
   public short getShortArg( String par, short dflt ) {
      if ( hasArg( par ) ) return getShortArg( par ); else putArg( par, "" + dflt ); return dflt; }
   public long getLongArg( String par, long dflt ) {
      if ( hasArg( par ) ) return getLongArg( par ); else putArg( par, "" + dflt ); return dflt; }

   public void RemoveArg( ObjVector args, String arg ) {
      for ( int i = 0, size = args.size(); i < size; i ++ ) {
	 if ( args.elementAt( i ) instanceof String ) continue;
	 Object[] objs = (Object[]) args.elementAt( i );
	 String s = (String) objs[ 0 ];
	 if ( s.equals( arg ) ) { args.removeElementAt( i ); break; }
      }
   }

   public void ModifyDefaultArg( String arg, String newDefault ) {
      for ( int i = 0, size = args.size(); i < size; i ++ ) {
	 if ( args.elementAt( i ) instanceof String ) continue;
	 Object[] objs = (Object[]) args.elementAt( i );
	 String s = (String) objs[ 0 ];
	 if ( s.equals( arg ) ) { objs[ 2 ] = newDefault; break; }
      }
   }

   public void SetupArgs( String propsFile ) {
      if ( argsHash == null ) argsHash = new java.util.HashMap();
      if ( args == null ) args = new ObjVector();
      ResourceBundle props = null;
      try { props = ResourceBundle.getBundle( propsFile ); }
      catch( Exception e ) { props = null; }
      if ( props == null ) return;
      try {
	 String newPropsFile = MyUtils.ReplaceSubstring( propsFile, ".", "/" );
	 Enumeration e = MyUtils.ReadFileLines( newPropsFile + ".properties" ).elements();
	 while( e.hasMoreElements() ) {
	    String line = (String) e.nextElement();
	    if ( line == null || line.equals( "" ) || line.startsWith( "#" ) ) continue;
	    String toks[] = MyUtils.Tokenize( line, " =" );
	    if ( toks[ 0 ].startsWith( "args.include" ) ) {
	       SetupArgs( props.getString( toks[ 0 ] ) );
	    } else if ( toks[ 0 ].startsWith( "args.section" ) ) {
	       args.addElement( props.getString( toks[ 0 ] ) );
	    } else if ( toks[ 0 ].startsWith( "args.modifydefault" ) ) {
	       ModifyDefaultArg( toks[ 0 ].substring( toks[ 0 ].lastIndexOf( '.' ) + 1 ),
				 props.getString( toks[ 0 ] ) );
	    } else if ( toks[ 0 ].startsWith( "args.option" ) ) {
	       String opt = props.getString( toks[ 0 ] );
	       String type = props.getString( "args.type." + opt );
	       if ( "null".equals( type ) ) type = null;
	       String deflt = props.getString( "args.default." + opt );
	       String desc = props.getString( "args.description." + opt );
	       args.addElement( new Object[] { opt, type, deflt, desc } );
	    }
	 }
      } catch( Exception e ) { e.printStackTrace(); }
   }

   public void AddArg( String sectionHeader ) {
      if ( args == null ) args = new ObjVector();
      args.addElement( sectionHeader );
   }

   public void AddArg( String arg, String type, String dflt, String desc ) {
      if ( args == null ) args = new ObjVector();
      args.addElement( new Object[] { arg, type, dflt, desc } );
   }

   /**
    * Accept default args in the form of "arg.ARGNAME.default = value" from either a 
    * properties file or the System properties (given by "-Darg.ARGNAME.default=value")
    */
   public void ModifyDefaultArgs( String propsFile ) {
      if ( args == null ) args = new ObjVector();
      ResourceBundle props = null;
      try { ResourceBundle.getBundle( propsFile ); }
      catch( Exception e ) { 
	 System.err.println( "Warning: could not find arguments property file " + propsFile );
	 props = null; }
      if ( props == null ) return;
      for ( int i = 0, size = args.size(); i < size; i ++ ) {
	 if ( args.elementAt( i ) instanceof String ) continue;
	 Object[] objs = (Object[]) args.elementAt( i );
	 String s = (String) objs[ 0 ];
	 try {
	    String prop = props != null ? props.getString( "arg." + s + ".default" ) :
	       System.getProperty( "arg." + s + ".default" );
	    if ( prop != null ) objs[ 2 ] = prop;
	 } catch( Exception e ) { };
      }
   }

   protected String getRealFlag( String flag ) {
      if ( flag.startsWith( "--" ) ) flag = flag.substring( 2 );
      else if ( flag.startsWith( "-" ) ) flag = flag.substring( 1 );
      if ( flag.startsWith( "no-" ) ) flag = flag.substring( 3 );
      if ( flag.indexOf( '=' ) > 0 ) flag = flag.substring( 0, flag.indexOf( '=' ) );
      return flag;
   }

   public boolean ProcessArgs( ArgProcessorInterface obj, String argv[] ) {
      this.obj = obj;
      String usageString = null;
      if ( argv == null || argv.length == 0 ) usageString = "";
      for ( int i = 0; i < argv.length; i ++ ) finalArgs += argv[ i ] + " ";
      SetupArgs( obj != null ? ReflectUtils.getFullClassName( obj ) : "" );
      if ( obj != null ) {
	 SetupArgs( ReflectUtils.getFullClassName( obj ) + "_args" );
	 obj.SetupArgs( this );
      }
      int realSize = 0;
      for ( int i = 0, size = args.size(); i < size; i ++ ) 
	 if ( args.elementAt( i ) instanceof Object[] ) realSize ++;
      LongOpt opts[] = new LongOpt[ realSize ];
      for ( int i = 0, j = 0, size = args.size(); i < size; i ++ ) {
	 if ( args.elementAt( i ) instanceof String ) continue;
	 Object[] objs = (Object []) args.elementAt( i );
	 opts[ j ++ ] = new LongOpt( (String) objs[ 0 ], objs[1] == null ? 
				  LongOpt.NO_ARGUMENT : LongOpt.REQUIRED_ARGUMENT,
				  null, i );
      }
      if ( usageString == null ) {
	 Getopt g = new Getopt( argv[ 0 ], argv, "-", opts, true );
	 g.setOpterr( false ); // We'll do our own error handling
	 int c;
	 String flag = null;
	 while ( ( c = g.getopt() ) != -1 ) {	   
	    String value = g.getOptarg();
	    int ind = g.getLongind();
	    if ( ind < 0 ) {
	       if ( value == null ) {
		  String boolVal = "true";
		  flag = argv[ g.getOptind() - 1 ];
		  if ( flag.startsWith( "--" ) ) flag = flag.substring( 2 );
		  else if ( flag.startsWith( "-" ) ) flag = flag.substring( 1 );
		  if ( flag.startsWith( "no-" ) ) {
		     boolVal = "false"; flag = flag.substring( 3 ); }
		  if ( flag.indexOf( '=' ) > 0 ) {
		     boolVal = flag.substring( flag.indexOf( '=' ) + 1 );
		     flag = flag.substring( 0, flag.indexOf( '=' ) ); }
		  argsHash.put( flag, boolVal );
	       } else {
		  String curr = (String) argsHash.get( flag );
		  if ( curr != null && ! "true".equals( curr ) && ! "false".equals( curr ) ) 
		     argsHash.put( flag, curr + " " + value );
		  else argsHash.put( flag, value );
	       }
	       continue;
	    }
	    Object[] obs = (Object []) args.elementAt( opts[ ind ].getVal() );
	    flag = (String) obs[ 0 ];
	    if ( ind >= 0 ) {
	       if ( value != null && 
		    ! flag.equals( getRealFlag( argv[ g.getOptind() - 2 ] ) ) ) 
		  usageString = "ERROR: Invalid option: " + argv[ g.getOptind() - 2 ];
	       else if ( value == null && 
		    ! flag.equals( getRealFlag( argv[ g.getOptind() - 1 ] ) ) ) 
		  usageString = "ERROR: Invalid option: " + argv[ g.getOptind() - 1 ];
	    }
	    if ( "h".equals( flag ) || "?".equals( flag ) ) usageString = "";
	    if ( value == null || value.equals( "null" ) ) value = "true";
	    try { 
	       argsHash.put( flag, value );
	       obs[ 2 ] = value;
	    } catch( NumberFormatException e ) { 
	       System.err.println( "ERROR: value '" + value + "' for option '-" + flag + 
				   "' is invalid." );
	    }
	 }
      }
      if ( usageString != null ) {
	 PrintUsage( usageString );
	 return false;
      }
      for ( int i = 0, size = args.size(); i < size; i ++ ) {
	 if ( args.elementAt( i ) instanceof String ) continue;
	 Object[] obs = (Object []) args.elementAt( i );
	 String flag = (String) obs[ 0 ];
	 String value = (String) obs[ 2 ];
	 argsHash.put( flag, value );
	 if ( ! "h".equals( flag ) && ! "?".equals( flag ) ) {
	    if ( ! "".equals( value ) ) finalParams += "\t" + flag + " = " + value + "\n";
	    else finalParams += "\t" + flag + " = ''\n";
	 }
      }
      if ( obj != null ) obj.SetArgs( this );
      return true;
   }

   public String getFinalArgs() { return finalArgs; }
   public String getFinalParams() { return finalParams; }

   public String printArgs() {
      String out = "";
      Iterator i = getArgsMap().keySet().iterator();
      while( i.hasNext() ) {
	 String key = (String) i.next();
	 out += MyUtils.SPrintf( "\t%20s", key ) + 
	    MyUtils.SPrintf( "\t%40s\n", (String) getArgsMap().get( key ) );
      }
      return out;
   }

   public void PrintUsage( String usageString ) {
      if ( ! MyUtils.IsNullString( usageString ) ) { 
	 System.err.println( "\n" + usageString + "\n" );
	 MyUtils.Sleep( 3 );	
      }
      String className = obj != null ? ReflectUtils.getFullClassName( obj ) : "";
      System.err.println( "\n" + className + " version " + MyUtils.version );
      System.err.println( "David Reiss (ISB): reiss@uw.edu" );
      System.err.print( "USAGE: java " + className + " " );
      String output = "";
      int len = 0;
      for ( int i = 0, size = args.size(); i < size; i ++ ) {
	 if ( args.elementAt( i ) instanceof String ) {
	    System.err.print( "\n" + args.elementAt( i ) + "\n\t" );
	    len = 0;
	    continue;
	 }
	 Object[] objs = (Object []) args.elementAt( i );
	 output = "-" + objs[0] + " " + ( objs[1] != null ? (objs[1]+" ") : "" );
	 len += output.length();
	 if ( len > 80 ) { System.err.print( "\n\t" ); len = 0; }
	 System.err.print( output );
      }
      System.err.println( "\n" );
      for ( int i = 0, size = args.size(); i < size; i ++ ) {
	 if ( args.elementAt( i ) instanceof String ) {
	    System.err.print( "\n" + args.elementAt( i ) + "\n" );
	    continue;
	 }
	 Object[] objs = (Object []) args.elementAt( i );
	 Format.print( System.err, "%15s ", "-" + (String) objs[ 0 ] );
	 Format.print( System.err, "%15s ", (String) ( objs[1] != null ? objs[1] : "" ) );
	 System.err.print( objs[3] );
	 if ( objs[2] != null && ! "".equals( objs[2] ) )
	    System.err.println( " (default=" + objs[2] + ")" );
	 else System.err.println();
      }
   }
}
