package djr.util;

import java.lang.reflect.*;

/**
 * Class <code>ReflectUtils</code>
 *
 * @author <a href="mailto:dreiss@systemsbiology.org">David Reiss</a>
 * @version 1.9978 (Fri Nov 07 05:56:26 PST 2003)
 */
public class ReflectUtils {
   public static String getFullClassName( Object obj ) {
      return obj.getClass().getName();
   }

   public static String getClassName( Object obj ) {
      String name = getFullClassName( obj );
      return name.substring( name.lastIndexOf( '.' ) + 1 );
   }

   public static void setField( Object obj, String fieldName, Object value ) {
      try {
         Class c = obj.getClass();
         Field f = c.getField( fieldName );
         f.set( obj, value );
      } catch( Exception e ) { };
   }
   
   public static void setField( Object obj, String fieldName, double value ) {
      try {
         Class c = obj.getClass();
         Field f = c.getField( fieldName );
         f.setDouble( obj, value );
      } catch( Exception e ) { };
   }
   
   public static void setField( Object obj, String fieldName, int value ) {
      try {
         Class c = obj.getClass();
         Field f = c.getField( fieldName );
         f.setInt( obj, value );
      } catch( Exception e ) { };
   }
   
   public static void setField( Object obj, String fieldName, long value ) {
      try {
         Class c = obj.getClass();
         Field f = c.getField( fieldName );
         f.setLong( obj, value );
      } catch( Exception e ) { };
   }
   
   public static void setField( Object obj, String fieldName, short value ) {
      try {
         Class c = obj.getClass();
         Field f = c.getField( fieldName );
         f.setShort( obj, value );
      } catch( Exception e ) { };
   }
   
   public static void setField( Object obj, String fieldName, byte value ) {
      try {
         Class c = obj.getClass();
         Field f = c.getField( fieldName );
         f.setByte( obj, value );
      } catch( Exception e ) { };
   }
   
   public static void setField( Object obj, String fieldName, char value ) {
      try {
         Class c = obj.getClass();
         Field f = c.getField( fieldName );
         f.setChar( obj, value );
      } catch( Exception e ) { };
   }
   
   public static void setField( Object obj, String fieldName, boolean value ) {
      try {
         Class c = obj.getClass();
         Field f = c.getField( fieldName );
         f.setBoolean( obj, value );
      } catch( Exception e ) { };
   }
   
   public static Object getField( Object obj, String fieldName ) {
      try {
         Class c = obj.getClass();
         Field f = c.getField( fieldName );
         return f.get( obj );
      } catch( Exception e ) { };
      return null;
   }

   public static void shallowCopy( Object from, Object to ) {
      try {
	 Class cfrom = from.getClass(), cto = to.getClass();
	 Field fields[] = cfrom.getFields();
	 for ( int i = 0; i < fields.length; i ++ ) {
	    Field f2 = cto.getField( fields[ i ].getName() );
	    f2.set( to, fields[ i ].get( from ) );
	 }
      } catch( Exception e ) { e.printStackTrace(); }
   }

   public static final boolean ExecuteStaticMethod( String className, String method, Object[] args ) {
      Method meth = null;
      try {
	 Class c = Class.forName( className );
	 Class types[] = new Class[ args.length ];
	 for ( int i = 0; i < types.length; i ++ ) types[ i ] = args[ i ].getClass();
	 meth = c.getDeclaredMethod( method, types );
      } catch( Exception e ) { System.err.println( e ); e.printStackTrace(); }
      if ( meth != null ) {
	 try { meth.invoke( null, args ); return true; }
	 catch( Exception e ) { System.err.println( e ); e.printStackTrace(); }
      } else {
	 System.err.println( "ERROR: NO METHOD " + method + " FOUND IN CLASS " + className );
      }
      return false;
   }

   public static final boolean ExecuteMethod( Object obj, String method, Object[] args ) {
      Method meth = null;
      try {
	 Class c = obj.getClass();
	 Class types[] = new Class[ args.length ];
	 for ( int i = 0; i < types.length; i ++ ) types[ i ] = args[ i ].getClass();
	 meth = c.getDeclaredMethod( method, types );
      } catch( Exception e ) { System.err.println( e ); e.printStackTrace(); }
      if ( meth != null ) {
	 try { meth.invoke( obj, args ); return true; }
	 catch( Exception e ) { System.err.println( e ); e.printStackTrace(); }
      } else {
	 System.err.println( "ERROR: NO METHOD " + method + " FOUND IN CLASS " + obj.getClass().getName() );
      }
      return false;
   }

   public static final Object GetField( Object obj, String field ) {
      Field fld = null;
      try {
	 Class c = obj.getClass();
	 fld = c.getField( field );
      } catch( Exception e ) { System.err.println( e ); e.printStackTrace(); }
      if ( fld != null ) {
	 try { return fld.get( obj ); }
	 catch( Exception e ) { System.err.println( e ); e.printStackTrace(); }
      } else {
	 System.err.println( "ERROR: NO FIELD " + field + " FOUND IN CLASS " + 
			     obj.getClass().getName() );
      }
      return null;
   }

   public static final Object CallConstructor( String className, Object args[] ) {
      Constructor constr = null;
      try {
	 Class c = Class.forName( className );
	 Class types[] = new Class[ args.length ];
	 for ( int i = 0; i < types.length; i ++ ) types[ i ] = args[ i ].getClass();
	 constr = c.getDeclaredConstructor( types );
      } catch( Exception e ) { System.err.println( e ); e.printStackTrace(); }
      if ( constr != null ) {
	 try { return constr.newInstance( args ); }
	 catch( Exception e ) { System.err.println( e ); e.printStackTrace(); }
      } else {
	 System.err.println( "ERROR: REQUESTED CONSTRUCTOR NOT FOUND IN CLASS " + className );
      }
      return null;
   }
}
