package djr.util;

/**
 *
 * HTML contains static methods and variables used to generate a commonly 
 * used subset of HTML formatting sequences.
 *
 * @author <a href="mailto:reiss@uw.edu">David Reiss</a>
 * @version 1.9978 (Fri Nov 07 05:56:26 PST 2003)
 */
public class HTML {
   protected static boolean IN_TABLE = false;
   protected static boolean IN_FONT = false;
   public static int CELL_PADDING = 0;

   // colors
   public static final String TABLE_START = "<table border=0 frame=void cellspacing=0 cellpadding=" + CELL_PADDING + "><tr>";
   public static final String TABLE_END = "</tr></table>";
   public static final String NEW_ROW   = "</tr><tr>";
   public static final String BR 	= "<br>";
   public static final String HR 	= "<hr>";
   public static final String BOLD 	= "<b>";
   public static final String SPACE	= "&nbsp";

   public static final String FGBLACK 	= "black";
   public static final String FGRED	= "#FF5555";
   public static final String FGGREEN 	= "#55FF55";
   public static final String FGYELLOW	= "#FFFF55";
   public static final String FGBLUE    = "#5555FF";
   public static final String FGMAGENTA	= "#FF55FF";
   public static final String FGCYAN	= "#55FFFF";
   public static final String FGWHITE 	= "white";

   public static final String BGBLACK 	= "black";
   public static final String BGRED	= "#AA1111";
   public static final String BGGREEN 	= "#11AA11";
   public static final String BGYELLOW	= "#AAAA11";
   public static final String BGBLUE    = "#1111AA";
   public static final String BGMAGENTA	= "#AA11AA";
   public static final String BGCYAN	= "#11AAAA";
   public static final String BGWHITE 	= "white";

   public static final String begin( String title ) {
      return "<html>\n<head>\n<title>" + title + "</title>\n</head>\n"; }
   public static final String body( String bg ) {
      return "<body bgcolor=" + bg + ">\n"; }
   public static final String end() {
      String out = ""; if ( IN_TABLE ) { IN_TABLE = false; out += TABLE_END; }
      if ( IN_FONT ) { IN_FONT = false; out += "</font>"; }
      out += "</body>\n</html>\n"; return out; }
   public static final String font( String font ) {
      IN_FONT = true; return "<font face=" + font + ">"; }
   public static final String font( String font, int size ) {
      IN_FONT = true; return "<font face=" + font + " size=" + size + ">"; }
   public static final String name( String label ) { 
      return "<a name=\"" + label + "\">"; }
   public static final String bold( String s ) { 
      return "<b>" + s + "</b>"; }
   public static final String underscore( String s ) { 
      return "<u>" + s + "</u>"; }
   public static final String blink( String s ) { 
      return "<blink>" + s + "</blink>"; }
   public static final String colorStart( String c ) {
      IN_FONT = true; return "<font color=" + c + ">"; }
   public static final String colorStop() { return "</font>"; }
   public static final String fgcolor( String c, String s ) {
      return colorStart( c ) + s + "</font>"; }
   public static final String tableStart() {
      String out = "";
      if ( ! IN_TABLE ) { out += TABLE_START; IN_TABLE = true; } return out; }
   public static final String cellStart() {
      String out = tableStart(); out += "<td>"; return out; }
   public static final String cellStartA( String align ) {
      String out = tableStart(); out += "<td align = " + align + ">"; return out; }
   public static final String cellStart( String c ) {
      String out = tableStart(); IN_FONT = true;
      out += "<td bgcolor=" + c + "><font color=black>"; return out; }
   public static final String cellStart( String c, String align ) {
      String out = tableStart(); IN_FONT = true;
      out += "<td bgcolor=" + c + " align=" + align + "><font color=black>"; return out; }
   public static final String cellEnd() { return "</td>"; }
   public static final String bgcolor( String c, String s ) {
      return cellStart( c ) + s + "</font>"; }
   public static final String reset() {
      if ( IN_TABLE ) { IN_TABLE = false; return TABLE_END; } return ""; }

   /*public static void main( String argv[] ) {
      System.out.println( HTML.begin( "HTML TEST" ) );
      System.out.println( HTML.body( HTML.BGRED ) );
      System.out.println( HTML.font( "courier", 18 ) );
      System.out.println( HTML.fgcolor( HTML.FGBLUE, "THIS IS A TEST!!!<br>" ) );
      System.out.println( "Really it is...<br>" );
      System.out.println( HTML.bgcolor( HTML.BGBLUE, "WILL THIS WORK???" ) );
      System.out.println( HTML.bgcolor( HTML.BGGREEN, "?????????" ) + 
			  HTML.reset() );
      System.out.println( "I dunno." );
      System.out.println( HTML.end() );
   }
   */
}
