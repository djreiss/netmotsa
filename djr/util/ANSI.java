package djr.util;

/**
 *
 * ANSI contains static methods and variables used to generate a commonly 
 * used subset of ANSI X3.64 terminal control escape sequences.  Does not 
 * check for valid arguments; invalid command sequences are generally 
 * ignored by the ANSI display driver.
 *
 * @author <a href="mailto:reiss@uw.edu">David Reiss</a>
 * @version 1.9978 (Fri Nov 07 05:56:26 PST 2003)
 */
public class ANSI {
   // attributes
   public static final String RESET		= "0";
   public static final String BOLD		= "1";
   public static final String UNDERSCORE       	= "4";
   public static final String BLINK		= "5";
   public static final String REVERSE 	= "7";
   public static final String INVISIBLE	= "8";

   // colors
   public static final String FGBLACK 	= "30";
   public static final String FGRED		= "31";
   public static final String FGGREEN 	= "32";
   public static final String FGYELLOW	= "33";
   public static final String FGBLUE		= "34";
   public static final String FGMAGENTA	= "35";
   public static final String FGCYAN		= "36";
   public static final String FGWHITE 	= "37";
   public static final String BGBLACK 	= "40";
   public static final String BGRED		= "41";
   public static final String BGGREEN 	= "42";
   public static final String BGYELLOW	= "43";
   public static final String BGBLUE		= "44";
   public static final String BGMAGENTA	= "45";
   public static final String BGCYAN		= "46";
   public static final String BGWHITE 	= "47";

   public static final String ESCAPE = "\u001b[";

   /**
      Sets display attributes.  Multiple attributes may be set in a single 
      command sequence by separating the attributes with semicolons, e.g., 
      <code>ANSI.FGRED + ';' + ANSI.BOLD</code> for bold red foreground text.  
      The ANSI X3.64 specification does not call for a trailing semicolon, and 
      most implementations will ignore the command sequence if one is present.
   */
   public static final String setAttr(String attr) {
      return(ESCAPE + attr + 'm');
   }

   /**
      Set cursor location, counting from the top left position (1,1)
   */
   public static final String setCursor(int row, int col) {
      return(ESCAPE + row + ';' + col + 'H');
   }

   /**
      Save the cursor location; see restoreCursor()
   */
   public static final String saveCursor() {
      return(ESCAPE + "s");
   }

   /**
      Restore previously saved cursor position; see saveCursor()
   */
   public static final String restoreCursor() {
      return(ESCAPE + "u");
   }

   /**
    * Move the cursor n positions to the left. If n is greater or
    * equal to the current cursor column, the cursor is moved to the
    * first column.
    */
   public static final String leftCursor(int n) {
      return ESCAPE + n + "D";
   }

   /**
    * Move the cursor n positions to the right. If n plus the current
    * cursor column is greater than the rightmost column, the cursor
    * is moved to the rightmost column.
    */
   public static final String rightCursor(int n) {
      return ESCAPE + n + "C";
   }

   /**
    * Move the cursor n rows up without changing the current column.
    * If n is greater than or equal to the current row, the cursor is
    * placed in the first row.
    */
   public static final String upCursor(int n) {
      return ESCAPE + n + "A";
   }

   /**
    * Move the cursor n rows down. If n plus the current row is greater
    * than the bottom row, the cursor is moved to the bottom row.
    */
   public static final String downCursor(int n) {
      return ESCAPE + n + "B";
   }

   /**
      Clear the screen and reset the cursor position
   */
   public static final String clearScreen() {
      return(ESCAPE + "2J");
   }


   /**
      Clear line from cursor position to end of line
   */
   public static final String clearLine() {
      return(ESCAPE + "K");
   }
} // end class ANSI
