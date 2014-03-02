package djr.util.gui;

import java.util.*;
import java.awt.*;
import java.awt.event.*;
import javax.swing.event.*;
import java.beans.*;
import javax.swing.*;

import gnu.getopt.*;
import corejava.*;

import ucar.util.prefs.*;
import ucar.util.prefs.ui.*;

import djr.util.*;
import djr.util.array.*;

/**
 * Abstract class <code>ArgProcessorGUI</code>.
 *
 * @author <a href="mailto:dreiss@systemsbiology.org">David Reiss</a>
 * @version 1.9978 (Fri Nov 07 05:56:26 PST 2003)
 */
public class ArgProcessorGUI extends ArgProcessor {
   protected String argv[], newArgs[];
   TabbedPrefsDialog dialog = null;
   ArgProcessorInterface obj = null;

   public boolean ProcessArgs( ArgProcessorInterface obj, String argv[] ) {
      this.argv = argv;
      this.obj = obj;
      try { 
	 SetupArgs( ReflectUtils.getFullClassName( obj ) );
	 SetupArgs( ReflectUtils.getFullClassName( obj ) + "_args" );
	 obj.SetupArgs( this );
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

	 if ( argv != null && argv.length > 0 ) {
	    Getopt g = new Getopt( argv[ 0 ], argv, "-", opts, true );
	    g.setOpterr( false ); // We'll do our own error handling
	    int c;
	    while ( ( c = g.getopt() ) != -1 ) {
	       String value = g.getOptarg();
	       int ind = g.getLongind();
	       if ( ind < 0 ) continue;
	       Object[] obs = (Object []) args.elementAt( opts[ ind ].getVal() );
	       if ( value == null || value.equals( "null" ) ) value = "true";
	       obs[ 2 ] = value;
	    }
	 }

	 dialog = new TabbedPrefsDialog( null, true, "Options for " + 
				    ReflectUtils.getFullClassName( obj ), null, this );
	 PrefPanel pp = null;
	 int nargs = 0;

	 for ( int i = 0, size = args.size(); i < size; i ++ ) {
	    if ( args.elementAt( i ) instanceof String ) {
	       String name = (String) args.elementAt( i );
	       if ( name.endsWith( ":" ) ) name = name.substring( 0, name.length()-1 );
	       pp = dialog.newTab( name );
	       nargs = 0;
	       continue;
	    }
	    if ( nargs ++ >= 10 ) {
	       pp.newColumn();
	       nargs = 0;
	    }
	    if ( pp == null ) pp = dialog.getPrefPanel();
	    Object[] obs = (Object []) args.elementAt( i );
	    String flag = (String) obs[ 0 ];
	    if ( "h".equals( flag ) || "?".equals( flag ) ) continue;
	    String type = (String) obs[ 1 ];
	    String deflt = (String) obs[ 2 ];
	    String desc = flag;
	    String toolTip = (String) obs[ 3 ];
	    if ( "<string>".equals( type ) ) {
	       pp.addTextField( desc, desc, deflt ).setToolTipText( toolTip );
	    } else if ( "<int>".equals( type ) ) {
	       try { 
		  pp.addIntField( desc, desc, Integer.parseInt( deflt ) ).setToolTipText( toolTip );
	       } catch( Exception e ) {
		  try { 
		     pp.addIntField( desc, desc, (int) (Long.
							valueOf( deflt )).longValue() ).setToolTipText( toolTip );
		  } catch( Exception ee ) { };
	       }
	    } else if ( "<float>".equals( type ) ) {
	       pp.addDoubleField( desc, desc, (Double.valueOf( deflt )).
				  doubleValue() ).setToolTipText( toolTip );
	    } else if ( type == null ) { 
	       pp.addCheckBoxField( desc, desc, (Boolean.valueOf( deflt )).
				    booleanValue() ).setToolTipText( toolTip );
	    } else if ( "<fname>".equals( type ) ) {
	       //pp.addTextField( desc, desc, "" + deflt ).setToolTipText( toolTip );
	       pp.addField( new FilenameField( desc, desc, deflt ) ).setToolTipText( toolTip );
	    } else if ( type.startsWith( "<irange:" ) ) {
	       StringTokenizer t = new StringTokenizer( type.substring( 8, type.length()-1 ), ":" );
	       pp.addField( new IntRangeField( desc, desc, 
					       Integer.parseInt( deflt ),
					       Integer.parseInt( t.nextToken() ),
					       Integer.parseInt( t.nextToken() ) ) ).setToolTipText( toolTip );
	    } else if ( type.startsWith( "<frange:" ) ) {
	       StringTokenizer t = new StringTokenizer( type.substring( 8, type.length()-1 ), ":" );
	       pp.addField( new DoubleRangeField( desc, desc,
						  ( Double.valueOf( deflt ) ).doubleValue(),
						  ( Double.valueOf( t.nextToken() ) ).doubleValue(),
						  ( Double.valueOf( t.nextToken() ) ).doubleValue() ) ).setToolTipText( toolTip );
	    } else if ( type.indexOf( '|' ) >= 0 ) {
	       StringTokenizer t = new StringTokenizer( type.substring( 1, type.length()-1 ), "|" );
	       ArrayList list = new ArrayList();
	       int ind = 0, dind = 0;
	       while( t.hasMoreElements() ) {
		  String tok = t.nextToken();
		  list.add( tok );
		  if ( tok.equals( deflt ) ) dind = ind;
		  ind ++;
	       }
	       pp.addTextComboField( desc, desc, list, dind ).setToolTipText( toolTip );
	    }
	 }
	 dialog.finish();
	 dialog.setLocation( 50, 50 );
	 dialog.show();

	 //if ( dialog.wasAccepted() ) {
	 // doAccept();
	    return true;
	    //} else return false;
      } catch( Exception e ) { 
	 e.printStackTrace();
      }
      return false;
   }

   protected void doAccept() {
      ObjVector newArgv = new ObjVector( argv );
      for ( int i = 0, size = args.size(); i < size; i ++ ) {
	 if ( args.elementAt( i ) instanceof String ) continue;
	 Object[] obs = (Object []) args.elementAt( i );
	 String flag = (String) obs[ 0 ];
	 if ( "h".equals( flag ) || "?".equals( flag ) ) continue;
	 Field fld = dialog.getField( flag );
	 String value = "";
	 if ( fld instanceof FilenameField ) value = ( (FilenameField) fld ).getText();
	 else if ( fld instanceof DoubleRangeField ) value = ( (DoubleRangeField) fld ).getDouble() + "";
	 else if ( fld instanceof IntRangeField ) value = ( (IntRangeField) fld ).getInt() + "";
	 else if ( fld instanceof Field.Text ) value = ( (Field.Text) fld ).getText();
	 else if ( fld instanceof Field.TextCombo ) value = ( (Field.TextCombo) fld ).getText();
	 else if ( fld instanceof Field.CheckBox ) value = ( (Field.CheckBox) fld ).isSelected() + "";
	 else if ( fld instanceof Field.Double ) value = ( (Field.Double) fld ).getDouble() + "";
	 else if ( fld instanceof Field.Int ) value = ( (Field.Int) fld ).getInt() + "";
	 if ( value.indexOf( ' ' ) >= 0 ) value = "'" + value + "'";
	 argsHash.put( flag, value );
	 if ( ! MyUtils.IsNullString( value ) ) finalArgs += "-" + flag + " " + value + " ";
	 else finalArgs += "-" + flag + " '' ";
	 String deflt = (String) obs[ 2 ];
	 if ( value != null && ! deflt.equals( value ) && ! "".equals( value ) ) {
	    if ( ! "true".equals( value ) && ! "false".equals( value ) ) 
	       newArgv.addElement( "-" + flag + " " + value );
	    else newArgv.addElement( "-" + flag );
	 }
      }
      
      newArgs = new String[ newArgv.size() ];
      for ( int i = 0, size = newArgv.size(); i < size; i ++ ) 
	 newArgs[ i ] = (String) newArgv.elementAt( i );

      if ( obj != null ) obj.SetArgs( this );
   }

   public String[] getNewArgs() { return newArgs; }

   public static class TabbedPrefsDialog extends JDialog implements ActionListener, 
								    java.beans.PropertyChangeListener {
      protected ObjVector prefPanels = null;
      protected PrefPanel pp;
      protected JTabbedPane tp;
      protected PreferencesExt substore = null;
      protected boolean accepted = false;
      protected ArgProcessorGUI proc = null;

      public TabbedPrefsDialog( RootPaneContainer parent, boolean modal, 
				String title, PreferencesExt store, ArgProcessorGUI proc ) {
	 super( ( parent != null ) && ( parent instanceof JFrame ) ? 
		(JFrame) parent : null );
	 this.proc = proc;
	 setModal( modal );
	 if ( title != null ) setTitle( title );
	 if ( store != null ) substore = (PreferencesExt) store.node( "Dialog" );

	 UIManager.addPropertyChangeListener( this );

	 tp = new JTabbedPane();
	 getContentPane().add( tp, BorderLayout.CENTER );
      }

      public PrefPanel newTab( String title ) {
	 pp = newPrefPanel( title );
	 tp.add( pp );
	 if ( prefPanels == null ) prefPanels = new ObjVector();
	 prefPanels.addElement( pp );
	 return pp;
      }

      protected PrefPanel newPrefPanel( String title ) {
	 pp = new PrefPanel( title, substore );
	 pp.addActionListener( this );

	 JButton dismiss = new JButton( "Close" );
	 dismiss.addActionListener( this );
	 pp.addButton( dismiss );

	 JButton cancel = new JButton( "Cancel" );
	 cancel.addActionListener( this );
	 pp.addButton( cancel );
	 return pp;
      }

      public void actionPerformed( ActionEvent e ) {
	 String cmd = e.getActionCommand();
	 if ( cmd.equals( "Close" ) ) accept();
	 else if ( cmd.equals( "Accept" ) ) proc.doAccept();
	 if ( cmd.equals( "Close" ) || cmd.equals( "Cancel" ) ) { setVisible( false ); dispose(); }
      }

      public void propertyChange( PropertyChangeEvent e) {
	 if ( e.getPropertyName().equals( "lookAndFeel" ) )
	    SwingUtilities.updateComponentTreeUI( TabbedPrefsDialog.this );
      }

      public PrefPanel getPrefPanel() { 
	 if ( pp == null ) pp = newPrefPanel( "" );
	 return pp; 
      }
      
      public void accept() {
	 /*if ( prefPanels != null ) {
	    for ( int i = 0, sz = prefPanels.size(); i < sz; i ++ ) {
	       PrefPanel p = (PrefPanel) prefPanels.elementAt( i );
	       p.accept();
	    }
	    } else pp.accept();*/
	 accepted = true;
      }

      public boolean wasAccepted() { return accepted; }

      public Field getField( String name ) {
	 Field out = pp.getField( name );
	 if ( out != null ) return out;
	 if ( prefPanels != null ) {
	    for ( int i = 0, sz = prefPanels.size(); i < sz; i ++ ) {
	       PrefPanel p = (PrefPanel) prefPanels.elementAt( i );
	       out = p.getField( name );
	       if ( out != null ) return out;
	    }
	 }
	 return null;
      }

      public void finish() {
	 if ( prefPanels != null ) {
	    for ( int i = 0, sz = prefPanels.size(); i < sz; i ++ ) {
	       PrefPanel p = (PrefPanel) prefPanels.elementAt( i );
	       p.finish();
	    }
	 } else {
	    getContentPane().add( pp );
	    pp.finish();
	 }

	 pack();
      }
   }

   public static class FilenameField extends Field.Text implements ActionListener {
      protected JPanel panel = null;
      protected JButton button = null;
      
      public FilenameField( String name, String label, String deflt ) {
	 super( name, label, deflt, null );
	 panel = new JPanel( new FlowLayout( FlowLayout.LEFT, 0, 0 ) );
	 button = new JButton( "..." );
	 button.setMargin( new Insets( -6, 8, 0, 8 ) );
	 button.addActionListener( this );
	 panel.add( tf );
	 panel.add( button );
      }

      public JComponent getEditComponent() { return panel; }      

      public void setToolTipText( String tip ) {
	 super.setToolTipText( tip );
	 tf.setToolTipText( tip );
	 panel.setToolTipText( tip );
	 button.setToolTipText( tip );
      }

      public void actionPerformed( ActionEvent e ) {
	 String currDir = System.getProperty( "user.dir" );
	 JFileChooser fd = new JFileChooser( currDir );
	 fd.setDialogTitle( "Choose your file:" );
	 fd.setDialogType( JFileChooser.OPEN_DIALOG );
	 int returnVal = fd.showOpenDialog( button ); 
	 if ( returnVal == JFileChooser.APPROVE_OPTION ) {
	    String file = fd.getSelectedFile().getAbsolutePath();
	    tf.setText( file );
	 }
      }
   }

   public static class IntRangeField extends Field.Text implements ChangeListener {
      protected JPanel panel = null;
      protected JLabel jlabel = null;
      protected JSlider slider = null;
      protected String format;
      
      public IntRangeField( String name, String label,
			    int deflt, int min, int max ) {
	 super( name, label, null );
	 value = getStoreValue( new Integer( deflt ) );
	 panel = new JPanel( new FlowLayout( FlowLayout.LEFT, 0, 0 ) );
	 int mwidth = ( "" + max ).length();
	 format = "%" + mwidth + "d";
	 jlabel = new JLabel( IntUtils.SPrintf( format, deflt ) );
	 slider = new JSlider( min, max, ( (Integer) value ).intValue() );
	 slider.addChangeListener( this );
	 panel.add( slider );
	 panel.add( jlabel );
      }

      public JComponent getEditComponent() { return panel; }

      protected Object getEditValue() { return new Integer( slider.getValue() ); }

      protected void setEditValue( Object value ) { 
	 slider.setValue( ( (Integer) value ).intValue() ); }

      public int getInt() { return ( (Integer) value ).intValue(); }

      public void setToolTipText( String tip ) {
	 jlabel.setToolTipText( tip );
	 panel.setToolTipText( tip );
	 slider.setToolTipText( tip );
      }

      public void stateChanged( ChangeEvent e ) {
	 accept(); next(); 
	 jlabel.setText( IntUtils.SPrintf( format, ( (JSlider) e.getSource() ).getValue() ) );
      }
   }

   public static class DoubleRangeField extends IntRangeField {
      public DoubleRangeField( String name, String label,
			       double deflt, double min, double max ) {
	 super( name, label, (int) (deflt*1000), (int) (min*1000), (int) (max*1000) );
	 int mwidth = ( "" + max ).length();
	 format = "%" + mwidth + ".3f";
	 jlabel.setText( DoubleUtils.SPrintf( format, 
					      ( (double) slider.getValue() ) / 1000.0 ) );
	 slider.removeChangeListener( slider.getChangeListeners()[ 0 ] );
	 slider.addChangeListener( this );
      }

      public void stateChanged( ChangeEvent e ) {
	 accept(); next(); 
	 jlabel.setText( DoubleUtils.SPrintf( format, ( (double) ( (JSlider) e.getSource() ).getValue() ) / 1000.0 ) );
      }

      protected Object getEditValue() { 
	 return new java.lang.Double( ( (double) slider.getValue() ) / 1000.0 ); }

      protected void setEditValue( Object value ) { 
	 slider.setValue( (int) ( (java.lang.Double) value ).doubleValue() * 1000 ); }

      public double getDouble() { return ( (java.lang.Double) value ).doubleValue(); }
   }
}

