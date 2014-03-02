package djr.util;

import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.util.jar.*;

/**
 * This utility makes a JAR runnable by inserting a <code>Main-Class</code>
 * attribute into the Manifest.
 *
 * @author Shawn Silverman
 */
public class MakeJarRunnable {
    public static void main(String[] args) throws IOException {
        if (args.length != 3) {
            System.out.println("Usage: MakeJarRunnable "
                               + "<jar file> <Main-Class> <output>");
            System.exit(0);
        }

        // Create the JarInputStream object, and get its Manifest

        JarInputStream jarIn = new JarInputStream(new FileInputStream(args[0]));
        Manifest manifest = jarIn.getManifest();
        if (manifest == null) {
            // This will happen if there is no Manifest

            manifest = new Manifest();
        }

        Attributes a = manifest.getMainAttributes();
        String oldMainClass = a.putValue("Main-Class", args[1]);

        // If there was an old value there, tell the user about it and exit

        if (oldMainClass != null) {
            System.out.println("Warning: old Main-Class value is: "
                               + oldMainClass);
            System.exit(1);
        }

        System.out.println("Writing to " + args[2] + "...");
        JarOutputStream jarOut = new JarOutputStream(new FileOutputStream(args[2]),
                                                     manifest);

        // Create a read buffer to be used for transferring data from the input

        byte[] buf = new byte[4096];

        // Iterate the entries

        JarEntry entry;
        while ((entry = jarIn.getNextJarEntry()) != null) {
            // Exclude the Manifest file from the old JAR

            if ("META-INF/MANIFEST.MF".equals(entry.getName())) continue;

            // Write out the entry to the output JAR

            jarOut.putNextEntry(entry);
            int read;
            while ((read = jarIn.read(buf)) != -1) {
                jarOut.write(buf, 0, read);
            }

            jarOut.closeEntry();
        }

        // Flush and close all the streams

        jarOut.flush();
        jarOut.close();

        jarIn.close();
    }
}
