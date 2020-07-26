import java.io.*;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.FilenameFilter;
import java.util.Date;


public class Main {
    public static void main(String[] args) {
        BufferedWriter bufferedWriter = null;
        File f = new File("out.txt");
        if (f.exists())
            f.delete();

        FilenameFilter filter = new FilenameFilter() {
        public boolean accept(File dir, String name) {
        return name.endsWith("_result");
            }
        };

        File folder = new File(".");
        File[] files = folder.listFiles(filter);
        try {
            //Construct the BufferedWriter object
            bufferedWriter = new BufferedWriter(new FileWriter("out.txt"));
            //Start writing to the output stream
            File file = new File("test.fa_result");
            BufferedReader reader = null;

            for(File currentfile : files)
  {
     try {
            reader = new BufferedReader(new FileReader(currentfile));
            String text = null;

            // repeat until all lines is read
            while ((text = reader.readLine()) != null) {
                bufferedWriter.write(text+" ");
            }
            bufferedWriter.write("\n");
        } catch (FileNotFoundException e) {
            e.printStackTrace();
        } catch (IOException e) {
            e.printStackTrace();
        } finally {
            try {
                if (reader != null) {
                    reader.close();
                }
            } catch (IOException e) {
                e.printStackTrace();
            }
  //      }

    }    // show file contents here
            }     //System.out.println(contents.toString());

} catch (FileNotFoundException ex) {
            ex.printStackTrace();
        } catch (IOException ex) {
            ex.printStackTrace();
        } finally {
            //Close the BufferedWriter
            try {
                if (bufferedWriter != null) {
                    bufferedWriter.flush();
                    bufferedWriter.close();
                }
            } catch (IOException ex) {
                ex.printStackTrace();
            }
        }
    }

}
