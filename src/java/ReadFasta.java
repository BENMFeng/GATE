import java.io.*;
import java.util.*;
import org.biojava.bio.*;
import org.biojava.bio.symbol.*;
import org.biojavax.SimpleNamespace;
import org.biojavax.bio.seq.*;

public class ReadFasta {
/**
* The program takes two args: the first is the file name of the Fasta file.
* The second is the name of the Alphabet. Acceptable names are 
* DNA RNA or PROTEIN.
*/
  public static void main(String[] args) throws
      FileNotFoundException, BioException {
    String filename = args[0];
    BufferedReader br = new BufferedReader(new FileReader(filename));
    Alphabet alpha = AlphabetManager.alphabetForName(args[1]);
    SimpleNamespace ns = new SimpleNamespace("biojava");

    RichSequenceIterator iterator = RichSequence.IOTools.readFasta(br,
            alpha.getTokenization("token"), ns);
    while (iterator.hasNext()) {
        RichSequence rec = iterator.nextRichSequence();
        System.out.println(rec.getName());
        System.out.println(rec.length());
    }
  }
}