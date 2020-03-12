import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;


public class ClosedLoop_utils {

	public static void write_T0(File t0, double lambda, double df) {
		// System.out.println(lambda);
		// System.out.println(df);

		String complete_out = "";
		boolean changed = false;
		if (t0.exists() && !t0.isDirectory()) {

			try {
				final FileReader FR = new FileReader(t0);
				final BufferedReader br = new BufferedReader(FR); // ASCII STL öffnen
				String line;
				while (((line = br.readLine()) != null) && (!line.trim().equals("END"))) {
					changed = false;
					if (line.contains("variables")) {
						final String[] containing = line.split("t=");
						System.out.println("variables \"df=" + String.valueOf(df/1000) + ";lambda=" + String.valueOf(lambda) + ";t=" + containing[1] + "\n");
						complete_out = complete_out + "variables \"df=" + String.valueOf(df / 1000) + ";lambda=" + String.valueOf(lambda) + ";t=" + containing[1] + "\n";
						// changed = true;
					}
					else {
						complete_out = complete_out + line + "\n";
					}

				}
				br.close();
				System.out.println("controlDict changed");
			}
			catch (final FileNotFoundException e) {
				System.out.println("l");
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
			catch (final IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}

			try (BufferedWriter bw = new BufferedWriter(new FileWriter(t0))) {
				bw.write(complete_out);
			}
			catch (final IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}

		}

	}
	
}
