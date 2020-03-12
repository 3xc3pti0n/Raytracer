import java.io.BufferedReader;
import java.io.FileReader;
import java.io.PrintWriter;

public class stl_to_tri

{

	public static void main(String arg, PrintWriter ausgabe) {
		int counter_vertex = 0;
		int counter_triangle = 0;

		try {
			final BufferedReader br = new BufferedReader(new FileReader(arg));
			String line;
			while (((line = br.readLine()) != null) && (!line.trim().equals("END"))) {
				String bezeichnung = new String();
				bezeichnung = line;

				if (bezeichnung.contains("vertex")) {

					if (counter_vertex < 1) {
						ausgabe.println("TRIANGLE " + counter_triangle);
						counter_triangle += 1;
					}

					line = line.trim();

					//String[] splitresult = line.split(" ");
					String[] splitresult = line.split("\\s+");

					if (splitresult[1].isEmpty()) { // Falls ein Leerzeichen zuviel
						splitresult = line.split("  ");
						splitresult = splitresult[1].split(" ");
						// System.out.println(splitresult[0]);

						splitresult[0] = (Double.valueOf(splitresult[0]).toString());
						splitresult[1] = (Double.valueOf(splitresult[1]).toString());
						splitresult[2] = (Double.valueOf(splitresult[2]).toString());

						ausgabe.println(splitresult[0] + ",    " + splitresult[1] + ",    " + splitresult[2]);

						counter_vertex = counter_vertex + 1;

					}
					else {

						splitresult[1] = (Double.valueOf(splitresult[1]).toString());
						splitresult[2] = (Double.valueOf(splitresult[2]).toString());
						splitresult[3] = (Double.valueOf(splitresult[3]).toString());

						ausgabe.println(splitresult[1] + ",    " + splitresult[2] + ",    " + splitresult[3]);

						counter_vertex = counter_vertex + 1;

					}

				}
				else {
					counter_vertex = 0;
				}

			}
			br.close();
		}
		catch (final Exception e) {
			System.out.println(e.toString());
		}
	}
	
	public static void stl_to_detector(String arg, PrintWriter ausgabe) {
		int counter_vertex = 0;
		int counter_triangle = 0;
		try {
			final BufferedReader br = new BufferedReader(new FileReader(arg));
			String line;
			while (((line = br.readLine()) != null) && (!line.trim().equals("END"))) {
				String bezeichnung = new String();
				bezeichnung = line;

				if (bezeichnung.contains("vertex")) {

					if (counter_vertex < 1) {
						ausgabe.println("DETECT_TRIANGLE " + counter_triangle);
						counter_triangle += 1;
					}
					line = line.trim();
					String[] splitresult = line.split(" ");
					if (splitresult[1].isEmpty()) { // Falls ein Leerzeichen zuviel
						splitresult = line.split("  ");
						splitresult = splitresult[1].split(" ");
						splitresult[0] = (Double.valueOf(splitresult[0]).toString());
						splitresult[1] = (Double.valueOf(splitresult[1]).toString());
						splitresult[2] = (Double.valueOf(splitresult[2]).toString());
						ausgabe.println(splitresult[0] + ",    " + splitresult[1] + ",    " + splitresult[2]);
						counter_vertex = counter_vertex + 1;
					}
					else {
						splitresult[1] = (Double.valueOf(splitresult[1]).toString());
						splitresult[2] = (Double.valueOf(splitresult[2]).toString());
						splitresult[3] = (Double.valueOf(splitresult[3]).toString());
						ausgabe.println(splitresult[1] + ",    " + splitresult[2] + ",    " + splitresult[3]);
						counter_vertex = counter_vertex + 1;
					}
				}
				else {
					counter_vertex = 0;
				}
			}
			br.close();
		}
		catch (final Exception e) {
			System.out.println(e.toString());
		}
	}
}
