/*
Raytracer zur Berechnung der absorbierten Intensität
(c) Andreas Michalowski
Autor: Andreas Michalowski
       Robert Bosch GmbH
	   CR/APJ2, Schwieberdingen
	   0711/811-43423
	   andreas.michalowski@de.bosch.com
 */

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
//import java.util.Random;

//import javax.swing.JFrame;

public class Raytracer {
	String[] argv;
	public static TriangleSet[] global_triSet;
	public ByteBuffer byte_buffer_stl_for_server;
	public ByteBuffer binary_intensities;

	public static double GesamtleistunRT;
	public static double AbsorbierteLeistungRT;

	public static TriangleSet[] get_global_triSet() {

		return global_triSet;
	}

	public void set_global_triSet(TriangleSet[] global_triSet) {

		Raytracer.global_triSet = global_triSet;

	}

	public static void main(String argv[]) {
		final Raytracer raytracer = new Raytracer(argv);
	}

	public Raytracer(String argv[]) {
		this.argv = argv;
		calculate(argv);
	}

	public void calculate(String argv[]) // Added FFE
	{
		System.out.print("neu...\n");
		if (argv.length < 3) {
			System.out.println("Error! Not enough arguments.");
			System.out.println("Usage:");
			System.out.println("java -Xmx512M Raytracer <output> <parameter-file> <input1> {<input2> <input3> ...}");
			System.exit(1);
		}
		final String filename_output = argv[0];
		final Parameters params = new Parameters(argv[1]);
		String logdaten = params.toString();
		logdaten = logdaten + "-------------------------------------------------------------\n";
		// Zufallszahlengenerator initialisieren
		// Random zufallsgenerator;
		// zufallsgenerator = new Random();
		// Minimale Distanz die ein Lichtstrahl gehen muss um erneut aufzutreffen (Wegen numerischen Fehlern ist 0 oft nur 1e-15)
		// double mintranslation = params.beam_min_distance;

		// Material aus dem die Dreiecke bestehen
		final Material material = new Material(params.material_IOR_real, params.material_IOR_imag, params.material_name);

		// Verloren gegangene Leistung
		double power_lost = 0;
		// Anzahl der Strahlen mit Leistung 0 die dann nicht berechnet werden
		int zerobeams = 0;
		// Eingestrahlte Leistung
		double power_in = 0;
		// Array welches die Anzahl an Reflexionen enthält
		final int[] number_of_reflections = new int[params.number_of_reflections];
		// Array welches die Anzahl an Reflexionen enthält
		final double[] abs_power_of_reflections = new double[params.number_of_reflections];
		// Dreiecke der bestrahlten Struktur einlesen
		final TriangleSet[] tset = new TriangleSet[argv.length - 2];

		for (int i = 0; i < tset.length; i++) {
			final String input = argv[i + 2];
			logdaten = logdaten + "Triangle Input: " + input + "\n";
			final Objectlist tl = leseDreieckeEin(input, material);
			tset[i] = new TriangleSet(tl, argv[i + 2]);
			if (i > 0) {
				tset[i].setType(1);
			}

		}

		final double[] extremeT = getMinMaxTriangleParameters(tset);
		System.out.println("Kuerzeste Dreiecksseite: " + extremeT[0]);
		System.out.println("Laengste Dreiecksseite : " + extremeT[1]);
		logdaten = logdaten + "Kuerzeste Dreiecksseite: " + extremeT[0] + "\n";
		logdaten = logdaten + "Laengste Dreiecksseite : " + extremeT[1] + "\n";

		// Grenzen festlegen
		final double xmin = extremeT[2] - params.bounding_box_spacing;
		final double xmax = extremeT[3] + params.bounding_box_spacing;
		final double ymin = extremeT[4] - params.bounding_box_spacing;
		final double ymax = extremeT[5] + params.bounding_box_spacing;
		final double zmin = extremeT[6] - params.bounding_box_spacing;
		final double zmax = extremeT[7] + params.bounding_box_spacing;
		System.out.println("Bounding Box:");
		System.out.println("xmin: " + xmin + " xmax: " + xmax);
		System.out.println("ymin: " + ymin + " ymax: " + ymax);
		System.out.println("zmin: " + zmin + " zmax: " + zmax);
		logdaten = logdaten + "Bounding Box:\n";
		logdaten = logdaten + "xmin: " + xmin + " xmax: " + xmax + "\n";
		logdaten = logdaten + "ymin: " + ymin + " ymax: " + ymax + "\n";
		logdaten = logdaten + "zmin: " + zmin + " zmax: " + zmax + "\n";

		// Quaderbäume erzeugen
		final Quader[] q = new Quader[params.number_of_threads];
		for (int i = 0; i < params.number_of_threads; i++) {
			System.out.println("Erzeuge Quader fuer Thread " + i);
			q[i] = new Quader(new double[] { (xmin + xmax) / 2, (ymin + ymax) / 2, (zmin + zmax) / 2 }, (xmax - xmin) / 2, (ymax - ymin) / 2, (zmax - zmin) / 2);
			q[i].generateSublevels(params.number_of_sublevels);

		}
		final Objectlist qobjlist = q[0].getQuaderList();
		final Object[] quaders = qobjlist.getArray();
		final double[] extremeQ = getMinMaxQuaderParameters(quaders);
		System.out.println("Kuerzeste Quaderseite: " + extremeQ[0]);
		System.out.println("Laengste Quaderseite : " + extremeQ[1]);
		logdaten = logdaten + "Kuerzeste Quaderseite: " + extremeQ[0] + "\n";
		logdaten = logdaten + "Laengste Quaderseite : " + extremeQ[1] + "\n";
		for (int i = 0; i < tset.length; i++) {
			logdaten = logdaten + "Anzahl der Dreiecke (" + tset[i].name + ") : " + tset[i].triangles.length + "\n";
		}
		// Dreiecke der Struktur und der Detektoren einsortieren
		for (int j = 0; j < params.number_of_threads; j++) {
			System.out.println("Sortiere Dreiecke in die Quader fuer Thread " + j);
			for (int i = 0; i < tset.length; i++) {
				q[j].relateTriangle(tset[i].triangles);
			}
		}
		final String outname = filename_output;
		System.out.println("Ausgabedateiname: " + outname);

		//
		final long startzeit = System.currentTimeMillis();
		// long letztezeit = startzeit;
		// int percent = 0;
		final int NumberOfBeams = params.number_of_beams;

		final RaytracingLoop[] rtLoops = new RaytracingLoop[params.number_of_threads];
		final Thread[] rtThreads = new Thread[params.number_of_threads];
		for (int threadnumber = 0; threadnumber < rtLoops.length; threadnumber++) {
			rtLoops[threadnumber] = new RaytracingLoop(NumberOfBeams / params.number_of_threads, params, q[threadnumber], this, threadnumber);
			rtThreads[threadnumber] = new Thread(rtLoops[threadnumber]);
			rtThreads[threadnumber].start();
		}
		try {
			synchronized (this) {
				for (int threadnumber = 0; threadnumber < rtLoops.length; threadnumber++) {
					wait();
				}
			}
		}
		catch (final Exception e) {
			System.out.println(e.toString());
		}
		for (int threadnumber = 0; threadnumber < rtLoops.length; threadnumber++) {
			power_in += rtLoops[threadnumber].power_in;
			power_lost += rtLoops[threadnumber].power_lost;
			zerobeams += rtLoops[threadnumber].zerobeams;

			// Array welches die Anzahl an Reflexionen enthält
			for (int i = 0; i < number_of_reflections.length; i++) {
				number_of_reflections[i] += rtLoops[threadnumber].number_of_reflections[i];
				abs_power_of_reflections[i] += rtLoops[threadnumber].abs_power_of_reflections[i];
			}
		}
		double Leistung = 0;
		double Faktor = 0;
		double E_p = Double.valueOf(params.beamPower) / Double.valueOf(params.pulse_repetition_rate);
		double tau_p = Double.valueOf(params.pulse_duration);
		// Optionale Normierung der Leistung auf gewünschten Wert: -------------------------------------------------------------------------------------------------------
		// Bei Simulation im gepulsten Betrieb wird in Fluenz statt Leistung gerechnet
		if (params.beamPower != 0) {
				System.out.println("mainwindow Power: " + params.beamPower);
				try {
					Faktor = 2*E_p / (tau_p * power_in); 
					power_in = power_in * Faktor;
					power_lost = power_lost * Faktor;
					for (int o = 0; o < tset[0].getNumber(); o++) {
						tset[0].triangles[o].setPower(tset[0].triangles[o].getPower() * Faktor);
					}
				}
				catch (final Exception ee) {
				}
			
			/*else {
				System.out.println("mainwindow Power: " + params.beamPower);
				try {
					Leistung = Double.valueOf(params.beamPower); // Eingabewert der Mittleren Leistung aus Parameter file
					Faktor = Leistung / power_in; 
					power_in = power_in * Faktor;
					power_lost = power_lost * Faktor;
					for (int o = 0; o < tset[0].getNumber(); o++) {
						tset[0].triangles[o].setPower(tset[0].triangles[o].getPower() * Faktor);
					}
				}
				catch (final Exception ee) {
				}
			}*/
		}

		final long endzeit = System.currentTimeMillis();
		final long dauer = endzeit - startzeit;
		System.out.println("Dauer der Rechnung: " + dauer + " Millisekunden.");
		logdaten = logdaten + "Dauer der Rechnung: " + dauer + " Millisekunden.\n";
		System.out.println(zerobeams + " Strahlen nicht berechnet, da deren Leistung Null.");
		logdaten = logdaten + zerobeams + " Strahlen nicht berechnet, da deren Leistung Null.\n";

		double normfaktor = 1; // Normierungsfaktor um die Intensität auf das Maximun der eingestrahlten Intensitätsverteilung (im Fokus) zu skalieren
		final double w0 = params.beam_waist_radius; // Strahlradius im Fokus
		if (params.beam_intensity_distribution.equals("gauss")) {
			normfaktor = 2 / (w0 * w0 * Math.PI) * NumberOfBeams / (params.beam_source_radius * params.beam_source_radius * Math.PI);
		}
		else if (params.beam_intensity_distribution.equals("ring")) {
			normfaktor = 4 / (w0 * w0 * Math.PI * Math.exp(1)) * NumberOfBeams / (params.beam_source_radius * params.beam_source_radius * Math.PI);
		}
		else if (params.beam_intensity_distribution.equals("tophat")) {
			normfaktor = 1 / (w0 * w0 * Math.PI) * NumberOfBeams / (params.beam_source_radius * params.beam_source_radius * Math.PI);
		}
		else if (params.beam_intensity_distribution.equals("planewave")) {
			if (params.beam_source_type.equals("circle")) {
				normfaktor = 1 / (Math.PI * Math.pow(params.beam_source_radius, 2)) * NumberOfBeams / (Math.PI * Math.pow(params.beam_source_radius, 2));
			}
			else if (params.beam_source_type.equals("rectangle")) {
				normfaktor = 1 / (params.beam_source_width_x * params.beam_source_width_y) * NumberOfBeams / (params.beam_source_width_x * params.beam_source_width_y);
			}
		}
		if (params.target_csv) {
			for (int i = 0; i < tset.length; i++) {
				try {
					System.out.println("csv..");
					// File datei = new File(outname+"."+tset[i].name+".output.csv");
					final File datei = new File(outname + "_output.csv");
					final FileWriter ausgabestrom = new FileWriter(datei);
					final PrintWriter ausgabe = new PrintWriter(ausgabestrom);
					//schreibeDreieckIntensitaetCSVDatei(tset[i], ausgabe, normfaktor);
					schreibeDreieckIntensitaet_SP_Datei(tset[i], ausgabe, normfaktor);
					try {
						ausgabe.close();
					}
					catch (final Exception e) {
						System.out.println(e.toString());
					}
				}
				catch (final Exception e) {
					System.out.println(e.toString());
				}
			}
		}
		if (params.target_ply_log) {
			for (int i = 0; i < tset.length; i++) {
				try {
					// File datei = new File(outname+"."+tset[i].name+".output_log.ply");
					final File datei = new File(outname + "_output_log.ply");
					final FileWriter ausgabestrom = new FileWriter(datei);
					final PrintWriter ausgabe = new PrintWriter(ausgabestrom);
					schreibeDreieckIntensitaetPLYDatei(tset[i], "LOG", ausgabe);
					try {
						ausgabe.close();
					}
					catch (final Exception e) {
						System.out.println(e.toString());
					}
				}
				catch (final Exception e) {
					System.out.println(e.toString());
				}
			}
		}
		if (params.target_ply_lin) {
			for (int i = 0; i < tset.length; i++) {
				try {
					// File datei = new File(outname+"."+tset[i].name+".output_lin.ply");
					final File datei = new File(outname + "_intensity_out.ply");
					final FileWriter ausgabestrom = new FileWriter(datei);
					final PrintWriter ausgabe = new PrintWriter(ausgabestrom);
					schreibeDreieckIntensitaetPLYDatei(tset[i], "FETZ", ausgabe);
					try {
						ausgabe.close();
					}
					catch (final Exception e) {
						System.out.println(e.toString());
					}
				}
				catch (final Exception e) {
					System.out.println(e.toString());
				}
			}/*
			for (int i = 0; i < tset.length; i++) {
				try {
					// File datei = new File(outname+"."+tset[i].name+".output_lin.ply");
					final File datei = new File(outname + "_output_lin.ply");
					final FileWriter ausgabestrom = new FileWriter(datei);
					final PrintWriter ausgabe = new PrintWriter(ausgabestrom);
					schreibeDreieckIntensitaetPLYDatei(tset[i], "LIN", ausgabe);
					try {
						ausgabe.close();
					}
					catch (final Exception e) {
						System.out.println(e.toString());
					}
				}
				catch (final Exception e) {
					System.out.println(e.toString());
				}
			}*/
		}

		// / Ausgabe für Server generieren

		for (int i = 0; i < tset.length; i++) {
			System.out.println("Raytracer:calculate, triangle_sets: " + String.valueOf(tset.length));
			byte_buffer_stl_for_server = binary_stl_with_intensities(tset[i], "LIN", true); // achtung bei mehreren Sets!, true--> modified STL
		}
		for (int i = 0; i < tset.length; i++) {
			System.out.println("Raytracer: generating double intensities, triangle_sets: " + String.valueOf(tset.length));

			binary_intensities = double_binary_intensities(tset[i], "LIN"); // achtung bei mehreren Sets!

		}

		System.out.println("Gesamtleistung eingestrahlt: " + power_in);
		logdaten = logdaten + "Gesamtleistung eingestrahlt: " + power_in + "\n";

		for (int i = 0; i < tset.length; i++) {
			final double powges = tset[i].getPower();

			System.out.println("Gesamtleistung absorbiert (" + tset[i].name + "): " + powges + "  " + (100 * powges / power_in) + "%");

			// Output an System_Out:

			logdaten = logdaten + "Gesamtleistung absorbiert (" + tset[i].name + "): " + powges + "  " + (100 * powges / power_in) + "%\n";
			GesamtleistunRT = power_in;
			AbsorbierteLeistungRT = powges;
		}
		System.out.println("Verlorene Leistung: " + power_lost);

		logdaten = logdaten + "Verlorene Leistung: " + power_lost + "\n";
		logdaten = logdaten + "-------------------------------------------------------------\n";
		logdaten = logdaten + "Statistics: Number of reflections\nNumber\tquantity\tpower\n";
		for (int i = 0; i < params.number_of_reflections; i++) {
			logdaten = logdaten + (i + 1) + "\t" + number_of_reflections[i] + "\t" + abs_power_of_reflections[i] + "\n";
		}
		try {
			final File datei = new File(outname + "_settings.log");
			final FileWriter ausgabestrom = new FileWriter(datei);
			final PrintWriter ausgabe = new PrintWriter(ausgabestrom);
			ausgabe.println(logdaten);
			try {
				ausgabe.close();
			}
			catch (final Exception e) {
				System.out.println(e.toString());
			}
		}
		catch (final Exception e) {
			System.out.println(e.toString());
		}

		// Schreibe Intensität in Boundary-File für Wärmeleitung --------PB
		/* ab hier wieder rein

		Algorithmen.calculateSurfaceTension(tset[0], 1.8);

		// -----------------------------------------------------------------------------------------------------------------------------------------------

		final double[][] Matrix_dreiecke = new double[tset[0].getNumber() * 3][3];
		double[][] corners = new double[3][3];
		int anz;
		anz = tset[0].getNumber();

		for (int i = 1; i < anz + 1; i++) {
			corners = tset[0].triangles[i - 1].getCorners();
			for (int k = 1; k < 4; k++) {
				for (int j = 0; j < 3; j++) {

					Matrix_dreiecke[i * k - 1][j] = corners[k - 1][j];

				}

			}
		}

		// Algorithmen.readFoamTGradient(tset[0],45); // Findet jetzt im Hauptprogramm statt
*/
		this.set_global_triSet(tset);

		System.out.println("Done_Tracing\n\n");
		System.out.println("Skalierung: " + Faktor + "Pulsenergie: " + E_p + "Pulsdauer: " + tau_p);

	}
	
	public static void schreibeDreieckIntensitaet_SP_Datei(TriangleSet ts, PrintWriter ausgabe, double nf) throws IOException {
		System.out.println("schreibe...csw");
		if (ts != null) {
			ausgabe.println("Tri Number, SP X [ mu_m ], SP Y [ mu_m ], SP Z [ mu_m ], N X [ mu_m ], N Y [ mu_m ], N Z [ mu_m ], K X [ mu_m ], K Y [ mu_m ], K Z [ mu_m ], Leistung [ W/mu_m**2 ]");
			nf = 1;
			System.out.println("schreibe...csw2");
			for (int i = 0; i < ts.getNumber(); i++) {
				final double [][] corners = ts.triangles[i].getCorners();
				final double [] tris_norm = ts.triangles[i].getNormalizedFacetNormal();
				final double [] sumVec_k = ts.triangles[i].getSumVec_k();
					double SPx = (corners[0][0] + corners[0][1] + corners[0][2])/3.;
					double SPy = (corners[1][0] + corners[1][1] + corners[1][2])/3.;
					double SPz = (corners[2][0] + corners[2][1] + corners[2][2])/3.;
					double normx = tris_norm[0];
					double normy = tris_norm[1];
					double normz = tris_norm[2];
					double kx = sumVec_k[0];
					double ky = sumVec_k[1];
					double kz = sumVec_k[2];
					ausgabe.println(i + "\t" + (float) SPx + "\t" + (float) SPy + "\t" + (float) SPz + "\t"  + (float) normx + "\t"  + (float) normy + "\t" + (float) normz + "\t" + (float) kx + "\t"  + (float) ky + "\t" + (float) kz + "\t" + (ts.triangles[i].getIntensity() / nf));
			}
		}
	}

	public static void schreibeDreieckIntensitaetCSVDatei(TriangleSet ts, PrintWriter ausgabe, double nf) throws IOException {
		System.out.println("schreibe...csw");
		if (ts != null) {
			ausgabe.println("[Name]");
			ausgabe.println("raytracer");
			ausgabe.println();
			ausgabe.println("[Data]");
			ausgabe.println("Node Number, X [ m ], Y [ m ], Z [ m ], Leistung [ K ]");
			
			nf = 1;
			System.out.println("schreibe...csw2");
			for (int i = 0; i < ts.getNumber(); i++) {
				final double[][] corners = ts.triangles[i].getCorners();
				for (int j = 0; j < 3; j++) {
					ausgabe.println((3 * i + j) + "; " + (float) corners[0][j] + ";  " + (float) corners[1][j] + ";  " + (float) corners[2][j] + ";  " + (ts.triangles[i].getIntensity() / nf));
				}
			}
			ausgabe.println();
			ausgabe.println("[Faces]");
			for (int i = 0; i < ts.getNumber(); i++) {
				ausgabe.println(3 * i + ", " + (3 * i + 1) + ", " + (3 * i + 2));
			}
		}
	}

	// / Added FFE ---------------------------
	public static ByteBuffer double_binary_intensities(TriangleSet ts, String param) { // generate Bytebuffer containing the intensities for each triangle
		final ByteBuffer ibuffer = ByteBuffer.allocate(ts.getNumber() * Double.BYTES).order(ByteOrder.LITTLE_ENDIAN); // generate header
		for (int i = 0; i < ts.getNumber(); i++) {
			ibuffer.putDouble(ts.triangles[i].getIntensity());
		}
		System.out.println("In Raytracer.double_binary_intensities:  n_triangles: " + String.valueOf(ts.getNumber()));
		ibuffer.flip();
		return ibuffer;
	}

	public static ByteBuffer binary_stl_with_intensities(TriangleSet ts, String param, boolean modified_stl) {
		if (!modified_stl) {
			double maxI = 0; // get Maximum Intensity
			for (int i = 0; i < ts.getNumber(); i++) {
				final double inten = ts.triangles[i].getIntensity();
				if (inten > maxI) {
					maxI = inten;
				}
			}

			System.out.println("Raytraceer: binary_stl_with_intensities");
			final ByteBuffer header_buffer = ByteBuffer.allocate(80).order(ByteOrder.LITTLE_ENDIAN); // generate header

			final StringBuilder headerBuilder = new StringBuilder();
			headerBuilder.append("Some header"); // text im Header
			final int restlength = 40 - headerBuilder.length();
			for (int i = 0; i < restlength; i++) {
				headerBuilder.append(" ");
			}

			final String finalHeader = headerBuilder.toString();
			System.out.println(finalHeader + String.valueOf(headerBuilder.length()));
			final char[] header_array = finalHeader.toCharArray();

			for (int i = 0; i < header_array.length; i++) {
				header_buffer.putChar(header_array[i]); // header in Byte array
			}
			System.out.println("header_buffer built");
			// Daten erst die Orte
			final List<Float> DataList = new ArrayList<Float>();
			// System.out.println("Number of triangles: " + String.valueOf(ts.getNumber()));
			for (int i = 0; i < ts.getNumber(); i++) {
				DataList.add((float) ts.triangles[i].getNormalizedFacetNormal()[0]); // Add normal
				DataList.add((float) ts.triangles[i].getNormalizedFacetNormal()[1]);
				DataList.add((float) ts.triangles[i].getNormalizedFacetNormal()[2]);

				final double[][] edge = ts.triangles[i].getCorners();
				for (int j = 0; j < 3; j++) { // Add all vertexex
					for (int k = 0; k < 3; k++) {
						DataList.add((float) edge[k][j] / 1000.0f);
					}
				}
			}
			System.out.println("Data list built");
			final ByteBuffer data_buffer = ByteBuffer.allocate(DataList.size() * Float.BYTES + 2 * DataList.size() / 12).order(ByteOrder.LITTLE_ENDIAN);
			for (int i = 0; i < DataList.size(); i++) {
				data_buffer.putFloat(DataList.get(i));
				if ((i + 1) % 12 == 0) {
					data_buffer.putShort((short) (ts.triangles[i / 12].getIntensity() * Short.MAX_VALUE / maxI)); // UINT16 – Attribute byte count, am schluss - hier die intensität
					// System.out.println("Number of triangle: " + String.valueOf((short) (ts.triangles[(int) i/12].getIntensity() * Short.MAX_VALUE / maxI)));
				}
			}
			System.out.println("data_buffer filled");
			final ByteBuffer final_file = ByteBuffer.allocate(header_buffer.capacity() + data_buffer.capacity() + Integer.BYTES).order(ByteOrder.LITTLE_ENDIAN);
			header_buffer.flip();
			data_buffer.flip();
			// UINT32 – Number of triangles, nach header
			final_file.put(header_buffer).putInt(ts.getNumber()).put(data_buffer); // header und daten kombinieren

			final_file.flip();
			header_buffer.flip();
			data_buffer.flip();
			System.out.println("Len OutBuffer : " + String.valueOf(final_file.capacity()));

			return final_file;
		}
		else {
			double maxI = 0; // get Maximum Intensity
			for (int i = 0; i < ts.getNumber(); i++) {
				final double inten = ts.triangles[i].getIntensity();
				if (inten > maxI) {
					maxI = inten;
				}
			}
			final int n_triangles = ts.getNumber();
			final ByteBuffer bf = ByteBuffer.allocate(4 + n_triangles * (12 * 4 + 2)).order(ByteOrder.LITTLE_ENDIAN);
			bf.putInt(n_triangles);
			for (int i = 0; i < ts.getNumber(); i++) {
				final double[] normal = ts.triangles[i].getNormalizedFacetNormal();
				// System.out.println(i);
				for (int k = 0; k < 3; k++) {
					bf.putFloat((float) normal[k] / 1000.0f);
				}
				final double[][] edge = ts.triangles[i].getCorners();
				for (int j = 0; j < 3; j++) { // Add all vertexex
					for (int k = 0; k < 3; k++) {
						bf.putFloat((float) edge[k][j] / 1000.0f);
					}
				}
				bf.putShort((short) (ts.triangles[i].getIntensity() * Short.MAX_VALUE / maxI));
			}
			System.out.println("In Raytracer.binary_stl_with_intensities:  n_triangles: " + String.valueOf(n_triangles));
			bf.flip();
			return bf;
		}
	}

	public static float[] get_temperatures_of_stl_elements(ByteBuffer binary_stl, boolean modified_STL) { // reads binatry STL as BufferedBytes and returns the 2BYTE fields for colors - (might be temperatures)
		final List<Short> temperature_list = new ArrayList<Short>();
		if (!modified_STL) {
			binary_stl.position(80 + Integer.BYTES);
			while (binary_stl.remaining() > 12 * Float.BYTES + 2) {
				binary_stl.position(binary_stl.position() + 12 * Float.BYTES);
				temperature_list.add(binary_stl.getShort());
			}
		}
		else {
			binary_stl.position(Integer.BYTES); // no header
			while (binary_stl.remaining() >= 9 * Float.BYTES + 2) {
				binary_stl.position(binary_stl.position() + 9 * Float.BYTES);
				temperature_list.add(binary_stl.getShort());
			}
		}
		final float[] temperatures = new float[temperature_list.size()];
		for (int i = 0; i < temperature_list.size(); i++) {
			temperatures[i] = temperature_list.get(i);
			// System.out.println("element "+ String.valueOf(i) +"SHORT field value: " +String.valueOf(temperatures[i]));
		}
		return temperatures;
	}

	// / END Added FFE ---------------------------

	public static void schreibeDreieckIntensitaetPLYDatei(TriangleSet ts, // FFE changed
			String param, PrintWriter ausgabe) throws IOException {
		System.out.println("schreibeDreieckIntensitaetPLYDatei " + String.valueOf(ts.getNumber()));
		double maxI = 0;
		for (int i = 0; i < ts.getNumber(); i++) {
			final double inten = ts.triangles[i].getIntensity();
			if (inten > maxI) {
				maxI = inten;
			}
		}
		double[] is = new double[ts.getNumber()];
		for (int i = 0; i < ts.getNumber(); i++) {
			is[i]= ts.triangles[i].getIntensity();
		}
		Arrays.sort(is);
		List<Double> intensities = new ArrayList<Double>();
		for (int i = 0; i < ts.getNumber(); i++) {
			if (is[i]> 0){ // nur wenn getroffen
				intensities.add(is[i]);
			}
		}
		int first_hit_i = 0;
		while (intensities.get(first_hit_i)==0){ first_hit_i++;}
		System.out.println("first_hit_i "+ first_hit_i);
		List<Double[]> cmap = Algorithmen.read_colormap("cmap_jet_ncar.txt");
		System.out.println("schreibe..ply");
		if (ts != null) {
			ausgabe.println("ply");
			ausgabe.println("format ascii 1.0");
			ausgabe.println("comment author: IFSW Michalowski Raytracer changed FFE/JH");
			ausgabe.println("comment object: ");
			ausgabe.println("element vertex " + 3 * ts.getNumber());
			ausgabe.println("property float x");
			ausgabe.println("property float y");
			ausgabe.println("property float z");
			ausgabe.println("property uchar red");
			ausgabe.println("property uchar green");
			ausgabe.println("property uchar blue");
			ausgabe.println("element face " + ts.getNumber());
			ausgabe.println("property list uchar int vertex_index");
			ausgabe.println("element edge 0");
			ausgabe.println("property int vertex1");
			ausgabe.println("property int vertex2");
			ausgabe.println("property uchar red");
			ausgabe.println("property uchar green");
			ausgabe.println("property uchar blue");
			ausgabe.println("end_header");
			System.out.println("schreibe..ply");
			int R = 0;
			int G = 0;
			int B = 0;
			int cmap_i = 0;
			for (int i = 0; i < ts.getNumber(); i++) {
				System.out.println(i);
				final double[][] corners = ts.triangles[i].getCorners();
				cmap_i = 0;
				for (int j = 0; j < 3; j++) {
					final double inten = ts.triangles[i].getIntensity();
					if (param == "LOG") {
						R = ((int) (Math.log(inten + 1) * 255.0 / Math.log(maxI + 1)));
						G = (255 - (int) (Math.log(inten) * 255.0 / Math.log(maxI)));
						B = 0;
					}
					if (param == "LIN") {
						R = ((int) (inten * 255.0 / maxI));
						G = (255 - (int) (inten * 255.0 / maxI));
						B = 0;
					}
					if (param == "FETZ") {
							//cmap_i = (int) (99*( inten/maxI));
							if (intensities.indexOf(inten) < first_hit_i) cmap_i = 0;
							else{ cmap_i = (int) (99*( (float)intensities.indexOf(inten)/intensities.size()));}
							//System.out.println("cmap_i "+ cmap_i);
							R = (int) ( cmap.get(cmap_i)[0]*255.);
							G = (int) ( cmap.get(cmap_i)[1]*255.);
							B = (int) ( cmap.get(cmap_i)[2]*255.);
					}
					ausgabe.println((float) (corners[0][j] / 1000.0f) + " " // FFE: nochmal durch 1000 geteilt
							+ (float) (corners[1][j] / 1000.0f) + " " + (float) (corners[2][j] / 1000.0f) + " " + R + " " + G + " " + B);
					System.out.println(i + " " + (float) (corners[0][j] / 1000.0f) + " " + (float) corners[0][j] / 1000.0f);
				}
			}
			for (int i = 0; i < ts.getNumber(); i++) {
				ausgabe.println("3 " + 3 * i + " " + (3 * i + 1) + " " + (3 * i + 2));
			}
		}
	}

	public static double[] getMinMaxTriangleParameters(TriangleSet[] ts) {
		final double[] e = ts[0].getMinMaxTriangleParameters();
		for (int i = 1; i < ts.length; i++) {
			final double[] t = ts[i].bounds;
			e[0] = t[0] < e[0] ? t[0] : e[0];
			e[1] = t[1] > e[1] ? t[1] : e[1];
			e[2] = t[2] < e[2] ? t[2] : e[2];
			e[3] = t[3] > e[3] ? t[3] : e[3];
			e[4] = t[4] < e[4] ? t[4] : e[4];
			e[5] = t[5] > e[5] ? t[5] : e[5];
			e[6] = t[6] < e[6] ? t[6] : e[6];
			e[7] = t[7] > e[7] ? t[7] : e[7];
		}
		return e;
	}

	public static double[] getMinMaxQuaderParameters(Object[] qs) {
		double minL = 1e10;
		double maxL = -1;
		for (int i = 0; i < qs.length; i++) {
			final double ax = ((Quader) qs[i]).ax;
			final double ay = ((Quader) qs[i]).ay;
			final double az = ((Quader) qs[i]).az;
			minL = ax < minL ? ax : minL;
			minL = ay < minL ? ay : minL;
			minL = az < minL ? az : minL;
			maxL = ax > maxL ? ax : maxL;
			maxL = ay > maxL ? ay : maxL;
			maxL = az > maxL ? az : maxL;
		}
		return new double[] { minL, maxL };
	}

	public static void adjust_tempdep_n_k() { //

	}

	public static Objectlist leseDreieckeEin(String datname, Material material) {
		final Objectlist ol = new Objectlist();
		double[] p1;
		double[] p2;
		double[] p3;
		// double[] nv_p1;
		// double[] nv_p2;
		// double[] nv_p3;
		// double n_re;
		// double n_imag;
		// int type;
		try {
			final BufferedReader br = new BufferedReader(new FileReader(datname));
			String line;
			System.out.println("\n\nleseDreieckeEin " + datname);
			while (((line = br.readLine()) != null) && (!line.trim().equals("END"))) {
				String bezeichnung = new String();
				bezeichnung = line;
				if (bezeichnung.startsWith("TRIANGLE")) {
					line = br.readLine();
					final String[] p1s = line.split(",");
					line = br.readLine();
					final String[] p2s = line.split(",");
					line = br.readLine();
					final String[] p3s = line.split(",");
					p1 = new double[3];
					p2 = new double[3];
					p3 = new double[3];
					for (int i = 0; i < 3; i++) {
						p1[i] = (new Double(p1s[i])).doubleValue() * 1000; 
						p2[i] = (new Double(p2s[i])).doubleValue() * 1000;
						p3[i] = (new Double(p3s[i])).doubleValue() * 1000;
						// Konvertieren von mm zu µm für Raytracing
					}
					ol.add(new Triangle(p1, p2, p3, material, bezeichnung));
				}
				else if (bezeichnung.startsWith("DETECT_TRIANGLE")) {
					line = br.readLine();
					final String[] p1s = line.split(",");
					line = br.readLine();
					final String[] p2s = line.split(",");
					line = br.readLine();
					final String[] p3s = line.split(",");
					p1 = new double[3];
					p2 = new double[3];
					p3 = new double[3];
					for (int i = 0; i < 3; i++) {
						p1[i] = (new Double(p1s[i])).doubleValue() * 1000; 
						p2[i] = (new Double(p2s[i])).doubleValue() * 1000;
						p3[i] = (new Double(p3s[i])).doubleValue() * 1000;
						// Konvertieren von mm zu µm für Raytracing
					}
					Triangle t = new Triangle(p1, p2, p3, material, bezeichnung);
					t.type = 1;
					ol.add(t);
				}
				else if (bezeichnung.startsWith("PERIODIC")) {
					line = br.readLine();
					final String[] p1s = line.split(",");
					line = br.readLine();
					final String[] p2s = line.split(",");
					line = br.readLine();
					final String[] p3s = line.split(",");
					line = br.readLine();
					final String[] trs = line.split(",");
					p1 = new double[3];
					p2 = new double[3];
					p3 = new double[3];
					final double[] tr = new double[3];
					for (int i = 0; i < 3; i++) {
						p1[i] = (new Double(p1s[i])).doubleValue();
						p2[i] = (new Double(p2s[i])).doubleValue();
						p3[i] = (new Double(p3s[i])).doubleValue();
						tr[i] = (new Double(trs[i])).doubleValue();
					}
					final Triangle t = new Triangle(p1, p2, p3, material, bezeichnung);
					t.setTranslationvector(tr[0], tr[1], tr[2]); // Vektor der
					// auf den
					// Auftreffpunkt
					// addiert
					// wird um
					// periodische
					// Randbedingungen
					// zu
					// bewirken.
					t.type = 2; // Periodische Randbedingungs Dreieck
					ol.add(t);

				}

				/*
				 * //if (p1s.length == 3) if (true) { ol.add(new Triangle(p1,
				 * p2, p3, material, bezeichnung)); } else if (p1s.length == 9)
				 * { nv_p1 = new double[3]; nv_p2 = new double[3]; nv_p3 = new
				 * double[3]; for (int i=0; i<3; i++) { nv_p1[i] = (new
				 * Double(p1s[i+3])).doubleValue(); nv_p2[i] = (new
				 * Double(p2s[i+3])).doubleValue(); nv_p3[i] = (new
				 * Double(p3s[i+3])).doubleValue(); } n_re = (new
				 * Double(p1s[6])).doubleValue(); n_imag = (new
				 * Double(p1s[7])).doubleValue(); p1s[8] = p1s[8].trim(); type =
				 * Integer.valueOf(p1s[8]).intValue(); ol.add(new Triangle(p1,
				 * p2, p3, nv_p1, nv_p2, nv_p3, type, new Material(n_re, n_imag,
				 * "n"+n_re+" k"+n_imag), bezeichnung)); }
				 */
				else {
					System.out.println("Inputformat für Dreiecke nicht erkannt.");
					System.out.println("Moeglichkeit 1:");
					System.out.println("Deieck1");
					System.out.println("p1x, p1y, p1z");
					System.out.println("p2x, p2y, p2z");
					System.out.println("p3x, p3y, p3z");
					System.out.println("usw.");
					System.out.println("Moeglichkeit 2:");
					System.out.println("Deieck1");
					System.out.println("p1x, p1y, p1z, nv1x, nv1y, nv1z, Re(IOR), Im(IOR), TYPE");
					System.out.println("p2x, p2y, p2z, nv2x, nv2y, nv2z, Re(IOR), Im(IOR), TYPE");
					System.out.println("p3x, p3y, p3z, nv3x, nv3y, nv3z, Re(IOR), Im(IOR), TYPE");
					System.out.println("usw.");
					System.exit(1);
				}
			}
			br.close();
		}
		catch (final Exception e) {
			System.out.println("Problem Methode: _lese Dreieck ein " + e.toString());
		}
		return ol;
	}
}
