import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.Random;

public class RaytracingLoop implements Runnable {
	Random zufallsgenerator;
	Parameters params;
	Quader q;
	int NumberOfBeams;
	double power_lost;
	int zerobeams;
	double power_in;
	int[] number_of_reflections;
	double[] abs_power_of_reflections;
	Objectlist tlist;
	Raytracer rt;
	int threadnumber;

	public RaytracingLoop(int NumberOfBeams, Parameters params, Quader q, Raytracer rt, int tn) {
		zufallsgenerator = new Random();
		this.params = params;
		this.q = q;
		this.NumberOfBeams = NumberOfBeams;
		this.rt = rt;
		this.power_in = 0;
		this.power_lost = 0;
		this.zerobeams = 0;
		this.threadnumber = tn;
	}

	@Override
	public void run() {
		File beam_track_file = new File("./Raytracer_local_output/tracks.txt");
		int percent = 0;
		final long startzeit = System.currentTimeMillis();
		long letztezeit = startzeit;

		// Verloren gegangene Leistung
		power_lost = 0;
		// Anzahl der Strahlen mit Leistung 0 die dann nicht berechnet werden
		zerobeams = 0;
		// Eingestrahlte Leistung
		power_in = 0;
		// Array welches die Anzahl an Reflexionen enthält
		number_of_reflections = new int[params.number_of_reflections];
		// Array welches die Anzahl an Reflexionen enthält
		abs_power_of_reflections = new double[params.number_of_reflections];
		// Minimale Distanz die ein Lichtstrahl gehen muss um erneut aufzutreffen (Wegen numerischen Fehlern ist 0 oft nur 1e-15)
		final double mintranslation = params.beam_min_distance;
		double[][] Iprof = null;
		if (params.beam_intensity_distribution.equals("profile")) {
			Iprof = Ray_extend_tools.read_array_from_file("profile1.csv");
		}
		
		
		for (int strnr = 0; strnr < NumberOfBeams; strnr++) {
			boolean write_tracks = false;
			List<double[]> beam_tracks = new ArrayList<double[]>();
		    //Algorithmen.append_doubles_to_file(ray_log_file, starter, -1);// log strahlen
			if (((20 * strnr) % NumberOfBeams == 0) && (strnr > 0)) {
				System.out.println("Thread " + threadnumber + ":  " + percent + " % finished.");

				percent += 5;
				final long aktzeit = System.currentTimeMillis();
				final int restzeit = (int) ((20 - percent / 5) * (aktzeit - letztezeit) / 1000);
				System.out.println("Thread " + threadnumber + ":  Restzeit etwa " + restzeit + " Sekunden.");
				letztezeit = aktzeit;
			}
			// Lichtstrahl erzeugen
			// Zufallszahl fuer die x und y Koordinate erzeugen
			double zx = 0;
			double zy = 0;

			Ray strahl;

			// Falls Strahlquelle kreisförmig
			if (params.beam_source_type.equals("circle")) {
				final double KR = params.beam_source_radius; // Radius in dem Strahlen liegen dürfen
				boolean validbeam = false;
				while (!validbeam) {
					zx = 2 * KR * (zufallsgenerator.nextDouble() - 0.5);
					zy = 2 * KR * (zufallsgenerator.nextDouble() - 0.5);
					if ((Math.pow(zx, 2) + Math.pow(zy, 2)) < Math.pow(KR, 2)) {
						validbeam = true;
					}
				}
			}
			else if (params.beam_source_type.equals("rectangle")) {
				zx = params.beam_source_width_x * (zufallsgenerator.nextDouble() - 0.5);
				zy = params.beam_source_width_y * (zufallsgenerator.nextDouble() - 0.5);
			}
			else {
				System.out.println("No beam_source_type defined!");
				System.exit(1);
			}

			final double[] ov = params.incoming_beam_position_vector;
			final double[] rv = params.incoming_beam_direction_vector;
			final double w0 = params.beam_waist_radius;
			//final double ks = Math.sqrt(zx * zx + zy * zy) / w0;  // evtl. falsch!!
			final double[] nzs = new double[] { zx, zy, 0 };
			final double ztt = params.beam_distance_to_focus;
			final double M2 = params.beam_M2; // M2 des Strahls
			final double lambda = params.beam_wavelength; // Wellenlaenge des Strahls
			final double z0 = Math.PI * w0 * w0 / (lambda * M2);// 0.5;
			
			// ffetzer zur Korrektur  der Kaustik
			final double w_origin = w0*Math.sqrt(1 + ztt*ztt/(z0*z0) );
			final double ks = Math.sqrt(zx*zx + zy*zy)/w_origin;
			// ffetzer end korr
			
			//System.out.println("z0: "+z0);
			//System.out.println("ks: "+ks);
			double[] p1 = new double[] { 1, 0, 0 };
			double[] p2 = new double[] { 1, 0, 0 };
			if (params.polarization.equals("radial")) {
				p1 = new double[] { zx, zy, 0 };
				p2 = new double[] { zx, zy, 0 };// Radiale Polarisation
			}
			else if (params.polarization.equals("azimutal")) {
				p1 = new double[] { -zy, zx, 0 };
				p2 = new double[] { -zy, zx, 0 };// Azimutale Polarisation
			}
			else if (params.polarization.equals("linear_x")) {
				p1 = new double[] { 1, 0, 0 };
				p2 = new double[] { 1, 0, 0 };// Lineare Polarisation x
			}
			else if (params.polarization.equals("linear_y")) {
				p1 = new double[] { 0, 1, 0 };
				p2 = new double[] { 0, 1, 0 };// Lineare Polarisation y
			}
			else if (params.polarization.equals("zirkular")) {
				p1 = new double[] { 1, 0, 0 };
				p2 = new double[] { 0, 1, 0 };// Zirkulare Polarisation
			}
			else {
				System.out.println("Keine Polarisation wurde angegeben!");
				System.exit(1);
			}

			final double r = Math.sqrt(zx * zx + zy * zy);
			final double w2z = w0 * w0 * (1 + (ztt * ztt) / (z0 * z0)); // w(z)quadrat

			double I = 0.0;
			boolean kaustik = false;
			if (params.beam_intensity_distribution.equals("gauss")) {
				I = (2.0 / (w2z * Math.PI)) * Math.exp(-2 * r * r / w2z);// Gauss Grundmode
				kaustik = true;
			}
			else if (params.beam_intensity_distribution.equals("profile")) {
				I = Ray_extend_tools.lin_interpolate2d(zx, zy, Iprof);// profile1.csv
				kaustik = true;
			}
			else if (params.beam_intensity_distribution.equals("ring")) {
				I = 16 * r * r / (Math.PI * w2z * w2z) * Math.exp(-4 * r * r / w2z);// Ringmode
				kaustik = true;
			}
			else if (params.beam_intensity_distribution.equals("tophat")) {
				I = (r * r < w2z) ? 1.0 / (w2z * Math.PI) : 0.0; // Flat top(Hutprofil)
				kaustik = true;
			}
			else if (params.beam_intensity_distribution.equals("planewave")) {
				if (params.beam_source_type.equals("circle")) {
					I = 1.0 / (Math.PI * Math.pow(params.beam_source_radius, 2));
				}
				else if (params.beam_source_type.equals("rectangle")) {
					I = 1.0 / (params.beam_source_width_x * params.beam_source_width_y);
				}
				kaustik = false;
			}
			else {
				System.out.println("Kein Intensitaetsprofil wurde angegeben!");
				System.out.println(params.polarization);
				System.exit(1);
			}

			power_in += I;

			if (kaustik == true) {
				strahl = new Ray(ov, rv, ks, w0, nzs, ztt, z0, I, p1, p2, "ray");
		        //Algorithmen.append_doubles_to_file(ray_log_file, start_point, 0); // ray log
			}
			else {
				strahl = new Ray(new double[] { ov[0] + zx, ov[1] + zy, ov[2] }, rv, p1, p2, "ray");
				strahl.setPower(I);
				
			}

			if ((Algorithmen.norm2(strahl.p1) == 0) || (Algorithmen.norm2(strahl.p2) == 0)) {
				zerobeams++;
			}
			else {
				double[] start_point = new double[]{ov[0]+zx,ov[1]+zy,ov[2], strahl.I};// ray log
				beam_tracks.add(start_point);
				for (int anz_reflexe = 0; anz_reflexe < params.number_of_reflections; anz_reflexe++) {
					// Liste von Dreiecken erzeugen
					tlist = new Objectlist();
					// Aus allen geschnittenen Quadern die Dreiecke einer Liste anhaengen
					tlist = q.getQuaderIntersectionTriangleList_noduplicates(strahl);

					final Object[] tris = tlist.getArray();
					// System.out.println("Anzahl der Kandidaten: "+tris.length);
					Triangle target = null;
					double[] schnittpunkt = new double[3];
					double xnear = -1;
					double[] inbase = new double[3];
					// Alle moeglicherweise getroffenen Dreiecke durchgehen

					if (tris == null) {
						System.out.println(strahl.toString());
					}
					else {
						for (int i = 0; i < tris.length; i++) {
							final Triangle T = (Triangle) tris[i];
							/**
							 * Ab hier wird der Schnitt aus Gerade mit Dreiecksebene berechnet Hier muss eine Variante mit Strahl mit Kaustik ergaenzt werden.
							 */

							if (strahl.kaustik == false) {
								// Bestimme das Produkt aus Richtungsvektor und Normalenvektor der Ebene. Ist das Null, so kann es keinen Schnitt geben

								final double tnv = Algorithmen.product(strahl.rv, T.nv);
								if (tnv == 0) {
									continue;
								}
								;
								// Welchen Parameter muss man in die Gerade bis zum Schnittpunkt einsetzen?
								final double x = Algorithmen.paramx_intersection_plane_line(T.p1, T.nv, strahl.ov, Algorithmen.product(strahl.rv, T.nv));
								// Nur wenn es geradeaus ist (x>1e-10 bzw. mintranslation) gibt es einen Schnitt, die untere Grenze mintranslation kommt daher, weil
								// wegen Rundungsfehlern das selbe Dreieck evtl. gleich zweimal hintereinander geschnitten wird

								if (x > mintranslation) {
									// Schnittpunkt der Geraden mit der Dreiecksebene bestimmen
									final double[] p = Algorithmen.point(strahl.ov, strahl.rv, x);

									// Vom Schnittpunkt die Basis des Dreiecks abziehen
									final double[] p_t = new double[3];
									for (int j = 0; j < 3; j++) {
										p_t[j] = p[j] - T.p1[j];
									}
									// Nur wenn der Schnittpunkt innerhalb der Dreiecksebene liegt, Test 1
									inbase = Algorithmen.product(T.transform, p_t);
									// System.out.println("Entwicklungskoeffizienten\n"+inbase[0]+" "+inbase[1]+" "+inbase[2]+"\n");
									if ((inbase[0] >= 0) && (inbase[1] >= 0)) {
										// Nur wenn der Schnittpunkt innerhalb der Dreiecksebene liegt, Test 2
										final double inout = inbase[0] + inbase[1];
										// System.out.println("inout "+inout);
										if (inout <= 1) {
											// Falls es das erste Dreieck ist, ist es zunaechst mal das naechste
											if ((xnear < 0) || (x < xnear)) {
												xnear = x;
												target = T;
												schnittpunkt = p;

											}
										}
									}
								}
							}
							else {
								// Welchen Parameter muss man in die Gerade bis zum Schnittpunkt, wenn es einen gibt, einsetzen?
								final double[] xparams = Algorithmen.paramx_intersection_plane_gaussianray(T.p1, T.nv, strahl);

								// System.out.println("Xpara1 "+xparams[0]+" xpara2 "+xparams[1]+" xpara3 "+xparams[2]); // Hier liegt der Fehler, Methode mit solvePolynom funktioniert

								// Kein Schnittpunkt vorhanden
								if (xparams[0] == 0) {
									continue;
								}
								// Zwei Schnittpunkte vorhanden, sortiere nach Distanz
								if (xparams[0] == 2) { // Zwei Schnittpunkte vorhanden, sortiere nach Distanz
									if (xparams[1] > xparams[2]) {
										final double merk = xparams[2];
										xparams[2] = xparams[1];
										xparams[1] = merk;
									}
								}

								// Nur wenn es geradeaus ist (x>1e-10) gibt es einen Schnitt, die untere Grenze 1e-10 kommt daher, weil
								// wegen Rundungsfehlern das selbe Dreieck evtl. gleich zweimal hintereinander geschnitten wird
								for (int k = 0; k < (int) xparams[0]; k++) {

									if (xparams[k + 1] > mintranslation) {
										// Schnittpunkt der Geraden mit der Dreiecksebene bestimmen
										final double[] p = Algorithmen.point(strahl, xparams[k + 1]);
										// Vom Schnittpunkt die Basis des Dreiecks abziehen
										final double[] p_t = new double[3];
										for (int j = 0; j < 3; j++) {
											p_t[j] = p[j] - T.p1[j];
										}
										// Nur wenn der Schnittpunkt innerhalb der Dreiecksebene liegt, Test 1
										inbase = Algorithmen.product(T.transform, p_t);
										// System.out.println("Entwicklungskoeffizienten\n"+inbase[0]+" "+inbase[1]+" "+inbase[2]+"\n");
										if ((inbase[0] >= 0) && (inbase[1] >= 0)) {
											// Nur wenn der Schnittpunkt innerhalb der Dreiecksebene liegt, Test 2
											final double inout = inbase[0] + inbase[1];
											// System.out.println("inout "+inout);

											if (inout <= 1) {
												// Falls es das erste Dreieck ist, ist es zunaechst mal das naechste
												if (((xnear < 0) || (xparams[k + 1] < xnear))) {
													xnear = xparams[k + 1];
													target = T;
													schnittpunkt = p;
												}
											}
										}
									}
								}
							}
						}
					}
					// Diese Routine wird aufgerufen, wenn feststeht welches Dreieck getroffen wurde
					// Bestimme die neue Richtung des Strahls
					// Erhoehe die absorbierte Intensitaet des Dreiecks
					// Verringere die Intensitaet des Strahls um die absorbierte Intensitaet
					if (target != null) {

						number_of_reflections[anz_reflexe]++;
						// Falls es ein Strahl mit Kaustik ist der auftrifft, umwandeln in geraden Lichtstrahl
						if (strahl.kaustik == true) {
							strahl.rv = Algorithmen.direction(strahl, xnear);
							final double p1_length = Algorithmen.norm2(strahl.p1);
							final double p2_length = Algorithmen.norm2(strahl.p2);
							if ((Double.isNaN(p1_length)) || (p1_length == 0)) {
								System.out.println("Problem! " + p1_length);
							}
							strahl.p1 = Algorithmen.getNormTox(Algorithmen.perp_to_v1_in_v1v2_plane(strahl.rv, strahl.p1), p1_length);
							strahl.p2 = Algorithmen.getNormTox(Algorithmen.perp_to_v1_in_v1v2_plane(strahl.rv, strahl.p2), p2_length);
							strahl.kaustik = false;
						}
						strahl.ov = schnittpunkt;
						double[]  hitpoint = new double[] {schnittpunkt[0], schnittpunkt[1],schnittpunkt[2], strahl.I};// ray log
						beam_tracks.add(hitpoint);
						//Algorithmen.append_doubles_to_file(ray_log_file, schnittpunkt, strahl.I); // ray log

						//
						// Hier wird der Reflektionsstrahl berechnet, target.nv ist der Normalenvektor der Ebene
						// inbase[0] ist der Entwicklungskoeffizient des Vektors von p1 hin zu p2
						// inbase[1] ist der Entwicklungskoeffizient des Vektors von p1 hin zu p3
						//
						double[] md;
						if (target.interpolation) {
							inbase = Algorithmen.product(target.transform, schnittpunkt);
							double[] nv_eff = new double[3];
							for (int i = 0; i < 3; i++) {
								nv_eff[i] = target.nv_p1[i] + inbase[0] * (target.nv_p2[i] - target.nv_p1[i]) + inbase[1] * (target.nv_p3[i] - target.nv_p1[i]);
							}
							nv_eff = Algorithmen.getNormTox(nv_eff, 1);
							md = Algorithmen.getMirrorDirection(nv_eff, strahl.rv);
							// Hier kann man Test tun
						}
						else {
							md = Algorithmen.getMirrorDirection(target.nv, strahl.rv);
						}

						// Falls das getroffene Dreieck ein normales Dreieck ist, wird absorbiert und reflektiert
						if (target.type == 0) {
							// Berechnung des Winkels unter dem das Dreieck getroffen wird, der Winkel ist zwischen
							// der Flaechennormalen und dem Strahl definiert
							final double auftreffwinkel = 0.5 * (Math.PI - Math.acos(Algorithmen.product(strahl.rv, md) / (Algorithmen.norm2(strahl.rv) * Algorithmen.norm2(md))));
							final double winkelingrad = auftreffwinkel * 180.0 / Math.PI;
							final double[] rsrp = target.material.getrsrp(auftreffwinkel);
							double[] s;
							if (auftreffwinkel == 0) {
								// ausgabe.print(" Null Grad");
								s = Algorithmen.getNormalVector(strahl.rv); // Senkrechter Einfall
							}
							else {
								s = Algorithmen.xproduct(strahl.rv, md);// Richtung senkrechte Polarisation
							}
							s = Algorithmen.getNormTox(s, 1);
							double[] pe = Algorithmen.xproduct(strahl.rv, s);// Parallele Richtung des einfallenden Strahles
							pe = Algorithmen.getNormTox(pe, 1);
							double[] pr = Algorithmen.xproduct(md, s);// Parallele Richtung des reflektierten Strahles
							pr = Algorithmen.getNormTox(pr, 1);

							final double power_vor = strahl.I;
							strahl.calculateIntensity();

							final double detN = s[0] * pe[1] - s[1] * pe[0];
							double a = (strahl.p1[0] * pe[1] - strahl.p1[1] * pe[0]) / detN;
							double b = (strahl.p1[1] * s[0] - strahl.p1[0] * s[1]) / detN;

							strahl.p1 = Algorithmen.summe(Algorithmen.product(a * rsrp[0], s), Algorithmen.product(b * rsrp[1], pr));
							a = (strahl.p2[0] * pe[1] - strahl.p2[1] * pe[0]) / detN;
							b = (strahl.p2[1] * s[0] - strahl.p2[0] * s[1]) / detN;
							strahl.p2 = Algorithmen.summe(Algorithmen.product(a * rsrp[0], s), Algorithmen.product(b * rsrp[1], pr));
							strahl.calculateIntensity();
							target.incVec_k(strahl.rv);

							strahl.rv = md;

							target.incPower(power_vor - strahl.I);
							abs_power_of_reflections[anz_reflexe] += power_vor - strahl.I;
						}
						// Falls das getroffene Dreieck ein Power Detektor ist, wird nur die Leistung aufgezeichnet
						else if (target.type == 1) {
							target.incPower(strahl.I);
							write_tracks =true;
						}
					}
					else {
						power_lost += strahl.I;
						break;
					}
				}
				// ausgabe.println();
			}
			if (write_tracks){ // hit detector
				Algorithmen.write_beam_tracks_to_file(beam_track_file, beam_tracks);
			}
		}
		synchronized (rt) {
			rt.notify();
		}
	}
}
