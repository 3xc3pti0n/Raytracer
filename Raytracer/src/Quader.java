import java.util.Hashtable;

/**
 * Objekte dieser Klasse stellen Quader im Raum dar, welche wiederum Quader oder auch Dreiecke enthalten koennen
 */
public class Quader {
	double[] center;
	double[][] base; // die Basisvektoren stehen Spaltenweise
	double ax; // Halbe Ausdehnung des Quaders in X Richtung
	double ay; // Halbe Ausdehnung des Quaders in Y Richtung
	double az; // Halbe Ausdehnung des Quaders in Z Richtung
	double[] cbottom; // Koordinaten des Mittelpunktes der unteren xy Ebene (rechtshaendiges Koordinatensystem)
	double[] ctop; // Koordinaten des Mittelpunktes der oberen xy Ebene
	double[] cleft; // Koordinaten des Mittelpunktes der linken yz Ebene
	double[] cright; // Koordinaten des Mittelpunktes der rechten yz Ebene
	double[] cback; // Koordinaten des Mittelpunktes der hinteren xz Ebene
	double[] cfront; // Koordinaten des Mittelpunktes der vorderen xz Ebene

	Quader q1; // Quader 1 der diesen Quader aufteilt
	Quader q2; // Quader 2 der diesen Quader aufteilt

	Objectlist triangles; // Liste von Dreiecken, welche in diesem Quader enthalten sind
	int number_of_triangles; // Integervariable enthaelt die Anzahl aller in diesem Quader enthaltenden Dreiecke

	int level; // Hierachiestufe, 0 bedeutet die unterste Ebene

	static int count = 0; // Zaehler, wieviele Tests gemacht wurden, fuer Performancetest

	/**
	 * Konstruktor erzeugt Quader aus einem Mittelpunkt und der halben Ausdehnung in x, y und z Richtung
	 */
	public Quader(double[] c, double ax, double ay, double az) {
		center = new double[3];
		for (int i = 0; i < 3; i++) {
			center[i] = c[i];
		}
		this.ax = ax;
		this.ay = ay;
		this.az = az;
		triangles = new Objectlist();
		number_of_triangles = 0;
		cbottom = new double[] { center[0], center[1], center[2] - az };
		ctop = new double[] { center[0], center[1], center[2] + az };
		cleft = new double[] { center[0] - ax, center[1], center[2] };
		cright = new double[] { center[0] + ax, center[1], center[2] };
		cfront = new double[] { center[0], center[1] - ay, center[2] };
		cback = new double[] { center[0], center[1] + ay, center[2] };
	}

	/**
	 * Methode erzeugt die Baumstruktur der aufteilenden Quader mit angegebener Tiefe
	 * 
	 * @depth Tiefe der Hierachie
	 */
	public void generateSublevels(int depth) {
		level = depth;
		if (level > 0) {
			if ((ax >= ay) && (ax >= az)) { // x-Ausdehnung am groessten
				q1 = new Quader(new double[] { center[0] - 0.5 * ax, center[1], center[2] }, ax / 2, ay, az);
				q1.generateSublevels(depth - 1);
				q2 = new Quader(new double[] { center[0] + 0.5 * ax, center[1], center[2] }, ax / 2, ay, az);
				q2.generateSublevels(depth - 1);
			}
			else if ((ay >= ax) && (ay >= az)) { // y-Ausdehnung am groessten
				q1 = new Quader(new double[] { center[0], center[1] - 0.5 * ay, center[2] }, ax, ay / 2, az);
				q1.generateSublevels(depth - 1);
				q2 = new Quader(new double[] { center[0], center[1] + 0.5 * ay, center[2] }, ax, ay / 2, az);
				q2.generateSublevels(depth - 1);
			}
			else { // z-Ausdehnung am groessten
				q1 = new Quader(new double[] { center[0], center[1], center[2] - 0.5 * az }, ax, ay, az / 2);
				q1.generateSublevels(depth - 1);
				q2 = new Quader(new double[] { center[0], center[1], center[2] + 0.5 * az }, ax, ay, az / 2);
				q2.generateSublevels(depth - 1);
			}
		}
	}

	/**
	 * Methode gibt eine Objectlist mit allen an den Blaettern des Baumes vorhandenen Quadern zurueck
	 */
	public Objectlist getQuaderList() {
		final Objectlist list = new Objectlist();
		if (level == 0) {
			list.add(this);
		}
		else {
			list.add(q1.getQuaderList());
			list.add(q2.getQuaderList());
		}
		return list;
	}

	/**
	 * Methode gibt eine Objectlist mit allen an den Blaettern des Baumes vorhandenen Dreiecken zurueck
	 */
	public Objectlist getTriangleList() {
		final Objectlist list = new Objectlist();
		if (level == 0) {
			list.add(triangles);
		}
		else {
			list.add(q1.getTriangleList());
			list.add(q2.getTriangleList());
		}
		return list;
	}

	/**
	 * Methode sortiert ein uebergebendes Dreieck in die entsprechenden Blatt-Quader, welche das Dreieck enthalten
	 */
	public void relateTriangle(Triangle t) {
		if (this.inside(t)) {
			if (level == 0) {
				triangles.add(t);
				number_of_triangles++;
				// System.out.println(t.name+" einsortiert.");
			}
			else {
				q1.relateTriangle(t);
				q2.relateTriangle(t);
				number_of_triangles = q1.number_of_triangles + q2.number_of_triangles;
			}
		}
	}

	/**
	 * Methode sortiert eine Liste von Dreiecken in die entsprechenden Blatt-Quader
	 */
	public void relateTriangle(Triangle[] t) {
		for (int i = 0; i < t.length; i++) {
			this.relateTriangle(t[i]);
		}
	}

	/**
	 * Methode gibt eine Objectlist mit allen Quadern zurueck, die von der Geraden geschnitten werden und Blaetter des Baumes sind
	 */
	public Objectlist getQuaderIntersectionList(double[] s, double[] t) {
		final Objectlist list = new Objectlist();
		final boolean inter = intersect(s, t);
		if (level == 0) {
			if (inter) {
				list.add(this);
			}
		}
		else if (inter) {
			final Objectlist l1 = q1.getQuaderIntersectionList(s, t);
			final Objectlist l2 = q2.getQuaderIntersectionList(s, t);
			list.add(l1);
			list.add(l2);
		}
		return list;
	}

	/**
	 * Methode gibt eine Objectlist mit allen Quadern zurueck, die vom Strahl geschnitten werden und Blaetter des Baumes sind
	 */
	public Objectlist getQuaderIntersectionList(Ray ray) {
		final Objectlist list = new Objectlist();
		final boolean inter = intersect(ray);
		if (level == 0) {
			if (inter) {
				list.add(this);
			}
		}
		else if (inter) {
			final Objectlist l1 = q1.getQuaderIntersectionList(ray);
			final Objectlist l2 = q2.getQuaderIntersectionList(ray);
			list.add(l1);
			list.add(l2);
		}
		return list;
	}

	/**
	 * Methode gibt eine Objectlist mit allen Quadern zurueck, die von der Geraden geschnitten werden und Blaetter des Baumes sind UND AUCH DREIECKE ENTHALTEN
	 */
	public Objectlist getQuaderIntersectionList_NotEmpty(double[] s, double[] t) {
		final Objectlist list = new Objectlist();
		if (number_of_triangles > 0) // Quader brauchen nur getestet werden, wenn sie ueberhaupt Dreiecke enthalten
		{
			final boolean inter = intersect(s, t);
			if (level == 0) {
				if (inter) {
					list.add(this);
				}
			}
			else if (inter) {
				final Objectlist l1 = q1.getQuaderIntersectionList_NotEmpty(s, t);
				final Objectlist l2 = q2.getQuaderIntersectionList_NotEmpty(s, t);
				list.add(l1);
				list.add(l2);
			}
		}
		return list;
	}

	/**
	 * Methode gibt eine Objectlist mit allen Quadern zurueck, die von der Geraden geschnitten werden und Blaetter des Baumes sind UND AUCH DREIECKE ENTHALTEN
	 */
	public Objectlist getQuaderIntersectionList_NotEmpty(Ray ray) {
		final Objectlist list = new Objectlist();
		if (number_of_triangles > 0) // Quader brauchen nur getestet werden, wenn sie ueberhaupt Dreiecke enthalten
		{
			final boolean inter = intersect(ray);
			if (level == 0) {
				if (inter) {
					list.add(this);
				}
			}
			else if (inter) {
				final Objectlist l1 = q1.getQuaderIntersectionList_NotEmpty(ray);
				final Objectlist l2 = q2.getQuaderIntersectionList_NotEmpty(ray);
				list.add(l1);
				list.add(l2);
			}
		}
		return list;
	}

	/**
	 * Methode gibt eine Objectlist mit allen Dreiecken zurueck, die in Quadern liegen, welche von der Geraden geschnitten werden und Blaetter des Baumes sind.
	 */
	public Objectlist getQuaderIntersectionTriangleList(double[] s, double[] t) {
		final Objectlist list = new Objectlist();
		if (number_of_triangles > 0) // Quader brauchen nur getestet werden, wenn sie ueberhaupt Dreiecke enthalten
		{
			final boolean inter = intersect(s, t);
			if (level == 0) {
				if (inter) {
					list.add(triangles);
				}
			}
			else if (inter) {
				final Objectlist l1 = q1.getQuaderIntersectionTriangleList(s, t);
				final Objectlist l2 = q2.getQuaderIntersectionTriangleList(s, t);
				list.add(l1);
				list.add(l2);
			}
		}
		return list;
	}

	/**
	 * Methode gibt eine Objectlist mit allen Dreiecken zurueck, die in Quadern liegen, welche von der Geraden geschnitten werden und Blaetter des Baumes sind.
	 */
	public Objectlist getQuaderIntersectionTriangleList(Ray ray) {
		final Objectlist list = new Objectlist();
		if (number_of_triangles > 0) // Quader brauchen nur getestet werden, wenn sie ueberhaupt Dreiecke enthalten
		{
			final boolean inter = intersect(ray);
			if (level == 0) {
				if (inter) {
					list.add(triangles);
				}
			}
			else if (inter) {
				final Objectlist l1 = q1.getQuaderIntersectionTriangleList(ray);
				final Objectlist l2 = q2.getQuaderIntersectionTriangleList(ray);
				list.add(l1);
				list.add(l2);
			}
		}
		return list;
	}

	/**
	 * Methode gibt eine Objectlist mit allen Dreiecken zurueck, die in Quadern liegen, welche von der Geraden geschnitten werden und Blaetter des Baumes sind.
	 */
	public Objectlist getQuaderIntersectionTriangleList_noduplicates(double[] s, double[] t) {
		final Objectlist list = getQuaderIntersectionTriangleList(s, t);
		final Objectlist list_clean = new Objectlist();
		final int count = 0;
		if (list.N > 0) {
			final Hashtable<String, Integer> hash = new Hashtable<String, Integer>(100);
			// Hashtable hash = new Hashtable(100);
			list.rewind();
			Triangle tri;
			for (int i = 0; i < list.N; i++) {
				tri = (Triangle) list.getActual();
				if (!hash.containsKey(tri.name)) {
					hash.put(tri.name, 1);
					list_clean.add(tri);
				}
				list.next();
			}
		}
		return list_clean;
	}

	/**
	 * Methode gibt eine Objectlist mit allen Dreiecken zurueck, die in Quadern liegen, welche von der Geraden geschnitten werden und Blaetter des Baumes sind.
	 */
	public Objectlist getQuaderIntersectionTriangleList_noduplicates(Ray ray) {
		final Objectlist list = getQuaderIntersectionTriangleList(ray);
		final Objectlist list_clean = new Objectlist();
		final int count = 0;
		if (list.N > 0) {
			final Hashtable<String, Integer> hash = new Hashtable<String, Integer>(100);
			// Hashtable hash = new Hashtable(100);
			list.rewind();
			Triangle tri;
			for (int i = 0; i < list.N; i++) {
				tri = (Triangle) list.getActual();
				if (!hash.containsKey(tri.name)) {
					hash.put(tri.name, 1);
					list_clean.add(tri);
				}
				list.next();
			}
		}
		return list_clean;
	}

	/**
	 * Methode liefert true zurueck, falls der Quader von einem Strahl geschnitten wird. Es wird nur die "vorwaertsrichtung" des Strahls beruecksichtigt.
	 * 
	 * @ray Strahl
	 */
	public boolean intersect(Ray ray) {
		boolean inter = false;
		if (ray.kaustik == false) {
			inter = intersect(ray.ov, ray.rv);
		}
		else {
			double[] sp; // Koordinaten des Schnittpunktes
			double[] xparams;

			// Teste den Schnitt eines Strahls mit der Bottom Ebene
			// Gehe alle Parameter durch, die in den Strahl bis zum Schnittpunkt eingesetzt werden muessen
			xparams = Algorithmen.paramx_intersection_plane_gaussianray(cbottom, new double[] { 0, 0, 1 }, ray);
			for (int i = 0; i < (int) xparams[0]; i++) {
				if ((!inter) && (xparams[i + 1] >= 0)) {
					// Schritt 1: Bestimme Schnittpunkt mit der Ebene
					sp = Algorithmen.point(ray, xparams[i + 1]);// i+1 weil 0tes Arrayelement Anzahl der Loesungen im Array
					// Schritt 2: Pruefe, ob x-Koordinate zwischen cleft und cright und y-Koodinate zwischen cfront und cback
					if ((sp[0] >= cleft[0]) && (sp[0] <= cright[0]) && (sp[1] >= cfront[1]) && (sp[1] <= cback[1])) {
						// Aha, Schnittpunkt vorhanden!
						inter = true;
					}
				}
			}
			if (!inter) {
				// Teste den Schnitt eines Strahls mit der Top Ebene
				// Gehe alle Parameter durch, die in den Strahl bis zum Schnittpunkt eingesetzt werden muessen
				xparams = Algorithmen.paramx_intersection_plane_gaussianray(ctop, new double[] { 0, 0, 1 }, ray);
				for (int i = 0; i < (int) xparams[0]; i++) {
					if ((!inter) && (xparams[i + 1] >= 0)) {
						// Schritt 1: Bestimme Schnittpunkt mit der Ebene
						sp = Algorithmen.point(ray, xparams[i + 1]);// i+1 weil 0tes Arrayelement Anzahl der Loesungen im Array
						// Schritt 2: Pruefe, ob x-Koordinate zwischen cleft und cright und y-Koodinate zwischen cfront und cback
						if ((sp[0] >= cleft[0]) && (sp[0] <= cright[0]) && (sp[1] >= cfront[1]) && (sp[1] <= cback[1])) {
							// Aha, Schnittpunkt vorhanden!
							inter = true;
						}
					}
				}
			}
			if (!inter) {
				// Teste den Schnitt eines Strahls mit der linken Ebene
				// Gehe alle Parameter durch, die in den Strahl bis zum Schnittpunkt eingesetzt werden muessen
				xparams = Algorithmen.paramx_intersection_plane_gaussianray(cleft, new double[] { 1, 0, 0 }, ray);// !
				for (int i = 0; i < (int) xparams[0]; i++) {
					if ((!inter) && (xparams[i + 1] >= 0)) {
						// Schritt 1: Bestimme Schnittpunkt mit der Ebene
						sp = Algorithmen.point(ray, xparams[i + 1]);// i+1 weil 0tes Arrayelement Anzahl der Loesungen im Array
						// Schritt 2: Pruefe, ob Koordinaten innerhalb der Grenzen liegen
						if ((sp[1] >= cfront[1]) && (sp[1] <= cback[1]) && (sp[2] >= cbottom[2]) && (sp[2] <= ctop[2])) {
							// Aha, Schnittpunkt vorhanden!
							inter = true;
						}
					}
				}
			}
			if (!inter) {
				// Teste den Schnitt eines Strahls mit der rechten Ebene
				// Gehe alle Parameter durch, die in den Strahl bis zum Schnittpunkt eingesetzt werden muessen
				xparams = Algorithmen.paramx_intersection_plane_gaussianray(cright, new double[] { 1, 0, 0 }, ray);// !
				for (int i = 0; i < (int) xparams[0]; i++) {
					if ((!inter) && (xparams[i + 1] >= 0)) {
						// Schritt 1: Bestimme Schnittpunkt mit der Ebene
						sp = Algorithmen.point(ray, xparams[i + 1]);// i+1 weil 0tes Arrayelement Anzahl der Loesungen im Array
						// Schritt 2: Pruefe, ob Koordinaten innerhalb der Grenzen liegen
						if ((sp[1] >= cfront[1]) && (sp[1] <= cback[1]) && (sp[2] >= cbottom[2]) && (sp[2] <= ctop[2])) {
							// Aha, Schnittpunkt vorhanden!
							inter = true;
						}
					}
				}
			}
			if (!inter) {
				// Teste den Schnitt eines Strahls mit der vorderen Ebene
				// Gehe alle Parameter durch, die in den Strahl bis zum Schnittpunkt eingesetzt werden muessen
				xparams = Algorithmen.paramx_intersection_plane_gaussianray(cfront, new double[] { 0, 1, 0 }, ray);// !
				for (int i = 0; i < (int) xparams[0]; i++) {
					if ((!inter) && (xparams[i + 1] >= 0)) {
						// Schritt 1: Bestimme Schnittpunkt mit der Ebene
						sp = Algorithmen.point(ray, xparams[i + 1]);// i+1 weil 0tes Arrayelement Anzahl der Loesungen im Array
						// Schritt 2: Pruefe, ob Koordinaten innerhalb der Grenzen liegen
						if ((sp[0] >= cleft[0]) && (sp[0] <= cright[0]) && (sp[2] >= cbottom[2]) && (sp[2] <= ctop[2])) {
							// Aha, Schnittpunkt vorhanden!
							inter = true;
						}
					}
				}
			}
			if (!inter) {
				// Teste den Schnitt eines Strahls mit der hinteren Ebene
				// Gehe alle Parameter durch, die in den Strahl bis zum Schnittpunkt eingesetzt werden muessen
				xparams = Algorithmen.paramx_intersection_plane_gaussianray(cback, new double[] { 0, 1, 0 }, ray);// !
				for (int i = 0; i < (int) xparams[0]; i++) {
					if ((!inter) && (xparams[i + 1] >= 0)) {
						// Schritt 1: Bestimme Schnittpunkt mit der Ebene
						sp = Algorithmen.point(ray, xparams[i + 1]);// i+1 weil 0tes Arrayelement Anzahl der Loesungen im Array
						// Schritt 2: Pruefe, ob Koordinaten innerhalb der Grenzen liegen
						if ((sp[0] >= cleft[0]) && (sp[0] <= cright[0]) && (sp[2] >= cbottom[2]) && (sp[2] <= ctop[2])) {
							// Aha, Schnittpunkt vorhanden!
							inter = true;
						}
					}
				}
			}
		}
		return inter;
	}

	/**
	 * Methode liefert true zurueck, falls der Quader von einer Geraden geschnitten wird. Es wird nur die "vorwaertsrichtung" der Gerade beruecksichtigt.
	 * 
	 * @s Ortsvektor der Geraden
	 * @t Richtungsvektor der Geraden
	 */
	public boolean intersect(double[] s, double[] t) {
		count++;
		// System.out.println("Intersect Test "+count+" ausgefuehrt.");
		boolean inter = false;
		double x;
		final double[] sp = new double[] { 0, 0, 0 }; // Koordinaten des Schnittpunktes
		// Teste den Schnitt einer Geraden s+x*t mit der Bottom Ebene (s,t sind Vektoren)
		// Schritt 1: Ist die Gerade parallel zur xy Ebene? => kein einzelner SchnittPUNKT moeglich
		if (t[2] != 0) {
			// Schritt 2: gibt es einen Schnitt in vorwaertsrichtung?
			// Bestimme den Parameter x
			x = (cbottom[2] - s[2]) / t[2];
			if (x >= 0) {
				// Schritt 3: Bestimme den Schnittpunkt mit der Bottomebene, zunaechst x-Koordinate
				sp[0] = s[0] + x * t[0];
				// Schritt 4a: Pruefe ob der x-Wert zwischen der x-Koordinate der linken Ebene und der rechten Ebene liegt
				if ((sp[0] >= cleft[0]) && (sp[0] <= cright[0])) {
					// Schritt 4b: Pruefe ob der y-Wert zwischen der y-Koordinate der vorderen Ebene und der hinteren Ebene liegt
					sp[1] = s[1] + x * t[1]; // y Koordinate bestimmen
					if ((sp[1] >= cfront[1]) && (sp[1] <= cback[1])) {
						// Aha, Schnittpunkt vorhanden!
						inter = true;
					}
				}
			}
		}
		if (!inter) { // Falls noch kein Schnittpunkt
			// Teste den Schnitt einer Geraden s+x*t mit der Top Ebene (s,t sind Vektoren)
			// Schritt 1: Ist die Gerade parallel zur xy Ebene? => kein einzelner SchnittPUNKT moeglich
			if (t[2] != 0) {
				// Schritt 2: gibt es einen Schnitt in vorwaertsrichtung?
				// Bestimme den Parameter x
				x = (ctop[2] - s[2]) / t[2];
				if (x >= 0) {
					// Schritt 3: Bestimme den Schnittpunkt mit der Topebene, zunaechst x-Koordinate
					sp[0] = s[0] + x * t[0];
					// Schritt 4a: Pruefe ob der x-Wert zwischen der x-Koordinate der linken Ebene und der rechten Ebene liegt
					if ((sp[0] >= cleft[0]) && (sp[0] <= cright[0])) {
						// Schritt 4b: Pruefe ob der y-Wert zwischen der y-Koordinate der vorderen Ebene und der hinteren Ebene liegt
						sp[1] = s[1] + x * t[1]; // y Koordinate bestimmen
						if ((sp[1] >= cfront[1]) && (sp[1] <= cback[1])) {
							// Aha, Schnittpunkt vorhanden!
							inter = true;
						}
					}
				}
			}
		}
		if (!inter) { // Falls noch kein Schnittpunkt
			// Teste den Schnitt einer Geraden s+x*t mit der linken Ebene (s,t sind Vektoren)
			// Schritt 1: Ist die Gerade parallel zur yz Ebene? => kein einzelner SchnittPUNKT moeglich
			if (t[0] != 0) {
				// Schritt 2: gibt es einen Schnitt in vorwaertsrichtung?
				// Bestimme den Parameter x
				x = (cleft[0] - s[0]) / t[0];
				if (x >= 0) {
					// Schritt 3: Bestimme den Schnittpunkt mit der linken ebene, zunaechst y-Koordinate
					sp[1] = s[1] + x * t[1];
					// Schritt 4a: Pruefe ob der y Wert zwischen der y-Koordinate der vorderen Ebene und der hinteren Ebene liegt
					if ((sp[1] >= cfront[1]) && (sp[1] <= cback[1])) {
						// Schritt 4b: Pruefe ob der z Wert zwischen der z-Koordinate der unteren Ebene und der oberen Ebene liegt
						sp[2] = s[2] + x * t[2]; // z Koordinate bestimmen
						if ((sp[2] >= cbottom[2]) && (sp[2] <= ctop[2])) {
							// Aha, Schnittpunkt vorhanden!
							inter = true;
						}
					}
				}
			}
		}
		if (!inter) { // Falls noch kein Schnittpunkt
			// Teste den Schnitt einer Geraden s+x*t mit der rechten Ebene (s,t sind Vektoren)
			// Schritt 1: Ist die Gerade parallel zur yz Ebene? => kein einzelner SchnittPUNKT moeglich
			if (t[0] != 0) {
				// Schritt 2: gibt es einen Schnitt in vorwaertsrichtung?
				// Bestimme den Parameter x
				x = (cright[0] - s[0]) / t[0];
				if (x >= 0) {
					// Schritt 3: Bestimme den Schnittpunkt mit der rechten ebene, zunaechst y-Koordinate
					sp[1] = s[1] + x * t[1];
					// Schritt 4a: Pruefe ob der y Wert zwischen der y-Koordinate der vorderen Ebene und der hinteren Ebene liegt
					if ((sp[1] >= cfront[1]) && (sp[1] <= cback[1])) {
						// Schritt 4b: Pruefe ob der z Wert zwischen der z-Koordinate der unteren Ebene und der oberen Ebene liegt
						sp[2] = s[2] + x * t[2]; // z Koordinate bestimmen
						if ((sp[2] >= cbottom[2]) && (sp[2] <= ctop[2])) {
							// Aha, Schnittpunkt vorhanden!
							inter = true;
						}
					}
				}
			}
		}
		if (!inter) { // Falls noch kein Schnittpunkt
			// Teste den Schnitt einer Geraden s+x*t mit der vorderen Ebene (s,t sind Vektoren)
			// Schritt 1: Ist die Gerade parallel zur xz Ebene? => kein einzelner SchnittPUNKT moeglich
			if (t[1] != 0) {
				// Schritt 2: gibt es einen Schnitt in vorwaertsrichtung?
				// Bestimme den Parameter x
				x = (cfront[1] - s[1]) / t[1];
				if (x >= 0) {
					// Schritt 3: Bestimme den Schnittpunkt mit der vorderen ebene, zunaechst x-Koordinate
					sp[0] = s[0] + x * t[0];
					// Schritt 4a: Pruefe ob der x Wert zwischen der x-Koordinate der linken Ebene und der rechten Ebene liegt
					if ((sp[0] >= cleft[0]) && (sp[0] <= cright[0])) {
						// Schritt 4b: Pruefe ob der z Wert zwischen der z-Koordinate der unteren Ebene und der oberen Ebene liegt
						sp[2] = s[2] + x * t[2]; // z Koordinate bestimmen
						if ((sp[2] >= cbottom[2]) && (sp[2] <= ctop[2])) {
							// Aha, Schnittpunkt vorhanden!
							inter = true;
						}
					}
				}
			}
		}
		if (!inter) { // Falls noch kein Schnittpunkt
			// Teste den Schnitt einer Geraden s+x*t mit der hinteren Ebene (s,t sind Vektoren)
			// Schritt 1: Ist die Gerade parallel zur xz Ebene? => kein einzelner SchnittPUNKT moeglich
			if (t[1] != 0) {
				// Schritt 2: gibt es einen Schnitt in vorwaertsrichtung?
				// Bestimme den Parameter x
				x = (cback[1] - s[1]) / t[1];
				if (x >= 0) {
					// Schritt 3: Bestimme den Schnittpunkt mit der hinteren ebene, zunaechst x-Koordinate
					sp[0] = s[0] + x * t[0];
					// Schritt 4a: Pruefe ob der x Wert zwischen der x-Koordinate der linken Ebene und der rechten Ebene liegt
					if ((sp[0] >= cleft[0]) && (sp[0] <= cright[0])) {
						// Schritt 4b: Pruefe ob der z Wert zwischen der z-Koordinate der unteren Ebene und der oberen Ebene liegt
						sp[2] = s[2] + x * t[2]; // z Koordinate bestimmen
						if ((sp[2] >= cbottom[2]) && (sp[2] <= ctop[2])) {
							// Aha, Schnittpunkt vorhanden!
							inter = true;
						}
					}
				}
			}
		}
		return inter;
	}

	/**
	 * Methode liefert true zurueck, falls der Quader von einer Geraden geschnitten wird. Es wird nur die "vorwaertsrichtung" der Gerade beruecksichtigt. Weiterhin wird ein Array mit den Schnittpunkten erzeugt. Ist zum Debuggen gedacht.
	 * 
	 * @s Ortsvektor der Geraden
	 * @t Richtungsvektor der Geraden
	 */
	public boolean intersect(double[] s, double[] t, double[][] schnittpunkte) {
		boolean inter = false;
		double x;
		final double[] sp = new double[] { 0, 0, 0 }; // Koordinaten des Schnittpunktes
		// Teste den Schnitt einer Geraden s+x*t mit der Bottom Ebene (s,t sind Vektoren)
		// Schritt 1: Ist die Gerade parallel zur xy Ebene? => kein einzelner SchnittPUNKT moeglich
		if (t[2] != 0) {
			// Schritt 2: gibt es einen Schnitt in vorwaertsrichtung?
			// Bestimme den Parameter x
			x = (cbottom[2] - s[2]) / t[2];
			if (x >= 0) {
				// Schritt 3: Bestimme den Schnittpunkt mit der Bottomebene, zunaechst x-Koordinate
				sp[0] = s[0] + x * t[0];
				// Schritt 4a: Pruefe ob der x-Wert zwischen der x-Koordinate der linken Ebene und der rechten Ebene liegt
				if ((sp[0] >= cleft[0]) && (sp[0] <= cright[0])) {
					// Schritt 4b: Pruefe ob der y-Wert zwischen der y-Koordinate der vorderen Ebene und der hinteren Ebene liegt
					sp[1] = s[1] + x * t[1]; // y Koordinate bestimmen
					if ((sp[1] >= cfront[1]) && (sp[1] <= cback[1])) {
						// Aha, Schnittpunkt vorhanden!
						inter = true;
						sp[2] = cbottom[2];
						schnittpunkte[0] = sp.clone();
						System.out.println("Schnittpunkt mit unterer Ebene\n" + Algorithmen.ArrayToString(sp));
					}
				}
			}
		}
		// if (!inter){ //Falls noch kein Schnittpunkt
		// Teste den Schnitt einer Geraden s+x*t mit der Top Ebene (s,t sind Vektoren)
		// Schritt 1: Ist die Gerade parallel zur xy Ebene? => kein einzelner SchnittPUNKT moeglich
		if (t[2] != 0) {
			// Schritt 2: gibt es einen Schnitt in vorwaertsrichtung?
			// Bestimme den Parameter x
			x = (ctop[2] - s[2]) / t[2];
			if (x >= 0) {
				// Schritt 3: Bestimme den Schnittpunkt mit der Topebene, zunaechst x-Koordinate
				sp[0] = s[0] + x * t[0];
				// Schritt 4a: Pruefe ob der x-Wert zwischen der x-Koordinate der linken Ebene und der rechten Ebene liegt
				if ((sp[0] >= cleft[0]) && (sp[0] <= cright[0])) {
					// Schritt 4b: Pruefe ob der y-Wert zwischen der y-Koordinate der vorderen Ebene und der hinteren Ebene liegt
					sp[1] = s[1] + x * t[1]; // y Koordinate bestimmen
					if ((sp[1] >= cfront[1]) && (sp[1] <= cback[1])) {
						// Aha, Schnittpunkt vorhanden!
						inter = true;
						sp[2] = ctop[2];
						schnittpunkte[1] = sp.clone();
						System.out.println("Schnittpunkt mit oberer Ebene\n" + Algorithmen.ArrayToString(sp));
					}
				}
			}
		}
		// if (!inter){ //Falls noch kein Schnittpunkt
		// Teste den Schnitt einer Geraden s+x*t mit der linken Ebene (s,t sind Vektoren)
		// Schritt 1: Ist die Gerade parallel zur yz Ebene? => kein einzelner SchnittPUNKT moeglich
		if (t[0] != 0) {
			// Schritt 2: gibt es einen Schnitt in vorwaertsrichtung?
			// Bestimme den Parameter x
			x = (cleft[0] - s[0]) / t[0];
			if (x >= 0) {
				// Schritt 3: Bestimme den Schnittpunkt mit der linken ebene, zunaechst y-Koordinate
				sp[1] = s[1] + x * t[1];
				// Schritt 4a: Pruefe ob der y Wert zwischen der y-Koordinate der vorderen Ebene und der hinteren Ebene liegt
				if ((sp[1] >= cfront[1]) && (sp[1] <= cback[1])) {
					// Schritt 4b: Pruefe ob der z Wert zwischen der z-Koordinate der unteren Ebene und der oberen Ebene liegt
					sp[2] = s[2] + x * t[2]; // z Koordinate bestimmen
					if ((sp[2] >= cbottom[2]) && (sp[2] <= ctop[2])) {
						// Aha, Schnittpunkt vorhanden!
						inter = true;
						sp[0] = cleft[0];
						schnittpunkte[2] = sp.clone();
						System.out.println("Schnittpunkt mit linker Ebene\n" + Algorithmen.ArrayToString(sp));
					}
				}
			}
		}
		// if (!inter){ //Falls noch kein Schnittpunkt
		// Teste den Schnitt einer Geraden s+x*t mit der rechten Ebene (s,t sind Vektoren)
		// Schritt 1: Ist die Gerade parallel zur yz Ebene? => kein einzelner SchnittPUNKT moeglich
		if (t[0] != 0) {
			// Schritt 2: gibt es einen Schnitt in vorwaertsrichtung?
			// Bestimme den Parameter x
			x = (cright[0] - s[0]) / t[0];
			if (x >= 0) {
				// Schritt 3: Bestimme den Schnittpunkt mit der rechten ebene, zunaechst y-Koordinate
				sp[1] = s[1] + x * t[1];
				// Schritt 4a: Pruefe ob der y Wert zwischen der y-Koordinate der vorderen Ebene und der hinteren Ebene liegt
				if ((sp[1] >= cfront[1]) && (sp[1] <= cback[1])) {
					// Schritt 4b: Pruefe ob der z Wert zwischen der z-Koordinate der unteren Ebene und der oberen Ebene liegt
					sp[2] = s[2] + x * t[2]; // z Koordinate bestimmen
					if ((sp[2] >= cbottom[2]) && (sp[2] <= ctop[2])) {
						// Aha, Schnittpunkt vorhanden!
						inter = true;
						sp[0] = cright[0];
						schnittpunkte[3] = sp.clone();
						System.out.println("Schnittpunkt mit rechter Ebene\n" + Algorithmen.ArrayToString(sp));
					}
				}
			}
		}
		// if (!inter){ //Falls noch kein Schnittpunkt
		// Teste den Schnitt einer Geraden s+x*t mit der vorderen Ebene (s,t sind Vektoren)
		// Schritt 1: Ist die Gerade parallel zur xz Ebene? => kein einzelner SchnittPUNKT moeglich
		if (t[1] != 0) {
			// Schritt 2: gibt es einen Schnitt in vorwaertsrichtung?
			// Bestimme den Parameter x
			x = (cfront[1] - s[1]) / t[1];
			if (x >= 0) {
				// Schritt 3: Bestimme den Schnittpunkt mit der vorderen ebene, zunaechst x-Koordinate
				sp[0] = s[0] + x * t[0];
				// Schritt 4a: Pruefe ob der x Wert zwischen der x-Koordinate der linken Ebene und der rechten Ebene liegt
				if ((sp[0] >= cleft[0]) && (sp[0] <= cright[0])) {
					// Schritt 4b: Pruefe ob der z Wert zwischen der z-Koordinate der unteren Ebene und der oberen Ebene liegt
					sp[2] = s[2] + x * t[2]; // z Koordinate bestimmen
					if ((sp[2] >= cbottom[2]) && (sp[2] <= ctop[2])) {
						// Aha, Schnittpunkt vorhanden!
						inter = true;
						sp[1] = cfront[1];
						schnittpunkte[4] = sp.clone();
						System.out.println("Schnittpunkt mit vorderer Ebene\n" + Algorithmen.ArrayToString(sp));
					}
				}
			}
		}
		// if (!inter){ //Falls noch kein Schnittpunkt
		// Teste den Schnitt einer Geraden s+x*t mit der hinteren Ebene (s,t sind Vektoren)
		// Schritt 1: Ist die Gerade parallel zur xz Ebene? => kein einzelner SchnittPUNKT moeglich
		if (t[1] != 0) {
			// Schritt 2: gibt es einen Schnitt in vorwaertsrichtung?
			// Bestimme den Parameter x
			x = (cback[1] - s[1]) / t[1];
			if (x >= 0) {
				// Schritt 3: Bestimme den Schnittpunkt mit der hinteren ebene, zunaechst x-Koordinate
				sp[0] = s[0] + x * t[0];
				// Schritt 4a: Pruefe ob der x Wert zwischen der x-Koordinate der linken Ebene und der rechten Ebene liegt
				if ((sp[0] >= cleft[0]) && (sp[0] <= cright[0])) {
					// Schritt 4b: Pruefe ob der z Wert zwischen der z-Koordinate der unteren Ebene und der oberen Ebene liegt
					sp[2] = s[2] + x * t[2]; // z Koordinate bestimmen
					if ((sp[2] >= cbottom[2]) && (sp[2] <= ctop[2])) {
						// Aha, Schnittpunkt vorhanden!
						inter = true;
						sp[1] = cback[1];
						schnittpunkte[5] = sp.clone();
						System.out.println("Schnittpunkt mit hinterer Ebene\n" + Algorithmen.ArrayToString(sp));
					}
				}
			}
		}
		return inter;
	}

	/**
	 * Methode liefert true zurueck, falls der Quader von einer Linie geschnitten wird. Der Unterschied zu der Methode intersect besteht darin, dass nur dann true zurueckgegeben wird, wenn der Geradenparameter x der geraden s+x*t zwischen 0 und 1 (einschliesslich) liegt. Das bedeutet, der Schnittpunkt befindet sich zwischen dem Punkt s und dem Punkt s+t.
	 * 
	 * @s Ortsvektor der Geraden
	 * @t Richtungsvektor der Geraden
	 */
	public boolean intersectLine(double[] s, double[] t) {
		boolean inter = false;
		double x;
		final double[] sp = new double[] { 0, 0, 0 }; // Koordinaten des Schnittpunktes
		// Teste den Schnitt einer Geraden s+x*t mit der Bottom Ebene (s,t sind Vektoren)
		// Schritt 1: Ist die Gerade parallel zur xy Ebene? => kein einzelner SchnittPUNKT moeglich
		if (t[2] != 0) {
			// Schritt 2: gibt es einen Schnitt in vorwaertsrichtung?
			// Bestimme den Parameter x
			x = (cbottom[2] - s[2]) / t[2];
			if ((x >= 0) && (x <= 1)) {
				// Schritt 3: Bestimme den Schnittpunkt mit der Bottomebene, zunaechst x-Koordinate
				sp[0] = s[0] + x * t[0];
				// Schritt 4a: Pruefe ob der x-Wert zwischen der x-Koordinate der linken Ebene und der rechten Ebene liegt
				if ((sp[0] >= cleft[0]) && (sp[0] <= cright[0])) {
					// Schritt 4b: Pruefe ob der y-Wert zwischen der y-Koordinate der vorderen Ebene und der hinteren Ebene liegt
					sp[1] = s[1] + x * t[1]; // y Koordinate bestimmen
					if ((sp[1] >= cfront[1]) && (sp[1] <= cback[1])) {
						// Aha, Schnittpunkt vorhanden!
						inter = true;
					}
				}
			}
		}
		if (!inter) { // Falls noch kein Schnittpunkt
			// Teste den Schnitt einer Geraden s+x*t mit der Top Ebene (s,t sind Vektoren)
			// Schritt 1: Ist die Gerade parallel zur xy Ebene? => kein einzelner SchnittPUNKT moeglich
			if (t[2] != 0) {
				// Schritt 2: gibt es einen Schnitt in vorwaertsrichtung?
				// Bestimme den Parameter x
				x = (ctop[2] - s[2]) / t[2];
				if ((x >= 0) && (x <= 1)) {
					// Schritt 3: Bestimme den Schnittpunkt mit der Topebene, zunaechst x-Koordinate
					sp[0] = s[0] + x * t[0];
					// Schritt 4a: Pruefe ob der x-Wert zwischen der x-Koordinate der linken Ebene und der rechten Ebene liegt
					if ((sp[0] >= cleft[0]) && (sp[0] <= cright[0])) {
						// Schritt 4b: Pruefe ob der y-Wert zwischen der y-Koordinate der vorderen Ebene und der hinteren Ebene liegt
						sp[1] = s[1] + x * t[1]; // y Koordinate bestimmen
						if ((sp[1] >= cfront[1]) && (sp[1] <= cback[1])) {
							// Aha, Schnittpunkt vorhanden!
							inter = true;
						}
					}
				}
			}
		}
		if (!inter) { // Falls noch kein Schnittpunkt
			// Teste den Schnitt einer Geraden s+x*t mit der linken Ebene (s,t sind Vektoren)
			// Schritt 1: Ist die Gerade parallel zur yz Ebene? => kein einzelner SchnittPUNKT moeglich
			if (t[0] != 0) {
				// Schritt 2: gibt es einen Schnitt in vorwaertsrichtung?
				// Bestimme den Parameter x
				x = (cleft[0] - s[0]) / t[0];
				if ((x >= 0) && (x <= 1)) {
					// Schritt 3: Bestimme den Schnittpunkt mit der linken ebene, zunaechst y-Koordinate
					sp[1] = s[1] + x * t[1];
					// Schritt 4a: Pruefe ob der y Wert zwischen der y-Koordinate der vorderen Ebene und der hinteren Ebene liegt
					if ((sp[1] >= cfront[1]) && (sp[1] <= cback[1])) {
						// Schritt 4b: Pruefe ob der z Wert zwischen der z-Koordinate der unteren Ebene und der oberen Ebene liegt
						sp[2] = s[2] + x * t[2]; // z Koordinate bestimmen
						if ((sp[2] >= cbottom[2]) && (sp[2] <= ctop[2])) {
							// Aha, Schnittpunkt vorhanden!
							inter = true;
						}
					}
				}
			}
		}
		if (!inter) { // Falls noch kein Schnittpunkt
			// Teste den Schnitt einer Geraden s+x*t mit der rechten Ebene (s,t sind Vektoren)
			// Schritt 1: Ist die Gerade parallel zur yz Ebene? => kein einzelner SchnittPUNKT moeglich
			if (t[0] != 0) {
				// Schritt 2: gibt es einen Schnitt in vorwaertsrichtung?
				// Bestimme den Parameter x
				x = (cright[0] - s[0]) / t[0];
				if ((x >= 0) && (x <= 1)) {
					// Schritt 3: Bestimme den Schnittpunkt mit der rechten ebene, zunaechst y-Koordinate
					sp[1] = s[1] + x * t[1];
					// Schritt 4a: Pruefe ob der y Wert zwischen der y-Koordinate der vorderen Ebene und der hinteren Ebene liegt
					if ((sp[1] >= cfront[1]) && (sp[1] <= cback[1])) {
						// Schritt 4b: Pruefe ob der z Wert zwischen der z-Koordinate der unteren Ebene und der oberen Ebene liegt
						sp[2] = s[2] + x * t[2]; // z Koordinate bestimmen
						if ((sp[2] >= cbottom[2]) && (sp[2] <= ctop[2])) {
							// Aha, Schnittpunkt vorhanden!
							inter = true;
						}
					}
				}
			}
		}
		if (!inter) { // Falls noch kein Schnittpunkt
			// Teste den Schnitt einer Geraden s+x*t mit der vorderen Ebene (s,t sind Vektoren)
			// Schritt 1: Ist die Gerade parallel zur xz Ebene? => kein einzelner SchnittPUNKT moeglich
			if (t[1] != 0) {
				// Schritt 2: gibt es einen Schnitt in vorwaertsrichtung?
				// Bestimme den Parameter x
				x = (cfront[1] - s[1]) / t[1];
				if ((x >= 0) && (x <= 1)) {
					// Schritt 3: Bestimme den Schnittpunkt mit der vorderen ebene, zunaechst x-Koordinate
					sp[0] = s[0] + x * t[0];
					// Schritt 4a: Pruefe ob der x Wert zwischen der x-Koordinate der linken Ebene und der rechten Ebene liegt
					if ((sp[0] >= cleft[0]) && (sp[0] <= cright[0])) {
						// Schritt 4b: Pruefe ob der z Wert zwischen der z-Koordinate der unteren Ebene und der oberen Ebene liegt
						sp[2] = s[2] + x * t[2]; // z Koordinate bestimmen
						if ((sp[2] >= cbottom[2]) && (sp[2] <= ctop[2])) {
							// Aha, Schnittpunkt vorhanden!
							inter = true;
						}
					}
				}
			}
		}
		if (!inter) { // Falls noch kein Schnittpunkt
			// Teste den Schnitt einer Geraden s+x*t mit der hinteren Ebene (s,t sind Vektoren)
			// Schritt 1: Ist die Gerade parallel zur xz Ebene? => kein einzelner SchnittPUNKT moeglich
			if (t[1] != 0) {
				// Schritt 2: gibt es einen Schnitt in vorwaertsrichtung?
				// Bestimme den Parameter x
				x = (cback[1] - s[1]) / t[1];
				if ((x >= 0) && (x <= 1)) {
					// Schritt 3: Bestimme den Schnittpunkt mit der hinteren ebene, zunaechst x-Koordinate
					sp[0] = s[0] + x * t[0];
					// Schritt 4a: Pruefe ob der x Wert zwischen der x-Koordinate der linken Ebene und der rechten Ebene liegt
					if ((sp[0] >= cleft[0]) && (sp[0] <= cright[0])) {
						// Schritt 4b: Pruefe ob der z Wert zwischen der z-Koordinate der unteren Ebene und der oberen Ebene liegt
						sp[2] = s[2] + x * t[2]; // z Koordinate bestimmen
						if ((sp[2] >= cbottom[2]) && (sp[2] <= ctop[2])) {
							// Aha, Schnittpunkt vorhanden!
							inter = true;
						}
					}
				}
			}
		}
		return inter;
	}

	/**
	 * Damit kann man klaeren, ob ein Punkt innerhalb oder ausserhalb dieses Quaders liegt.
	 */
	public boolean inside(double[] p) {
		boolean in = false;
		// zunaechst den Mittelpunktsvektor abziehen, also den Quader in den Ursprung verschieben und den Punkt mitschieben
		final double[] p_ = new double[p.length];
		for (int i = 0; i < 3; i++) {
			p_[i] = p[i] - center[i];
		}
		// Jetzt pruefen, ob einzelnen Koordinaten zwischen -a(xyz) und +a(xyz) (einschliesslich) liegen
		if (((p_[0] <= ax) && (p_[0] >= -ax)) && ((p_[1] <= ay) && (p_[1] >= -ay)) && ((p_[2] <= az) && (p_[2] >= -az))) {
			in = true;
		}
		return in;
	}

	/**
	 * Diese Methode gibt true zurueck, wenn das uebergebene Dreieck (repraesentiert durch die drei Punkte) innerhalb des Quaders liegt. Am einfachsten ist es, wenn bereits ein Dreieckspunkt innerhalb des Quaders liegt, dann ist es klar. Ansonsten ist es noch moeglich, dass eine Dreieicksseite innerhalb des Quaders liegt. Dann muss diese Seite mindestens eine der Quaderflaechen schneiden.
	 */
	public boolean inside(Triangle tri) {
		boolean in = false;
		final double[] p1 = tri.p1;
		final double[] p2 = tri.p2;
		final double[] p3 = tri.p3;
		double[] s = new double[3];
		final double[] t = new double[3];
		// Teste zunaechst die drei Punkte, ob sie bereits innerhalb liegen
		if (inside(p1) == true) {
			in = true;
		}
		else if (inside(p2) == true) {
			in = true;
		}
		else if (inside(p3) == true) {
			in = true;
		}
		else // falls das nicht der Fall ist, teste ob eine der Dreiecksseiten eine Quaderflaeche schneidet
		{
			if (!in) {
				s = p1;
				t[0] = p2[0] - p1[0];
				t[1] = p2[1] - p1[1];
				t[2] = p2[2] - p1[2];
				in = intersectLine(s, t);
			}
			if (!in) {
				s = p1;
				t[0] = p3[0] - p1[0];
				t[1] = p3[1] - p1[1];
				t[2] = p3[2] - p1[2];
				in = intersectLine(s, t);
			}
			if (!in) {
				s = p2;
				t[0] = p3[0] - p2[0];
				t[1] = p3[1] - p2[1];
				t[2] = p3[2] - p2[2];
				in = intersectLine(s, t);
			}
			if (!in) // falls das nicht der Fall ist, teste ob eine der Quaderseiten die Dreiecksflaeche schneidet
			{
				// System.out.println("Start neuer Test");
				// Erhalte alle 12 Quaderkanten, jeweils 2 Punkte bilden die jeweilige Kante
				final double[][][] edges = this.getEdges();
				// Bilde aus den Punktepaaren Geraden
				int j = 0;
				while ((j < 12) && (!in)) {
					final double[] rv = new double[] { edges[j][0][1] - edges[j][0][0], edges[j][1][1] - edges[j][1][0], edges[j][2][1] - edges[j][2][0] };
					final double[] ov = new double[] { edges[j][0][0], edges[j][1][0], edges[j][2][0] };
					final double tnv = Algorithmen.product(rv, tri.nv);
					if (tnv == 0) {
						j++;
						continue;
					}
					final double x = Algorithmen.paramx_intersection_plane_line(tri.p1, tri.nv, ov, Algorithmen.product(rv, tri.nv));
					if ((x >= 0) && (x <= 1)) {
						// Schnittpunkt der Geraden mit der Dreiecksebene bestimmen
						final double[] p = Algorithmen.point(ov, rv, x);
						// Vom Schnittpunkt die Basis des Dreiecks abziehen
						final double[] p_t = new double[3];
						for (int n = 0; n < 3; n++) {
							p_t[n] = p[n] - tri.p1[n];
						}
						// Nur wenn der Schnittpunkt innerhalb der Dreiecksebene liegt, Test 1
						final double[] inbase = Algorithmen.product(tri.transform, p_t);
						// System.out.println("Entwicklungskoeffizienten\n"+inbase[0]+" "+inbase[1]+" "+inbase[2]+"\n");
						if ((inbase[0] >= 0) && (inbase[1] >= 0)) {
							// Nur wenn der Schnittpunkt innerhalb der Dreiecksebene liegt, Test 2
							final double inout = inbase[0] + inbase[1];
							// System.out.println("inout "+inout);
							if (inout <= 1) {
								in = true;
							}
						}
					}
					j++;
					// System.out.println("Ergebnis neuer Test: "+in);
				}
			}
		}
		return in;
	}

	@Override
	public String toString() {
		String s = new String();
		s += "Quader\n";
		s += "Mittelpunkt " + center[0] + "," + center[1] + "," + center[2] + "\n";
		s += "ax " + ax + " ay " + ay + " az " + az + "\n";
		return s;
	}

	/**
	 * Methode gibt ein Array mit allen 8 Eckpunkten des Quaders zurueck In jeder Arrayspalte steht ein Eckpunkt
	 */
	public double[][] getCorners() {
		final double[][] c = new double[3][8];
		double[][] base = new double[3][3];
		base = new double[][] { { ax, 0, 0 }, { 0, ay, 0 }, { 0, 0, az } };
		for (int i = 0; i < 3; i++) {
			c[i][0] = center[i] + base[i][0] + base[i][1] + base[i][2];
		}
		for (int i = 0; i < 3; i++) {
			c[i][1] = center[i] - base[i][0] + base[i][1] + base[i][2];
		}
		for (int i = 0; i < 3; i++) {
			c[i][2] = center[i] - base[i][0] - base[i][1] + base[i][2];
		}
		for (int i = 0; i < 3; i++) {
			c[i][3] = center[i] + base[i][0] - base[i][1] + base[i][2];
		}
		for (int i = 0; i < 3; i++) {
			c[i][4] = center[i] + base[i][0] + base[i][1] - base[i][2];
		}
		for (int i = 0; i < 3; i++) {
			c[i][5] = center[i] - base[i][0] + base[i][1] - base[i][2];
		}
		for (int i = 0; i < 3; i++) {
			c[i][6] = center[i] - base[i][0] - base[i][1] - base[i][2];
		}
		for (int i = 0; i < 3; i++) {
			c[i][7] = center[i] + base[i][0] - base[i][1] - base[i][2];
		}
		return c;
	}

	/**
	 * Methode gibt ein Array mit allen 12 Kanten des Quaders zurueck
	 */
	public double[][][] getEdges() {
		final double[][] c = Algorithmen.transpose(getCorners());
		final double[][][] ed = new double[12][3][2];
		for (int i = 0; i < 3; i++) {
			ed[0][i][0] = c[0][i];
			ed[0][i][1] = c[1][i];
		}
		for (int i = 0; i < 3; i++) {
			ed[1][i][0] = c[0][i];
			ed[1][i][1] = c[3][i];
		}
		for (int i = 0; i < 3; i++) {
			ed[2][i][0] = c[0][i];
			ed[2][i][1] = c[4][i];
		}
		for (int i = 0; i < 3; i++) {
			ed[3][i][0] = c[2][i];
			ed[3][i][1] = c[1][i];
		}
		for (int i = 0; i < 3; i++) {
			ed[4][i][0] = c[2][i];
			ed[4][i][1] = c[3][i];
		}
		for (int i = 0; i < 3; i++) {
			ed[5][i][0] = c[2][i];
			ed[5][i][1] = c[6][i];
		}
		for (int i = 0; i < 3; i++) {
			ed[6][i][0] = c[4][i];
			ed[6][i][1] = c[5][i];
		}
		for (int i = 0; i < 3; i++) {
			ed[7][i][0] = c[4][i];
			ed[7][i][1] = c[7][i];
		}
		for (int i = 0; i < 3; i++) {
			ed[8][i][0] = c[6][i];
			ed[8][i][1] = c[5][i];
		}
		for (int i = 0; i < 3; i++) {
			ed[9][i][0] = c[6][i];
			ed[9][i][1] = c[7][i];
		}
		for (int i = 0; i < 3; i++) {
			ed[10][i][0] = c[3][i];
			ed[10][i][1] = c[7][i];
		}
		for (int i = 0; i < 3; i++) {
			ed[11][i][0] = c[1][i];
			ed[11][i][1] = c[5][i];
		}
		return ed;
	}

}
