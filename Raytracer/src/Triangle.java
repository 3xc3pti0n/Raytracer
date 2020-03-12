import java.util.ArrayList;

/**
 * Objekte dieser Klasse stellen Dreiecke im Raum dar, welche sich in der Regel innerhalb von Quadern befinden
 */
public class Triangle {
	double[] p1; // Eckpunkte des Dreiecks
	double[] p2; // Eckpunkte des Dreiecks
	double[] p3; // Eckpunkte des Dreiecks
	double[] v1; // Vektor zeigt von Punkt 1 zu Punkt 2
	double[] v2; // Vektor zeigt von Punkt 1 zu Punkt 3
	double[] nv; // Normalenvektor steht senkrecht auf v1 und v2, wird durch v1 x v2 berechnet
	double[] nv_p1; // "Normalenvektor" im Punkt P1, wird bei Interpolation des NV benötigt
	double[] nv_p2; // "Normalenvektor" im Punkt P2, wird bei Interpolation des NV benötigt
	double[] nv_p3; // "Normalenvektor" im Punkt P3, wird bei Interpolation des NV benötigt
	double[] translationvector; // Dieser Vektor ist nur wichtig, wenn das Dreieck zur Umsetzung von periodischen Randbedigungen verwendet werden soll.
								// Die Idee ist: Wenn der Dreieckstyp "2" ist, so beinhaltet dieser Vektor die notwendige Verschiebung zum gegenüberliegenden Dreieck auf der anderen Seite, damit der Strahl wieder auf der "anderen Seite" hineinkommt.
	int type; // Gibt an, wie sich das Dreieck bei Bestrahlung verhält.
				// 0: Das Dreieck absorbiert und reflektiert, Anteilig je nach Brechungsindex und Einfallswinkel
				// 1: Das Dreieck speichert die einfallende Leistung ab, läßt den Strahl aber einfach durch. Es wird nicht zwischen unterschieden, ob die Leistung "von oben" oder "von unten" auf das Dreieck trifft
				// 2: Das Dreieck absorbiert nichts. Es passiert nichts. Aber: Der "translationvector" dieses Dreiecks kann abgefragt werden, welcher dann zum Ortsvektor des Strahls addiert werden kann, so dass man periodische Randbedigungen realisieren kann
	double[][] transform; // Transformationsmatrix, welche einen Punkt in das Koordinatensystem transformiert, der von v1, v2 und nv aufgespannt wird
	volatile double P; // Einfallende Leistung
	double [] vec_k = new double [3]; // Einfallende k-Vektoren der Strahlen
	volatile double Temperature;
	volatile double delta_P;
	boolean interpolation; // Boolean gibt an, ob der Interpolationsmodus des Normalenvektors aktiviert ist, oder nicht

	Material material; // Material aus dem das Dreieck besteht
	String name; // Bezeichnung des Dreiecks, muss eindeutig sein, d.h. es darf keine 2 Dreiecke mit gleichem Namen geben
	double[] facet;
	double Power_OF;
	double Kruemmungsradius;
	ArrayList<Integer> neighbourTriangle_P1 = new ArrayList<Integer>();
	ArrayList<Integer> neighbourPoints_P1 = new ArrayList<Integer>();
	ArrayList<Integer> neighbourTriangle_P2 = new ArrayList<Integer>();
	ArrayList<Integer> neighbourPoints_P2 = new ArrayList<Integer>();
	ArrayList<Integer> neighbourTriangle_P3 = new ArrayList<Integer>();
	ArrayList<Integer> neighbourPoints_P3 = new ArrayList<Integer>();
	double[] sPunkt;
	double surfaceTension;
	double OFPressure;
	double DistanceInterface;
	double DeltaPressure;
	double RecoilPressure;
	double Verschiebefaktor;
	double Verschiebung;

	// for smoother

	public double[] Schwerpunkt() {
		sPunkt = new double[3];
		final double teiler = 1.0 / 3.0;
		sPunkt[0] = teiler * (this.p1[0] + this.p2[0] + this.p3[0]);
		sPunkt[1] = teiler * (this.p1[1] + this.p2[1] + this.p3[1]);
		sPunkt[2] = teiler * (this.p1[2] + this.p2[2] + this.p3[2]);
		return sPunkt;
	}

	public double[] getNormalizedFacetNormal() {
		final double corrNorm[] = new double[3];

		nv = Algorithmen.xproduct(v1, v2);

		final double normfaktor = Math.sqrt((this.nv[0] * this.nv[0]) + (this.nv[1] * this.nv[1]) + (this.nv[2] * this.nv[2]));
		// tri.facet[0]=nv[0]/normfaktor;
		// tri.facet[1]=nv[1]/normfaktor;
		// tri.facet[2]=nv[2]/normfaktor;
		corrNorm[0] = nv[0] / normfaktor;
		corrNorm[1] = nv[1] / normfaktor;
		corrNorm[2] = nv[2] / normfaktor;

		// System.out.println(Math.sqrt(corrNorm[0]*corrNorm[0]+corrNorm[1]*corrNorm[1]+corrNorm[2]*corrNorm[2]));

		return corrNorm;

	}

	public int[] vertex_indices = new int[3];

	public void set_vertex_indices(short k, int n) {
		vertex_indices[k] = n;
	}

	//

	public void setNeighbours_P1(int triangle, int point) {
		neighbourTriangle_P1.add(triangle);
		neighbourPoints_P1.add(point);
	}

	public ArrayList<Integer> getNeighbourTriangles_P1() {
		return neighbourTriangle_P1;
	}

	public ArrayList<Integer> getNeighbourPoints_P1() {
		return neighbourPoints_P1;
	}

	public void setNeighbours_P2(int triangle, int point) {
		neighbourTriangle_P2.add(triangle);
		neighbourPoints_P2.add(point);
	}

	public ArrayList<Integer> getNeighbourTriangles_P2() {
		return neighbourTriangle_P2;
	}

	public ArrayList<Integer> getNeighbourPoints_P2() {
		return neighbourPoints_P2;
	}

	public void setNeighbours_P3(int triangle, int point) {
		neighbourTriangle_P3.add(triangle);
		neighbourPoints_P3.add(point);
	}

	public ArrayList<Integer> getNeighbourTriangles_P3() {
		return neighbourTriangle_P3;
	}

	public ArrayList<Integer> getNeighbourPoints_P3() {
		return neighbourPoints_P3;
	}

	public void recalculate_triangle_vector() {

		for (int i = 0; i < 3; i++) {
			this.v1[i] = this.p2[i] - this.p1[i];
			this.v2[i] = this.p3[i] - this.p1[i];

		}
		this.nv = Algorithmen.xproduct(this.v1, this.v2);

	}

	public Triangle(double[] a, double[] b, double[] c, Material mat, String s) {
		p1 = a;
		p2 = b;
		p3 = c;
		v1 = new double[3];
		v2 = new double[3];
		for (int i = 0; i < 3; i++) {
			v1[i] = p2[i] - p1[i];
			v2[i] = p3[i] - p1[i];
		}
		nv = Algorithmen.xproduct(v1, v2);

		this.nv = Algorithmen.getNormTox(this.nv, 1);

		// nv = Algorithmen.product(-1,nv);
		final double[][] m = getBase();
		// transform = Algorithmen.transpose(Algorithmen.inversion(m));
		transform = Algorithmen.inversion(m);
		P = 0.0;
		material = mat;
		name = s;
		interpolation = false;
		type = 0;
	}

	public Triangle(double[] a, double[] b, double[] c, double[] nv_p1, double[] nv_p2, double[] nv_p3, int type, Material mat, String s) {
		this(a, b, c, mat, s);
		this.nv_p1 = nv_p1;
		this.nv_p2 = nv_p2;
		this.nv_p3 = nv_p3;
		interpolation = true;
		this.type = type;
	}

	public Triangle(double[] a, double[] b, double[] c, String s) {
		this(a, b, c, new Material(1, 0, "vacuum"), s);
	}

	/**
	 * Methode setzt das Material aus dem das Dreieck besteht
	 */

	public void setFacet_Normal(double[] fn) {
		this.facet = fn;
	}

	public double[] getFacet_Normal() {

		return this.facet;
	}

	public void setMaterial(Material m) {
		this.material = m;
	}

	public String getMaterial_data() {
		return this.material.getProperties();
	}

	/**
	 * Methode setzt den translationvector zur Realisierung von periodischen Randbedingungen
	 */
	public void setTranslationvector(double[] t) {
		this.translationvector = t;
	}

	// setzen der Triangle-Temperature

	public void setTemperature(double temp) {
		this.Temperature = temp;
	}

	public void setDistanceInterface(double temp) {
		this.DistanceInterface = temp;
	}

	public void setOFPressure(double temp) {
		this.OFPressure = temp;
	}

	public double getOFPressure() {

		return this.OFPressure;
	}

	// setzen der Verdampfungsleistung

	public void setDeltaP(double temp) {
		this.delta_P += temp;
	}

	public void setKruemmRadius(double temp) {
		this.Kruemmungsradius = temp;
	}

	public double getKruemmRadius() {

		return this.Kruemmungsradius;
	}

	public double getDistanceInterface() {

		return this.DistanceInterface;
	}

	public void setSurfaceTension(double temp) {
		this.surfaceTension = temp;
	}

	public double getSurfaceTension() {

		return this.surfaceTension;
	}

	/**
	 * Methode setzt den translationvector zur Realisierung von periodischen Randbedingungen
	 */
	public void setTranslationvector(double dx, double dy, double dz) {
		this.translationvector = new double[] { dx, dy, dz };
	}

	/**
	 * Erhoehe die Intensitaet um den uebergebenden Wert
	 */
	public void incPower(double in) {
		this.P += in;
	}
	
	public void incVec_k(double [] beam_vec_k) {
		this.vec_k[0] += beam_vec_k[0];
		this.vec_k[1] += beam_vec_k[1];
		this.vec_k[2] += beam_vec_k[2];
	}

	public double getPower() {

		return this.P;
	}
	
	public double [] getSumVec_k() {

		return this.vec_k;
	}
	
	public double getDeltaP() {

		return this.delta_P;
	}

	public void setPower(double in) {

		this.P = in;

	}

	public void setDeltaP_init(double in) {

		this.delta_P = in;

	}

	// Rückgelesene Leistung von OpenFOam
	public void setPower_OF(double Power_OF) {

		this.Power_OF = Power_OF;

	}

	public double getPower_OF() {

		return this.Power_OF;
	}

	/**
	 * gibt die Einfallende Intensitaet zurueck Einfallende Leistung geteilt durch Dreiecksflaeche
	 */
	public double getIntensity() {
		return this.P / (0.5 * Algorithmen.norm2(Algorithmen.xproduct(v1, v2)));
	}

	public double getTemperature() {
		return this.Temperature;
	}

	/**
	 * Methode gibt eine Matrix zurueck, welche in den ersten beiden Spalten die aufspannenden Vektoren des Dreiecks und in der dritten Spalte den Normalenvektor des Dreiecks beinhaltet.
	 */
	public double[][] getBase() {
		final double[][] m = new double[3][3];
		for (int zeile = 0; zeile < 3; zeile++) {
			m[zeile][0] = v1[zeile];
			m[zeile][1] = v2[zeile];
			m[zeile][2] = nv[zeile];
		}
		return m;
	}

	/**
	 * Methode gibt ein Array mit allen 3 Eckpunkten des Dreiecks zurueck In jeder Arrayspalte steht ein Eckpunkt
	 */
	public double[][] getCorners() {
		final double[][] c = new double[3][3];
		for (int i = 0; i < 3; i++) {
			c[i][0] = p1[i];
		}
		for (int i = 0; i < 3; i++) {
			c[i][1] = p2[i];
		}
		for (int i = 0; i < 3; i++) {
			c[i][2] = p3[i];
		}
		return c;
	}

	/**
	 * Methode gibt ein Array mit allen 3 Kanten des Dreiecks zurueck
	 */
	public double[][][] getEdges() {
		final double[][] c = Algorithmen.transpose(getCorners());
		final double[][][] ed = new double[3][3][2];
		for (int i = 0; i < 3; i++) {
			ed[0][i][0] = c[0][i];
			ed[0][i][1] = c[1][i];
		}
		for (int i = 0; i < 3; i++) {
			ed[1][i][0] = c[0][i];
			ed[1][i][1] = c[2][i];
		}
		for (int i = 0; i < 3; i++) {
			ed[2][i][0] = c[1][i];
			ed[2][i][1] = c[2][i];
		}
		return ed;
	}
	public void setDeltaPressure(double p){
		DeltaPressure = p;
	}
	public void setRecoilPressure(double p){
		RecoilPressure = p;
	}
	
	public void setVerschiebefaktor(double p){
		Verschiebefaktor = p;
	}
	
	public void setVerschiebung(double p){
		Verschiebung = p;
	}
	/**
	 * Methode gibt einen String mit wichtigen Eigenschaften des Dreiecks zurueck
	 */
	@Override
	public String toString() {
		String s = new String();
		s += "Dreieck: " + name + "\n";
		s += "P1 " + p1[0] + "," + p1[1] + "," + p1[2] + "\n";
		s += "P2 " + p2[0] + "," + p2[1] + "," + p2[2] + "\n";
		s += "P3 " + p3[0] + "," + p3[1] + "," + p3[2] + "\n";
		s += "v1 " + v1[0] + "," + v1[1] + "," + v1[2] + "\n";
		s += "v2 " + v2[0] + "," + v2[1] + "," + v2[2] + "\n";
		s += "nv " + nv[0] + "," + nv[1] + "," + nv[2] + "\n";
		return s;
	}
}
