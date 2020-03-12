/**
 * Objekte dieser Klasse stellen Lichtstrahlen dar, welche von Dreiecken reflektiert werden. Diese Strahlen besitzen eine bestimmte "Intensitaet" und eine Polarisationsrichtung. Die Strahlen koennen sowohl Geraden sein, als auch eine Kaustik aufweisen.
 *
 * Fuer gerade Strahlen lautet die Gleichung Punkt = ov + z * rv
 *
 * Fuer Strahlen mit Kaustik lautet die Gleichung Punkt = ov + z * rv + k * w0 * nzs * sqrt(1 + (z - ztt)**2/z0**2) Dabei sind: ov : Ortsvektor des Zentralstrahls rv : Normierter (auf Laenge 1) Richtungsvektor des Zentralstrahls k : Relative Distanz des Strahls vom Zentralstrahl, k=1 bedeutet Strahlradius-Strahl w0 : Radius der Strahltaille nzs: Auf Laenge 1 normierter Vektor, senkrecht auf rv
 * ztt: Distanz auf dem Zentralstrahl vom Startpunkt ov bis zur Strahltaille z0 : Rayleighlaenge
 */

public class Ray

{
	double[] ov; // Ortsvektor des Lichtstrahls
	double[] rv; // Richtungsvektor des Lichtstrahls
	double[] p1; // Polarisationsvektor 1 des Lichtstrahls, steht senkrecht auf der Ausbreitungsrichtung und hat die Laenge des 1. E-Feldvektors
	double[] p2; // Polarisationsvektor 2 des Lichtstrahls, steht senkrecht auf der Ausbreitungsrichtung und hat die Laenge des 2. E-Feldvektors
	double I; // Intensitaet des Lichtstrahls
	// Folgende Groessen sind speziell fuer Strahlen mit Kaustik
	double k; // Relative Distanz des Strahls vom Zentralstrahl, k=1 bedeutet Strahlradius-Strahl
	double w0; // Radius der Strahltaille
	double[] nzs;// Auf Laenge 1 normierter Vektor, senkrecht auf rv
	double ztt;// Distanz auf dem Zentralstrahl vom Startpunkt ov bis zur Strahltaille
	double z0; // Rayleighlaenge
	double F; // Fluenz des Strahls

	String name; // Bezeichnung des Strahls
	boolean kaustik; // true bedeutet, Strahl hat Kaustik, false bedeutet: normaler gerader Strahl

	// Folgender Konstruktor erzeugt Strahl ohne Kaustik
	public Ray(double[] ov, double[] rv, double[] p1, double[] p2, String s) {
		this.ov = ov;
		this.rv = rv;
		this.p1 = p1;
		this.p2 = p2;
		this.name = s;
		kaustik = false;
		calculateIntensity();
	}

	// Folgender Konstruktor erzeugt Strahl mit Kaustik
	// Vorsicht! Achte gegebenenfalls darauf, dass pol1 und pol2 auf 1 normiert sind
	public Ray(double[] ov, double[] rv, double k, double w0, double[] nzs, double ztt, double z0, double I, double[] pol1, double[] pol2, String s) {
		this.ov = ov;
		this.rv = Algorithmen.getNormTox(rv, 1);
		this.k = k;
		this.w0 = w0;
		this.nzs = Algorithmen.getNormTox(nzs, 1);
		this.ztt = ztt;
		this.z0 = z0;
		this.name = s;
		this.p1 = pol1;
		this.p2 = pol2;
		setPower(I);
		kaustik = true;
	}
	
	// Folgender Konstruktor erzeugt gepulste Laserstrahlen mit Kaustik (neuer Paramete ist die Fluenz)
	// Vorsicht! Achte gegebenenfalls darauf, dass pol1 und pol2 auf 1 normiert sind
	public Ray(double[] ov, double[] rv, double k, double w0, double[] nzs, double ztt, double z0, double I, double[] pol1, double[] pol2, double F, String s) 
	{
		this.ov = ov;
		this.rv = Algorithmen.getNormTox(rv, 1);
		this.k = k;
		this.w0 = w0;
		this.nzs = Algorithmen.getNormTox(nzs, 1);
		this.ztt = ztt;
		this.z0 = z0;
		this.name = s;
		this.p1 = pol1;
		this.p2 = pol2;
		setPower(I);
		this.F = F;
		kaustik = true;
	}
	
	
	/**
	 * Methode berechnet aus den Polarisationsvektoren die aktuelle Intensitaet
	 */
	public void calculateIntensity() {
		this.I = Algorithmen.norm2_square(p1) + Algorithmen.norm2_square(p2);
	}

	/**
	 * Methode setzt die Leistung des Strahls
	 */
	public void setPower(double I) {
		final double p1q_p2q = Algorithmen.norm2_square(p1) + Algorithmen.norm2_square(p2);
		p1 = Algorithmen.product(Math.sqrt(I / p1q_p2q), p1);
		p2 = Algorithmen.product(Math.sqrt(I / p1q_p2q), p2);
		this.I = I;
	}

	/**
	 * Methode gibt einen String mit wichtigen Eigenschaften des Strahls zurueck
	 */
	@Override
	public String toString() {
		String s = new String();
		s += "Strahl: " + name + "\n";
		s += "ov " + ov[0] + "," + ov[1] + "," + ov[2] + "\n";
		s += "rv " + rv[0] + "," + rv[1] + "," + rv[2] + "\n";
		s += "I  " + I + "\n";
		return s;
	}
}
