/**
 * Objekte dieser Klasse stellen Materialien dar, aus denen die Dreiecke bestehen und welche fuer die Reflektivitaetseigenschaften verantwortlich sind Es werden der Real- und Imaginaerteil des Brechungsindex benoetigt
 */
public class Material {
	public double n; // Realteil des Brechungsindex n - i*k
	public double k; // Imaginaerteil des Brechungsindex n - i*k
	private double kappa; // Hilfgroesse
	private double Ha; // Hilfgroesse
	private double Hb; // Hilfgroesse
	private double Hc; // Hilfgroesse
	// added FFE for multimaterial
	double mat_Hev = 6360000; // j/kg
	double mat_m_atom = 56 * 1.660538921e-27; // kg
	double mat_Tv = 3300; // K

	String name; // Bezeichnung des Materials

	public Material(double n, double k, String s) {
		this.n = n;
		this.k = k;
		name = s;
		calculateVariables();
	}

	public Material(double n, double k, String s, double mat1_Hev, double mat1_m_atom, double mat1_Tv) {
		this.n = n;
		this.k = k;
		this.mat_Hev = mat1_Hev;
		this.mat_m_atom = mat1_m_atom;
		this.mat_Tv = mat1_Tv;
		name = s;
		calculateVariables();
	}

	/*
	* Methode berechnet Hilfsvariablen, welche fuer die Berechnung
	* der winkelabhaengigen Reflektivitaet benoetigt werden.
	*/
	private void calculateVariables() {
		kappa = k / n;
		Ha = n * n - n * n * kappa * kappa;
		Hb = 4 * n * n * n * n * kappa * kappa;
		Hc = -n * n + n * n * kappa * kappa;
	}

	public String getProperties() {

		return ("n=" + this.n + " k=" + this.k + " name=" + this.name);

	}

	/*
	* Methode gibt die Reflektivitaet der senkrechten und parallelen Komponente
	* als 2-dimensionales Array zurueck. Alte Methode, die falsch ist.
	*/
	/*
	public double[] getRsRp(double theta) //Winkel in rad
	{
	  double sint2 = Math.pow(Math.sin(theta),2);
	  double Ha_sint2_2 = Math.pow(Ha-sint2,2);
	  double nt = Math.sqrt(0.5*(Ha+sint2+Math.sqrt(Hb+Ha_sint2_2)));
	  double kt = Math.sqrt(0.5*(Hc+sint2+Math.sqrt(Hb+Ha_sint2_2))/(n*n));
	  double k_t1 = nt*kt;
	  double n_t1 = nt;
	  double cos_t1 = Math.cos(theta);
	  double cos_t2 = Math.cos(Math.asin(Math.sin(theta)/n_t1));

	  double Rs = (Math.pow(cos_t1-n_t1*cos_t2,2) + Math.pow(k_t1*cos_t2,2))/(Math.pow(cos_t1+n_t1*cos_t2,2) + Math.pow(k_t1*cos_t2,2));
	  double Rp = (Math.pow(n_t1*cos_t1-cos_t2,2) + Math.pow(k_t1*cos_t1,2))/(Math.pow(n_t1*cos_t1+cos_t2,2) + Math.pow(k_t1*cos_t1,2));

	  return new double[]{Rs,Rp};
	}
	*/

	/*
	* Methode gibt die Amplitudenreflektivitaet der senkrechten und parallelen Komponente
	* als 2-dimensionales Array zurueck. Die Formel stammt aus "Born/Wolf: Optik" fuer Metallreflexion
	* In die Fresnelformel wird dazu einfach der komplexe Brechungsindex eingesetzt
	*/
	public double[] getrsrp(double x) // Winkel in rad
	{
		final double kappa = k / n;
		final double n2 = n * n;
		final double kappa2 = kappa * kappa;
		final double u2 = Math.sqrt((n2 * (1 - kappa2) - Math.pow(Math.sin(x), 2) + Math.sqrt(Math.pow(n2 * (1 - kappa2) - Math.pow(Math.sin(x), 2), 2) + 4 * n2 * n2 * kappa2)) / 2);
		final double v2 = Math.sqrt((-(n2 * (1 - kappa2) - Math.pow(Math.sin(x), 2)) + Math.sqrt(Math.pow(n2 * (1 - kappa2) - Math.pow(Math.sin(x), 2), 2) + 4 * n2 * n2 * kappa2)) / 2);
		return new double[] { Math.sqrt((Math.pow(Math.cos(x) - u2, 2) + v2 * v2) / (Math.pow(Math.cos(x) + u2, 2) + v2 * v2)), Math.sqrt((Math.pow(n2 * (1 - kappa2) * Math.cos(x) - u2, 2) + Math.pow(2 * n2 * kappa * Math.cos(x) - v2, 2)) / (Math.pow(n2 * (1 - kappa2) * Math.cos(x) + u2, 2) + Math.pow(2 * n2 * kappa * Math.cos(x) + v2, 2))) };
	}
}
