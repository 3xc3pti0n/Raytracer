import static java.nio.file.StandardCopyOption.REPLACE_EXISTING;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;
import java.io.PrintWriter;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.nio.file.StandardCopyOption;
import java.nio.file.StandardOpenOption;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.concurrent.TimeUnit;

import org.apache.commons.lang3.math.NumberUtils;

/**
 * In dieser Klasse sind einige nuetzliche Algorithmen.
 */
public class Algorithmen {
	// Schätzt den mittleren Krümmungsradius in m pro Patch ab
	public static void calculateSurfaceTension(TriangleSet triSet, double surfaceTensionCoeff) {
		// double Kruemmung;
		double meanRadius = 0;
		final ArrayList<Integer> neighbours = new ArrayList<Integer>();

		for (int i = 0; i < triSet.getNumber(); i++) {
			meanRadius = 0;
			neighbours.clear();

			for (int h = 0; h < triSet.getNumber(); h++) {

				if (((Arrays.equals(triSet.triangles[i].p1, triSet.triangles[h].p1) || Arrays.equals(triSet.triangles[i].p1, triSet.triangles[h].p2) || Arrays.equals(triSet.triangles[i].p1, triSet.triangles[h].p3)) && ((Arrays.equals(triSet.triangles[i].p2, triSet.triangles[h].p1) || Arrays.equals(triSet.triangles[i].p2, triSet.triangles[h].p2) || Arrays.equals(triSet.triangles[i].p2,
						triSet.triangles[h].p3)) || (Arrays.equals(triSet.triangles[i].p3, triSet.triangles[h].p1) || Arrays.equals(triSet.triangles[i].p3, triSet.triangles[h].p2) || Arrays.equals(triSet.triangles[i].p3, triSet.triangles[h].p3))))
						|| ((Arrays.equals(triSet.triangles[i].p2, triSet.triangles[h].p1) || Arrays.equals(triSet.triangles[i].p2, triSet.triangles[h].p2) || Arrays.equals(triSet.triangles[i].p2, triSet.triangles[h].p3)) && ((Arrays.equals(triSet.triangles[i].p1, triSet.triangles[h].p1) || Arrays.equals(triSet.triangles[i].p1, triSet.triangles[h].p2) || Arrays.equals(triSet.triangles[i].p1,
								triSet.triangles[h].p3)) || (Arrays.equals(triSet.triangles[i].p3, triSet.triangles[h].p1) || Arrays.equals(triSet.triangles[i].p3, triSet.triangles[h].p2) || Arrays.equals(triSet.triangles[i].p3, triSet.triangles[h].p3))))
						|| ((Arrays.equals(triSet.triangles[i].p3, triSet.triangles[h].p1) || Arrays.equals(triSet.triangles[i].p3, triSet.triangles[h].p2) || Arrays.equals(triSet.triangles[i].p3, triSet.triangles[h].p3)) && ((Arrays.equals(triSet.triangles[i].p1, triSet.triangles[h].p1) || Arrays.equals(triSet.triangles[i].p1, triSet.triangles[h].p2) || Arrays.equals(triSet.triangles[i].p1,
								triSet.triangles[h].p3)) || (Arrays.equals(triSet.triangles[i].p2, triSet.triangles[h].p1) || Arrays.equals(triSet.triangles[i].p2, triSet.triangles[h].p2) || Arrays.equals(triSet.triangles[i].p2, triSet.triangles[h].p3))))) {

					if (!(Arrays.equals(triSet.triangles[i].p1, triSet.triangles[h].p1) && Arrays.equals(triSet.triangles[i].p2, triSet.triangles[h].p2) && Arrays.equals(triSet.triangles[i].p3, triSet.triangles[h].p3))) {
						neighbours.add(h);
					}

				}

			}
			// Berechnen des Krümmungsradius pro patch
						
			double alpha = 0;
			double alpha_check = 0;
			double alpha_new = 0;
			double beta[] = new double[3];
			double dist = 0;
			double normfactor = 0;

			double helper;

			for (final int j : neighbours) {

				beta = Vector_add_substr(triSet.triangles[i].Schwerpunkt(), triSet.triangles[j].Schwerpunkt(), false); // Vektor von (Scherpunkt_i zu Schwerpunkt_j) und Flächennormale von j dürfen keinen Winkel grüßer 90° einschließen, sonst gilt für alpha = 180°-alpha;
				normfactor = Math.sqrt((beta[0] * beta[0]) + (beta[1] * beta[1]) + (beta[2] * beta[2]));
				beta[0] = beta[0] / normfactor;
				beta[1] = beta[1] / normfactor;
				beta[2] = beta[2] / normfactor;
				alpha_check = Algorithmen.product(triSet.triangles[i].getNormalizedFacetNormal(), beta);
				if (alpha_check > 1) {
					alpha_check = 1;
				}
				if (alpha_check < -1) {
					alpha_check = -1;
				}
				alpha_check = Math.acos(alpha_check);
				if (Math.abs(alpha_check) > (Math.PI / 2)) {
					helper = Algorithmen.product(triSet.triangles[i].getNormalizedFacetNormal(), triSet.triangles[j].getNormalizedFacetNormal());
					if (helper > 1) {
						helper = 1;
					}
					if (helper < -1) {
						helper = -1;
					}
					alpha_new = Math.acos(helper);
				}
				else {
					helper = Algorithmen.product(triSet.triangles[i].getNormalizedFacetNormal(), triSet.triangles[j].getNormalizedFacetNormal());
					if (helper > 1) {
						helper = 1;
					}
					if (helper < -1) {
						helper = -1;
					}
					alpha_new = -Math.acos(helper);
				}
				// alpha_check=Math.acos(Algorithmen.product(triSet.triangles[i].getNormalizedFacetNormal(),triSet.triangles[j].getNormalizedFacetNormal())

				alpha = alpha + alpha_new;
				dist = dist + Algorithmen.distance(triSet.triangles[i].Schwerpunkt(), triSet.triangles[j].Schwerpunkt());

				// meanRadius=meanRadius+ (2*Math.sin(alpha_new/2))/Algorithmen.distance(triSet.triangles[i].Schwerpunkt(),triSet.triangles[j].Schwerpunkt());

				if (Double.isInfinite(alpha)) {
					alpha = 1e5;
				}
				// if(alpha==0){alpha=1e-5;}

			}
			alpha = alpha / neighbours.size();
			dist = dist / neighbours.size();

			// meanRadius=meanRadius/neighbours.size();
			// meanRadius=1/meanRadius;
			meanRadius = dist / (2 * Math.sin(alpha / 2));

			// System.out.println(1/meanRadius);

			// Kruemmung=((double)1/neighbours.size())*(2*Math.sin(alpha)/dist);
			// meanRadius=(double)1.0/Kruemmung;
			// meanRadius=dist/(2*Math.sin(alpha/2));

			meanRadius = meanRadius / 1000D; // Umrechnen von mm nach m
			// triSet.triangles[i].setKruemmRadius(1/Kruemmung);

			triSet.triangles[i].setSurfaceTension(surfaceTensionCoeff * (2.0 / meanRadius));
			// System.out.println(triSet.triangles[i].getSurfaceTension());

		}

	}

	public static double[] Vector_add_substr(double[] V1, double[] V2, boolean add) {
		final double[] V3 = new double[3];

		if (add) {
			V3[0] = V1[0] + V2[0];
			V3[1] = V1[1] + V2[1];
			V3[2] = V1[2] + V2[2];
		}
		else {
			V3[0] = V2[0] - V1[0];
			V3[1] = V2[1] - V1[1];
			V3[2] = V2[2] - V1[2];
		}

		return V3;
	}

	public static double LeistungOpenFOAM;

	// Methode zur Korrektur der eingekoppelten Leistung

	// Methode kopiert FOAMDicts
	public static void copyFoamDicts(String fileToCopy, String TargetFile) {
		try {

			final java.nio.file.Path source = java.nio.file.Paths.get(fileToCopy);
			final java.nio.file.Path destination = java.nio.file.Paths.get(TargetFile);

			java.nio.file.Files.copy(source, destination, StandardCopyOption.COPY_ATTRIBUTES, StandardCopyOption.REPLACE_EXISTING);
		}
		catch (final Exception e) {
			System.out.println("Error in Algorithmen.copyFoamDicts, source:  + ");
		}

	}

	public static void copy_file(String fileToCopy, String TargetFile) {
		try {
			final java.nio.file.Path source = java.nio.file.Paths.get(fileToCopy);
			final java.nio.file.Path destination = java.nio.file.Paths.get(TargetFile);

			java.nio.file.Files.copy(source, destination, StandardCopyOption.COPY_ATTRIBUTES, StandardCopyOption.REPLACE_EXISTING);
		}
		catch (final IOException e) {
			System.out.println("Error in Algorithmen.copy_file");
		}

	}

	public static String get_OF_parameter(File file, String parameter) {
		String value = "";
		// System.out.println("get Parameter");
		try (BufferedReader br = new BufferedReader(new FileReader(file))) { // ASCII STL öffnen
			String line;
			while (((line = br.readLine()) != null) && (!line.trim().equals("END"))) { // Für jede line
				// System.out.println(line);
				line = line.trim();
				if (line.startsWith(parameter)) {
					// System.out.println(line);
					line = line.split(";")[0];
					// System.out.println(line);
					final String[] line2 = line.split(" ");
					// System.out.println(line2);
					value = line2[line2.length - 1];
					// System.out.println(value);
				}
			}
			br.close();
		}
		catch (final FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		catch (final IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}

		return value;
	}

	// Methode zum Speichern der OpenFoam Dicts

	public static boolean check_normal_orient(Triangle tri) {

		// Check der Orientierung der Normalen

		double[] nv;

		nv = xproduct(tri.v1, tri.v2);

		// normfaktor = Math.sqrt((tri.nv[0]*tri.nv[0])+(tri.nv[1]*tri.nv[1])+(tri.nv[2]*tri.nv[2]));

		if (Math.signum(nv[0]) == Math.signum(tri.facet[0]) && Math.signum(nv[1]) == Math.signum(tri.facet[1]) && Math.signum(nv[2]) == Math.signum(tri.facet[2])) {

			return false;

		}
		else {

			return true;

		}

	}

	public static void corrected_normal(Triangle tri, boolean normalen_invert) {

		double[] nv;
		double normfaktor;
		final double[] helper = new double[3];

		tri.recalculate_triangle_vector();

		if (!normalen_invert) {
			nv = xproduct(tri.v1, tri.v2);
			tri.nv = nv;
			normfaktor = Math.sqrt((tri.nv[0] * tri.nv[0]) + (tri.nv[1] * tri.nv[1]) + (tri.nv[2] * tri.nv[2]));
			// tri.facet[0]=nv[0]/normfaktor;
			// tri.facet[1]=nv[1]/normfaktor;
			// tri.facet[2]=nv[2]/normfaktor;
			helper[0] = nv[0] / normfaktor;
			helper[1] = nv[1] / normfaktor;
			helper[2] = nv[2] / normfaktor;
			tri.setFacet_Normal(helper);

		}
		else {
			nv = xproduct(tri.v2, tri.v1);
			tri.nv = nv;
			normfaktor = Math.sqrt((tri.nv[0] * tri.nv[0]) + (tri.nv[1] * tri.nv[1]) + (tri.nv[2] * tri.nv[2]));
			// tri.facet[0]=nv[0]/normfaktor;
			// tri.facet[1]=nv[1]/normfaktor;
			// tri.facet[2]=nv[2]/normfaktor;
			helper[0] = nv[0] / normfaktor;
			helper[1] = nv[1] / normfaktor;
			helper[2] = nv[2] / normfaktor;
			tri.setFacet_Normal(helper);
		}

	}

	public static double[] get_normalized_FacetNormal(Triangle tri, boolean normalen_invert) {

		double[] nv;
		double normfaktor;
		final double[] helper = new double[3];

		if (!normalen_invert) {
			nv = xproduct(tri.v1, tri.v2);
			tri.nv = nv;
			normfaktor = Math.sqrt((tri.nv[0] * tri.nv[0]) + (tri.nv[1] * tri.nv[1]) + (tri.nv[2] * tri.nv[2]));
			// tri.facet[0]=nv[0]/normfaktor;
			// tri.facet[1]=nv[1]/normfaktor;
			// tri.facet[2]=nv[2]/normfaktor;
			helper[0] = nv[0] / normfaktor;
			helper[1] = nv[1] / normfaktor;
			helper[2] = nv[2] / normfaktor;
			return helper;

		}
		else {
			nv = xproduct(tri.v2, tri.v1);
			tri.nv = nv;
			normfaktor = Math.sqrt((tri.nv[0] * tri.nv[0]) + (tri.nv[1] * tri.nv[1]) + (tri.nv[2] * tri.nv[2]));
			// tri.facet[0]=nv[0]/normfaktor;
			// tri.facet[1]=nv[1]/normfaktor;
			// tri.facet[2]=nv[2]/normfaktor;
			helper[0] = nv[0] / normfaktor;
			helper[1] = nv[1] / normfaktor;
			helper[2] = nv[2] / normfaktor;
			return helper;

		}

	}

	// Diese Methode bestimmt den Verschiebungsfaktor .. soll später durch Verdampfungsdruck berechnet werden

	public static double Verschiebungsfunktion(Triangle tri, double T_vapor, double T_max, double delta_global, double tol) {

		double delta_s;

		if (T_max > T_vapor) {
			delta_s = ((delta_global) / (T_max - T_vapor)) * (tri.getTemperature() - T_vapor);
		}
		else {
			delta_s = 0;
		}

		if (tri.getTemperature() < T_vapor * (1 + tol) && tri.getTemperature() > T_vapor * (1 - tol)) {
			delta_s = 0;
		}

		return delta_s;

	}

	public static double Verschiebungsfunktion_p(Triangle tri, double T_melt, double delta_Pmax, double delta_global, double tol) {

		//delta_Pmax= 1e5, Verschiebungsfaktor_global = 0.05
		
		// Verdampfungsmodell Philipp
		/*
		  double roh = 2;
		  double Hev=6213627;
		  double A_patch=0.5*Algorithmen.norm2(Algorithmen.xproduct(tri.v1,tri.v2))/1000000D;
		  double q_v=tri.getPower()-tri.getPower_OF();
		  double v=0;
		  double p_recoil=0;
		/*
		  if(q_v>=0){
			  v=Math.pow((Math.pow(((8*Math.pow(Hev,3))/27+Math.pow(q_v,2)/(Math.pow(A_patch,2)*Math.pow(roh,2))),(double)1/2)+q_v/(A_patch*roh)),(double)1/3)-(2*Hev)/(3*Math.pow((Math.pow(((8*Math.pow(Hev,3))/27+Math.pow(q_v,2)/(Math.pow(A_patch,2)*Math.pow(roh,2))),(double)1/2)+q_v/(A_patch*roh)),(double)1/3));
			   p_recoil=roh*Math.pow(v, 2);
		  }else{
		//	   q_v=q_v*(-1);
		//	   v=Math.pow((Math.pow(((8*Math.pow(Hev,3))/27+Math.pow(q_v,2)/(Math.pow(A_patch,2)*Math.pow(roh,2))),(double)1/2)+q_v/(A_patch*roh)),(double)1/3)-(2*Hev)/(3*Math.pow((Math.pow(((8*Math.pow(Hev,3))/27+Math.pow(q_v,2)/(Math.pow(A_patch,2)*Math.pow(roh,2))),(double)1/2)+q_v/(A_patch*roh)),(double)1/3));
		//	   p_recoil=-roh*Math.pow(v, 2);
			  p_recoil=0;
		  }

		//  double p_recoil=roh*Math.pow(v, 2);
		 */

		// Verdampfungsmodell Chen Wang

		final double Hev = 6360000;
		final double A_patch = 0.5 * Algorithmen.norm2(Algorithmen.xproduct(tri.v1, tri.v2)) / 1000000D;
		final double q_v = tri.getPower() - tri.getPower_OF();
		tri.delta_P = q_v;
		double p_recoil;
		final double m_atom = 56 * 1.660538921e-27;
		final double kB = 1.3806488e-23;
		final double Tv = 3300;

		final double relaxFactor = 0.5;

		p_recoil = (relaxFactor * q_v / (A_patch * Hev)) * Math.sqrt(Math.PI * kB * Tv / (2 * m_atom));

		if (p_recoil < 0) {
			p_recoil = 0;
		}

		// if ((tri.p1[2]>0.1 && tri.p2[2]>0.1 && tri.p3[2]>0.1 && tri.getPower_OF()<1e-10) || q_v<0 ){p_recoil=0;}

		double delta_p = 0;
		double delta_s = 0;

		// Berechnung der Druckbilanz über Temperaturgradient an Oberfläche + ChenWang + SurfaceTension
		//

		// delta_p= p_recoil-tri.getSurfaceTension()-tri.getOFPressure()*7000; //7000 kg/m^3 (z.B Stahl) - für scalartransportFoam "ab-tri... weg"

		delta_p = p_recoil - tri.getSurfaceTension();
		
		// Berechnung der Druckbilanz über intern berechneten Rückstoßdruck in OpenFOAM + SurfaceTension

		// delta_p= -tri.getSurfaceTension()-tri.getOFPressure()*7000; // Rückstoßdruck wird in aktueller Version im Solver direkt berechnet --> OFPressure sind dynamischer +- statische +- Rückstoßdruck

		// System.out.println("OFPressure:  "+tri.getOFPressure());
		// delta_p= tri.getSurfaceTension();

		delta_s = ((delta_global) / (delta_Pmax)) * (delta_p);
		
		System.out.println(delta_p);
		System.out.println(delta_global);
		System.out.println(delta_s);
		// ab hier: auskommentieren bei scalartransport...
		if (delta_s > tri.getDistanceInterface()) {
			delta_s = tri.getDistanceInterface();
		}
		// bis hier

		// if(tri.getTemperature()<T_melt && delta_s>0){delta_s=0;}

		// System.out.println("P_ST= "+tri.getSurfaceTension());
		// Bound of variables
		if (delta_s > delta_global) {
			delta_s = delta_global;
		}
		if (delta_s < -delta_global) {
			delta_s = -delta_global;
		}
		tri.setDeltaPressure(delta_p);
		tri.setRecoilPressure(p_recoil);
		tri.setVerschiebefaktor(delta_s);
		System.out.println(delta_s);
		return delta_s;

	}

	// Diese Methode setzt den temperaturabhängigen komplexen Brechungsindex

	public static void setRefractionIndex(TriangleSet triSet) {

		// Angabe von n Stützpunkten, dazwischen lineare Interpolation

		final double[] Stuetzp_n1 = new double[2];
		final double[] Stuetzp_n2 = new double[2];
		final double[] Stuetzp_n3 = new double[2];
		final double[] Stuetzp_n4 = new double[2];

		final double[] Stuetzp_k1 = new double[2];
		final double[] Stuetzp_k2 = new double[2];
		final double[] Stuetzp_k3 = new double[2];
		final double[] Stuetzp_k4 = new double[2];

		// n(T)
		Stuetzp_n1[0] = 0.3;
		Stuetzp_n2[0] = 0.5;
		Stuetzp_n3[0] = 1.2;
		Stuetzp_n4[0] = 1.6;

		// T
		Stuetzp_n1[1] = 300;
		Stuetzp_n2[1] = 1300;
		Stuetzp_n3[1] = 1400;
		Stuetzp_n4[1] = 3200;

		// k(T)
		Stuetzp_k1[0] = 6.7;
		Stuetzp_k2[0] = 6.5;
		Stuetzp_k3[0] = 7.2;
		Stuetzp_k4[0] = 6.5;

		// T
		Stuetzp_k1[1] = Stuetzp_n1[1];
		Stuetzp_k2[1] = Stuetzp_n2[1];
		Stuetzp_k3[1] = Stuetzp_n3[1];
		Stuetzp_k4[1] = Stuetzp_n4[1];

		double n_T;
		double k_T;

		double Temp;

		// lineare Interpolation mit Zuweisung der Werte zu den Triangles

		for (int k = 0; k < triSet.getNumber(); k++) {

			Temp = triSet.triangles[k].getTemperature();

			if (Temp < Stuetzp_n1[1]) {
				final Material m = new Material(Stuetzp_n1[0], Stuetzp_k1[0], "mat");
				triSet.triangles[k].setMaterial(m);
			}

			if (Temp > Stuetzp_n1[1] && Temp <= Stuetzp_n2[1]) {

				n_T = Stuetzp_n1[0] + ((Stuetzp_n2[0] - Stuetzp_n1[0]) / (Stuetzp_n2[1] - Stuetzp_n1[1])) * (Temp - Stuetzp_n1[1]);
				k_T = Stuetzp_k1[0] + ((Stuetzp_k2[0] - Stuetzp_k1[0]) / (Stuetzp_k2[1] - Stuetzp_k1[1])) * (Temp - Stuetzp_k1[1]);
				final Material m = new Material(n_T, k_T, "mat");
				triSet.triangles[k].setMaterial(m);
			}

			if (Temp > Stuetzp_n2[1] && Temp <= Stuetzp_n3[1]) {

				n_T = Stuetzp_n2[0] + ((Stuetzp_n3[0] - Stuetzp_n2[0]) / (Stuetzp_n3[1] - Stuetzp_n2[1])) * (Temp - Stuetzp_n2[1]);
				k_T = Stuetzp_k2[0] + ((Stuetzp_k3[0] - Stuetzp_k2[0]) / (Stuetzp_k3[1] - Stuetzp_k2[1])) * (Temp - Stuetzp_k2[1]);
				final Material m = new Material(n_T, k_T, "mat");
				triSet.triangles[k].setMaterial(m);
			}

			if (Temp > Stuetzp_n3[1] && Temp <= Stuetzp_n4[1]) {

				n_T = Stuetzp_n3[0] + ((Stuetzp_n4[0] - Stuetzp_n3[0]) / (Stuetzp_n4[1] - Stuetzp_n3[1])) * (Temp - Stuetzp_n3[1]);
				k_T = Stuetzp_k3[0] + ((Stuetzp_k4[0] - Stuetzp_k3[0]) / (Stuetzp_k4[1] - Stuetzp_k3[1])) * (Temp - Stuetzp_k3[1]);
				final Material m = new Material(n_T, k_T, "mat");
				triSet.triangles[k].setMaterial(m);
			}

			if (Temp > Stuetzp_n4[1]) {
				final Material m = new Material(Stuetzp_n4[0], Stuetzp_k4[0], "mat");
				triSet.triangles[k].setMaterial(m);
			}

			// System.out.println(triSet.triangles[k].getMaterial_data());

			// Material m = new Material(1,1,"mat");
			// triSet.triangles[k].setMaterial(m);

		}
	}

	public static void setRefractionIndex_Server_Window(TriangleSet triSet, double[] T_stuetz, double[][] nk_stuetz, float[] temperatures) {
		// lineare Interpolation mit Zuweisung der Werte zu den Triangles
		double Temp;
		double n_T;
		double k_T;
		for (int k = 0; k < triSet.getNumber(); k++) {
			Temp = triSet.triangles[k].getTemperature();

			if (Temp < T_stuetz[0]) {
				final Material m = new Material(nk_stuetz[0][0], nk_stuetz[0][1], "mat");
				triSet.triangles[k].setMaterial(m);
			}

			if (Temp > T_stuetz[0] && Temp <= T_stuetz[1]) {

				n_T = nk_stuetz[0][0] + ((nk_stuetz[1][0] - nk_stuetz[0][0]) / (T_stuetz[1] - T_stuetz[0])) * (Temp - T_stuetz[0]);
				k_T = nk_stuetz[0][1] + ((nk_stuetz[1][1] - nk_stuetz[0][1]) / (T_stuetz[1] - T_stuetz[0])) * (Temp - T_stuetz[0]);
				final Material m = new Material(n_T, k_T, "mat");
				triSet.triangles[k].setMaterial(m);
			}

			if (Temp > T_stuetz[1] && Temp <= T_stuetz[2]) {

				n_T = nk_stuetz[1][0] + ((nk_stuetz[2][0] - nk_stuetz[1][0]) / (T_stuetz[2] - T_stuetz[1])) * (Temp - T_stuetz[1]);
				k_T = nk_stuetz[1][1] + ((nk_stuetz[2][1] - nk_stuetz[1][1]) / (T_stuetz[2] - T_stuetz[1])) * (Temp - T_stuetz[1]);
				final Material m = new Material(n_T, k_T, "mat");
				triSet.triangles[k].setMaterial(m);
			}

			if (Temp > T_stuetz[2]) {

				n_T = nk_stuetz[2][0];
				k_T = nk_stuetz[2][1];
				final Material m = new Material(n_T, k_T, "mat");
				triSet.triangles[k].setMaterial(m);
			}

		}
	}

	public static void write_combined_output_as_ply(TriangleSet ts1, TriangleSet ts2, String param, double share_fresn, PrintWriter ausgabe) throws IOException {
		double maxI = 0;
		System.out.println("in write_combined_output_as_ply");
		for (int i = 0; i < ts1.getNumber(); i++) {
			ts1.triangles[i].setPower(share_fresn * ts1.triangles[i].getPower() + (1 - share_fresn) * ts2.triangles[i].getPower());
		}
		System.out.println("1: in write_combined_output_as_ply");
		for (int i = 0; i < ts1.getNumber(); i++) {
			final double inten = ts1.triangles[i].getIntensity();
			if (inten > maxI) {
				maxI = inten;
			}
		}
		System.out.println("2: in write_combined_output_as_ply");
		double[] is = new double[ts1.getNumber()];
		for (int i = 0; i < ts1.getNumber(); i++) {
			is[i]= ts1.triangles[i].getIntensity();
		}
		Arrays.sort(is);
		List<Double> intensities = new ArrayList<Double>();
		for (int i = 0; i < ts1.getNumber(); i++) {
			intensities.add(is[i]);
		}
		System.out.println("3: in write_combined_output_as_ply");
		
		List<Double[]> cmap = Algorithmen.read_colormap("cmap_jet_ncar.txt");
		if (ts1 != null) {
			ausgabe.println("ply");
			ausgabe.println("format ascii 1.0");
			ausgabe.println("comment author: IFSW Michalowski Raytracer");
			ausgabe.println("comment object: ");
			ausgabe.println("element vertex " + 3 * ts1.getNumber());
			ausgabe.println("property float x");
			ausgabe.println("property float y");
			ausgabe.println("property float z");
			ausgabe.println("property uchar red");
			ausgabe.println("property uchar green");
			ausgabe.println("property uchar blue");
			ausgabe.println("element face " + ts1.getNumber());
			ausgabe.println("property list uchar int vertex_index");
			ausgabe.println("element edge 0");
			ausgabe.println("property int vertex1");
			ausgabe.println("property int vertex2");
			ausgabe.println("property uchar red");
			ausgabe.println("property uchar green");
			ausgabe.println("property uchar blue");
			ausgabe.println("end_header");

			int R = 0;
			int G = 0;
			int B = 0;
			int cmap_i = 0;
			for (int i = 0; i < ts1.getNumber(); i++) {
				final double[][] corners = ts1.triangles[i].getCorners();

				for (int j = 0; j < 3; j++) {
					final double inten = ts1.triangles[i].getIntensity();
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
						cmap_i = (int) (100*((double) intensities.indexOf(inten))/intensities.size());
						R = (int) ( cmap.get(cmap_i)[0]*255.);
						G = (int) ( cmap.get(cmap_i)[1]*255.);
						B = (int) ( cmap.get(cmap_i)[2]*255.);
						
					}
					ausgabe.println((float) corners[0][j] / 1000000.0f + " " // FFE: nochmal durch 1000 geteilt
							+ (float) corners[1][j] / 1000000.0f + " " + (float) corners[2][j] / 1000000.0f + " " + R + " " + G + " " + B);

				}
			}
			for (int i = 0; i < ts1.getNumber(); i++) {
				ausgabe.println("3 " + 3 * i + " " + (3 * i + 1) + " " + (3 * i + 2));
			}
		}
		ausgabe.close();
	}

	public static String last_OF_timestep(File runDir) {
		final File[] in_case_dir = runDir.listFiles();
		final List<String> timestep_list = new ArrayList<String>();
		timestep_list.add("0");
		for (int kk = 0; kk < in_case_dir.length; kk++) {
			if (Algorithmen.isNumeric(in_case_dir[kk].getName()) > Algorithmen.isNumeric(timestep_list.get(timestep_list.size() - 1))) {
				timestep_list.add(in_case_dir[kk].getName());
			}

		}
		return timestep_list.get(timestep_list.size() - 1);
	}

	
	// Diese Methode liest die Temperatur der einzelnen Patches und ordnet sie den Triangles zu



	public static void delete_projectFolder(File dir) {

		final File[] files = dir.listFiles();

		if (files != null) {
			for (int i = 0; i < files.length; i++) {

				if (files[i].isDirectory()) {

					if (!files[i].getName().equalsIgnoreCase("oldTime") && !files[i].getName().equalsIgnoreCase("parallel") && !files[i].getName().equalsIgnoreCase("Kapillare") && files[i].getName().equalsIgnoreCase("0") && !files[i].getName().equalsIgnoreCase("system") && !files[i].getName().equalsIgnoreCase("OpenFOAM-2.1.x") && !files[i].getName().equalsIgnoreCase("constant")) {

						delete_projectFolder(files[i]);

					}

				}
				if (!files[i].getName().equalsIgnoreCase("laserweldfoam_noRec_2mat.bat") && !files[i].getName().equalsIgnoreCase("write_alpha1_setfields.py") && !files[i].getName().equalsIgnoreCase("generate_boundaries_setFields.bat") && !files[i].getName().equalsIgnoreCase("setFields.bat") && !files[i].getName().equalsIgnoreCase("rhoCentralFoam.bat")
						&& !files[i].getName().equalsIgnoreCase("chtMultiRegionSimpleFoam.bat") && !files[i].getName().equalsIgnoreCase("ptot.bat") && !files[i].getName().equalsIgnoreCase("foamtovtk_cht.bat") && !files[i].getName().equalsIgnoreCase("splitMeshRegions.bat") && !files[i].getName().equalsIgnoreCase("chtMultiRegionFoam.bat") && !files[i].getName().equalsIgnoreCase("setDir.bat")
						&& !files[i].getName().equalsIgnoreCase("patchAverage.bat") && !files[i].getName().equalsIgnoreCase("reconstructPar.bat") && !files[i].getName().equalsIgnoreCase("decomposePar.bat") && !files[i].getName().equalsIgnoreCase("gompi.bat") && !files[i].getName().equalsIgnoreCase("setvars.bat") && !files[i].getName().equalsIgnoreCase("DOS_Mode.bat")
						&& !files[i].getName().equalsIgnoreCase("snappyhexmesh.bat") && !files[i].getName().equalsIgnoreCase("transformpoints.bat") && !files[i].getName().equalsIgnoreCase("blockmesh.bat") && !files[i].getName().equalsIgnoreCase("potentialFoam.bat") && !files[i].getName().equalsIgnoreCase("scalartransportfoam.bat") && !files[i].getName().equalsIgnoreCase("foamtovtk.bat")
						&& !files[i].getName().equalsIgnoreCase("foamCalc.bat") && !files[i].getName().equalsIgnoreCase("buoyantpimplefoam.bat") && !files[i].getName().equalsIgnoreCase("buoyantboussinesqpimplefoam.bat")) {
					files[i].delete();

					/*
					  	   try{
					  	      org.apache.commons.io.FileUtils.deleteDirectory(files[i]);
					  	    }catch(Exception e){e.toString();}
					 */

				}
			}

		}
	}

	public static void copyVTK(File fileToCopy, File TargetFile) {
		try {
			org.apache.commons.io.FileUtils.copyDirectory(fileToCopy, TargetFile);
		}
		catch (final IOException e) {
			System.out.println("Error in Algorithmen.copyVTK, source: " + fileToCopy.getAbsolutePath() + "Target: " + TargetFile.getAbsolutePath());
		}
	}

	public static void copy_file(File fileToCopy, File TargetFile) {
		try {
			Files.copy(Paths.get(fileToCopy.getAbsolutePath()), Paths.get(TargetFile.getAbsolutePath()), REPLACE_EXISTING);
		}
		catch (final IOException e) {
			System.out.println("Error in Algorithmen.copy_file, source: " + fileToCopy.getAbsolutePath() + "Target: " + TargetFile.getAbsolutePath());
		}
	}

	public static void copy_dir(File fileToCopy, File TargetFile) {
		try {
			org.apache.commons.io.FileUtils.copyDirectory(fileToCopy, TargetFile);
		}
		catch (final Exception e) {
			System.out.println("Error in Algorithmen.copy_dir, source: " + fileToCopy.getAbsolutePath() + "  target: " + TargetFile.getAbsolutePath());
		}
	}

	// Boundaries

	public static void setDir(String outputfile, String dir) {

		try {

			final File datei = new File(outputfile);
			final FileWriter ausgabestrom = new FileWriter(datei);
			final PrintWriter ausgabe = new PrintWriter(ausgabestrom);

			ausgabe.println(dir);

			ausgabe.close();
			ausgabestrom.close();

		}
		catch (final Exception e) {
			System.out.println("Einlesen " + e.toString());
		}

	}

	public static void convert_to_multipart_stl(String inputfile, String outputfile) {
		String facet_normal = new String();
		String vertex1 = new String();
		String vertex2 = new String();
		String vertex3 = new String();

		int counter;
		try {

			final BufferedReader br = new BufferedReader(new FileReader(inputfile));
			final File datei = new File(outputfile);
			final FileWriter ausgabestrom = new FileWriter(datei);
			final PrintWriter ausgabe = new PrintWriter(ausgabestrom);

			String line;
			counter = 1;

			while (((line = br.readLine()) != null) && (!line.trim().equals("END"))) {
				String bezeichnung = new String();
				bezeichnung = line;

				if (bezeichnung.contains("facet normal")) {
					facet_normal = bezeichnung;
				}

				if (bezeichnung.contains("outer loop")) {
					vertex1 = br.readLine();
					vertex2 = br.readLine();
					vertex3 = br.readLine();

					ausgabe.println("solid bereich" + counter);
					ausgabe.println(facet_normal);
					ausgabe.println("    " + "outer loop");
					ausgabe.println("      " + vertex1);
					ausgabe.println("      " + vertex2);
					ausgabe.println("      " + vertex3);
					ausgabe.println("    " + "endloop");
					ausgabe.println("  " + "endfacet");
					ausgabe.println("endsolid bereich" + counter);

					counter = counter + 1;

				}

			}
			br.close();
			ausgabe.close();
			ausgabestrom.close();

		}

		catch (final Exception e) {
			System.out.println("Einlesen " + e.toString());
		}

	}

	// Ausführen der Batch-Skripte (OpenFoam)

	public static void callbatchfile(String name, boolean parallel) {
		final int time_limit = 3000000;
		int finished = -999;
		boolean wait_finished;

		// String[] command = {"cmd.exe", "/c","start","DOS_Mode.bat","start",name};
		final String[] command = { "cmd.exe", " /c", name };

		try {

			if (parallel) {
				final Process pb = Runtime.getRuntime().exec(command, null, new File("FEM/parallel/"));
				finished = pb.waitFor();
			}
			else {
				if (name.equals("transformpoints.bat")) {
					final Process pb = Runtime.getRuntime().exec(command, null, new File("FEM/"));
					wait_finished = pb.waitFor(30, TimeUnit.SECONDS);
				}
				final Process pb = Runtime.getRuntime().exec(command, null, new File("FEM/"));
				finished = pb.waitFor();
			}
			// Hier kommt die Rückmeldung vom Batchfile, sobald der CMD-Prompt geschlossen wird.

			if (finished == 0) {
				// System.out.println("Batchfile done");

			}
		}

		catch (final IOException e) {
			System.out.println(e.toString());
			// TODO Auto-generated catch block

		}
		catch (final InterruptedException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
			System.out.println(e.toString());
		}

		// final int finished = pb.waitFor();

		// if(finished==0){
		// pb = Runtime.getRuntime().exec("exit");

		// }

		// ProcessBuilder pb = new ProcessBuilder(command);
		// pb.directory(new File(System.getProperty("user.dir")+ "/FEM"));
		// System.out.println(pb.directory().getAbsolutePath());
		// System.out.println(name);
		// pb.start();

	}

	// Überladene CallBatchFile Methode

	public static void callbatchfile(String name, String option1, String option2, boolean parallel) {

		int finished;

		// String[] command = {"cmd.exe", "/c","start","DOS_Mode.bat","start",name};
		final String[] command = { "cmd.exe", " /c", name, option1, option2 };

		try {

			if (parallel) {
				final Process pb = Runtime.getRuntime().exec(command, null, new File("FEM/parallel/"));
				finished = pb.waitFor();
			}
			else {
				final Process pb = Runtime.getRuntime().exec(command, null, new File("FEM/"));
				finished = pb.waitFor();
			}

			// Hier kommt die Rückmeldung vom Batchfile, sobald der CMD-Prompt geschlossen wird.

			if (finished == 0) {
				// System.out.println("Batchfile done");

			}
		}

		catch (final IOException e) {
			System.out.println(e.toString());
			// TODO Auto-generated catch block

		}
		catch (final InterruptedException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
			System.out.println(e.toString());
		}

	}

	// Überladene CallBatchFile Methode

	public static void callbatchfile(String name, String option1, String option2, String option3, String option4, boolean parallel) {

		int finished;

		// String[] command = {"cmd.exe", "/c","start","DOS_Mode.bat","start",name};
		final String[] command = { "cmd.exe", " /c", name, option1, option2, option3, option4 };

		try {

			if (parallel) {
				final Process pb = Runtime.getRuntime().exec(command, null, new File("FEM/parallel/"));
				finished = pb.waitFor();
			}
			else {
				final Process pb = Runtime.getRuntime().exec(command, null, new File("FEM/"));
				finished = pb.waitFor();
			}

			// Hier kommt die Rückmeldung vom Batchfile, sobald der CMD-Prompt geschlossen wird.

			if (finished == 0) {
				// System.out.println("Batchfile done");

			}
		}

		catch (final IOException e) {
			System.out.println(e.toString());
			// TODO Auto-generated catch block

		}
		catch (final InterruptedException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
			System.out.println(e.toString());
		}

		// final int finished = pb.waitFor();

		// if(finished==0){
		// pb = Runtime.getRuntime().exec("exit");

		// }

		// ProcessBuilder pb = new ProcessBuilder(command);
		// pb.directory(new File(System.getProperty("user.dir")+ "/FEM"));
		// System.out.println(pb.directory().getAbsolutePath());
		// System.out.println(name);
		// pb.start();

	}

	public static void correctChtBoundary(String path, String sampleRegion, String samplePatch_eigen, String samplePatch_gegen) {

		final StringBuilder sb;
		sb = new StringBuilder();
		String Hilfsstring;
		final StringBuilder globalsb = new StringBuilder();
		String patchName;
		try {

			final BufferedReader br = new BufferedReader(new FileReader(path));

			String line;

			while (((line = br.readLine()) != null) && (!line.trim().equals("END"))) {
				String bezeichnung = new String();
				bezeichnung = line;
				globalsb.append(line);
				globalsb.append("\n");

				// Parser für magT --> Einlesen der Temperatur der Patches

				if (bezeichnung.contains("myfile")) {
					globalsb.append(" { ");
					patchName = bezeichnung;
					patchName = patchName.replace(samplePatch_eigen, samplePatch_gegen);

					do {
						sb.append(line);
						sb.append("\n");
						line = br.readLine();

					}
					while ((!line.contains("}")));

					Hilfsstring = sb.toString();
					Hilfsstring = Hilfsstring.replaceAll("\n", " ");
					sb.setLength(0);
					// System.out.println(Hilfsstring);

					final String Hilfsvariable[] = Hilfsstring.split(";");

					for (int i = 0; i < Hilfsvariable.length; i++) {

						if (Hilfsvariable[i].contains("type")) {
							Hilfsvariable[i] = "type mappedWall; sampleMode nearestPatchFace; sampleRegion " + sampleRegion + ";" + "samplePatch " + patchName;
						}

					}

					for (int i = 0; i < Hilfsvariable.length; i++) {

						globalsb.append(Hilfsvariable[i]);
						if (i < Hilfsvariable.length - 1) {
							globalsb.append(" ;");
						}
						globalsb.append("\n");

					}
					globalsb.append(" } ");
				}

			}
			br.close();

		}
		catch (final Exception e) {
			System.out.println("Fehler:" + e.toString());
		}

		try {

			final File datei = new File(path);
			final FileWriter ausgabestrom = new FileWriter(datei);
			final PrintWriter ausgabe = new PrintWriter(ausgabestrom);
			ausgabe.println(globalsb.toString());

			ausgabe.close();
			ausgabestrom.close();

		}
		catch (final Exception e) {
			System.out.println("Einlesen " + e.toString());
		}

		/*
		final StringBuilder sb;
		sb = new StringBuilder();
		String Hilfsstring;
		StringBuilder globalsb = new StringBuilder();
		String patchName;
		String Name;
		int anzPatch=0;
		int anz_gesamt=0;
			 try {


		      BufferedReader br = new BufferedReader(new FileReader(path));


		     String line;


		      while (((line = br.readLine())!=null)&&(!line.trim().equals("END"))){
		        String bezeichnung = new String();
		        bezeichnung = line;

		        // Einlesen von Header
				       if(!line.contains("myfile")){
				        globalsb.append(line);
		   		        globalsb.append("\n");
				       }



		        	if (bezeichnung.contains("myfile")){
		        		   Name = bezeichnung;
		        		   patchName=bezeichnung;
		        		   patchName=patchName.replace(samplePatch_eigen,samplePatch_gegen);
		        		   anz_gesamt+=1;

		        		   do {
		        			   sb.append(line);
		        			   sb.append("\n");
		        			   line=br.readLine();

		        		   }while((!line.contains("}")));

		        		   Hilfsstring = sb.toString();
		        		   Hilfsstring = Hilfsstring.replaceAll("\n"," ");
		        		   sb.setLength(0);
		        		   //  System.out.println(Hilfsstring);

		                   String Hilfsvariable[] = Hilfsstring.split(";");
		                   int nrFaces=0;
		                   //Check nach Zero Face-Patches
		                   for (int i=0;i<Hilfsvariable.length;i++){
		                	   if(Hilfsvariable[i].contains("nFaces")){String HilfsString2[]=Hilfsvariable[i].split("nFaces");nrFaces=Integer.parseInt(HilfsString2[1].replaceAll("[\\D]", "")); }
		                   }
		                   if(nrFaces>0){		//Wenn realer Patch dann koppeln
		                	globalsb.append(Name);
		       		        globalsb.append("\n");
		                	globalsb.append(" { ");
		                	anzPatch+=1;
			                   for (int i=0;i<Hilfsvariable.length;i++){

			                	   if(Hilfsvariable[i].contains("type")){Hilfsvariable[i]="type mappedWall; sampleMode nearestPatchFace; sampleRegion "+sampleRegion+";"+"samplePatch "+patchName;}

			                   }
			                   //Ansonsten Patch löschen
		                   }else{
		                	   		for (int i=0;i<Hilfsvariable.length;i++){

			                	//   if(Hilfsvariable[i].contains("type")){Hilfsvariable[i]="type empty";}
		                	   			Hilfsvariable[i]="";

			                   }
		                   }


		                   	for (int i=0;i<Hilfsvariable.length;i++){

		                	   globalsb.append(Hilfsvariable[i]);
		                	   if(Hilfsvariable[i]!=""&& i<Hilfsvariable.length-1){ globalsb.append(" ;");}
		                	   globalsb.append("\n");

		                   }
		                   	if(nrFaces>0){
		                   	globalsb.append(" } ");
		                   	}
		      }

		      }

		      br.close();
		     int index=globalsb.indexOf(String.valueOf(anz_gesamt+6));
		     String s=String.valueOf(anz_gesamt);
		     globalsb.delete(index,s.length()+index);
		     globalsb.insert(index,anzPatch+6);

		    	      } catch (Exception e) {System.out.println("Fehler:"+e.toString());}


		 try {


			  File datei = new File(path);
		      FileWriter ausgabestrom = new FileWriter(datei);
		      PrintWriter ausgabe= new PrintWriter(ausgabestrom);


		      ausgabe.println(globalsb.toString());


		      ausgabe.close();
		      ausgabestrom.close();



		}
		catch (Exception e) {System.out.println("Einlesen "+e.toString());}



		 */
	}

	public static void stl_out(TriangleSet ts, String param, PrintWriter ausgabe) throws IOException {

		if (ts != null) {
			ausgabe.println("solid kapillare");

			for (int i = 0; i < ts.getNumber(); i++) {
				final double[][] corners = ts.triangles[i].getCorners();

				ausgabe.println("  facet normal 0.000000e+00 0.000000e+00 0.000000e+00");
				ausgabe.println("    outer loop");

				for (int j = 0; j < 3; j++) {
					ausgabe.println("      vertex  " + (float) corners[0][j] + "  " + (float) corners[1][j] + "  " + (float) corners[2][j]);
				}
				ausgabe.println("    endloop");
				ausgabe.println("  endfacet");

			}
			ausgabe.println("endsolid kapillare");
		}
	}

	public static void stl_to_tri(String inputfile, String outputfile) {
		try {
			final File datei = new File(outputfile);
			final FileWriter ausgabestrom = new FileWriter(datei);
			final PrintWriter ausgabe = new PrintWriter(ausgabestrom);
			final stl_to_tri datawriter = new stl_to_tri();
			stl_to_tri.main(inputfile, ausgabe);
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
	
	public static void stl_to_tri_detector(String inputfile, String outputfile) {
		try {
			final File datei = new File(outputfile);
			final FileWriter ausgabestrom = new FileWriter(datei);
			final PrintWriter ausgabe = new PrintWriter(ausgabestrom);
			final stl_to_tri datawriter = new stl_to_tri();
			stl_to_tri.stl_to_detector(inputfile, ausgabe);
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

	/**
	 * Algorithmus zum invertieren einer Matrix, nach Bronstein "Matrizeninversion/Austauschverfahren" p7.2.1.4
	 *
	 * @matrix
	 */
	public static double[][] inversion(double[][] M) {
		double[][] matrix = new double[M.length][M.length];
		for (int i = 0; i < M.length; i++) {
			for (int j = 0; j < M.length; j++) {
				matrix[i][j] = M[i][j];
			}
		}
		final int N = matrix.length;
		final int[] zeilennummern = new int[N];
		for (int i = 0; i < N; i++) {
			zeilennummern[i] = i;
		}

		for (int zeile = 0; zeile < N; zeile++) {
			// Pivotelementsuche
			int maxindex = zeile;
			double summe = 0;
			for (int i = 0; i < N; i++) {
				summe += Math.abs(matrix[zeile][i]);
			}
			double maximum = Math.abs(matrix[zeile][zeile] / summe);
			for (int lauf = zeile + 1; lauf < N; lauf++) {
				summe = 0;
				for (int i = 0; i < N; i++) {
					summe += Math.abs(matrix[lauf][i]);
				}

				final double abszahl = Math.abs(matrix[lauf][zeile] / summe);
				if (abszahl > maximum) {
					maxindex = lauf;
					maximum = abszahl;
				}
			}
			// Zeile und Reihenfolgearray mit Pivotzeile tauschen
			final int merk = zeilennummern[zeile];
			zeilennummern[zeile] = zeilennummern[maxindex];
			zeilennummern[maxindex] = merk;
			tauschen(matrix[zeile], matrix[maxindex]);
			// Austauschschritt durchfuehren
			final double[][] matrix_ = new double[N][N];
			for (int i = 0; i < N; i++) {
				for (int j = 0; j < N; j++) {
					if ((i == zeile) && (j != zeile)) {
						matrix_[i][j] = -matrix[i][j] / matrix[zeile][zeile];
					}
					if ((i != zeile) && (j == zeile)) {
						matrix_[i][j] = matrix[i][j] / matrix[zeile][zeile];
					}
					if ((i != zeile) && (j != zeile)) {
						matrix_[i][j] = matrix[i][j] - matrix[i][zeile] * matrix[zeile][j] / matrix[zeile][zeile];
					}
				}
			}
			matrix_[zeile][zeile] = 1 / matrix[zeile][zeile];
			matrix = matrix_;
		}
		// Zuruecktauschen der Spalten in die richtige Reihenfolge

		final double[][] matrix_ = new double[N][N];
		for (int i = 0; i < N; i++) {
			for (int j = 0; j < N; j++) {
				matrix_[i][zeilennummern[j]] = matrix[i][j];
			}
		}
		matrix = matrix_;

		// System.out.println("ZNR:\n"+ArrayToString(zeilennummern));
		return matrix;
	}

	/**
	 * Methode tauscht zwei Arrayinhalte
	 */
	public static void tauschen(double[] a, double[] b) {
		double merk;
		for (int i = 0; i < a.length; i++) {
			merk = b[i];
			b[i] = a[i];
			a[i] = merk;
		}
	}

	/**
	 * Methode gibt den kleinsten Abstand zwischen einer Geraden s + x*t und einem Punkt P zurueck
	 */
	public static double distance(double[] s, double[] t, double[] p) {
		final int N = p.length;
		final double[] pms = new double[N];
		for (int i = 0; i < N; i++) {
			pms[i] = p[i] - s[i];
		}// Hilfsvariable erzeugen
		double z = 0, n = 0, k = 0;
		for (int i = 0; i < N; i++) {
			z += pms[i] * t[i];
			n += t[i] * t[i];
		}
		k = z / n;// Weitere Hilfsvariable
		double d2 = 0;
		for (int i = 0; i < N; i++) {
			d2 += (pms[i] - k * t[i]) * (pms[i] - k * t[i]);
		}
		final double d = Math.sqrt(d2);
		return d;
	}

	/**
	 * Methode gibt den Abstand zwischen zwei Punkten zurueck
	 */
	public static double distance(double[] p1, double[] p2) {
		final int N = p1.length;
		double d2 = 0;
		for (int i = 0; i < N; i++) {
			d2 += (p1[i] - p2[i]) * (p1[i] - p2[i]);
		}
		return Math.sqrt(d2);
	}

	/**
	 * Methode gibt den naechsten Punkt auf einer Geraden s + x*t zu einem Punkt P zurueck
	 */
	public static double[] nearest(double[] s, double[] t, double[] p) {
		final double x = nearest_param(s, t, p);
		return point(s, t, x);
	}

	/**
	 * Methode gibt den Punkt auf einer Geraden s + x*t zurueck, welcher zum Parameter x gehoert
	 *
	 * @s Ortsvektor
	 * @t Richtungsvektor
	 * @x Parameter der Geraden
	 */
	public static double[] point(double[] s, double[] t, double x) {
		final int N = s.length;
		final double[] result = new double[N];
		for (int i = 0; i < N; i++) {
			result[i] = s[i] + x * t[i];
		}
		return result;
	}

	/**
	 * Methode gibt den Punkt auf einem Strahl zurueck
	 *
	 * @ray Strahl
	 * @x Parameter des Strahls
	 */
	public static double[] point(Ray ray, double x) {
		final double[] result = new double[3];
		if (ray.kaustik == false) // Gerader Lichtstrahl
		{
			for (int i = 0; i < 3; i++) {
				result[i] = ray.ov[i] + x * ray.rv[i];
			}
		}
		else // Strahl mit Kaustik
		{
			final double k_wz = ray.k * ray.w0 * Math.sqrt(1 + (x - ray.ztt) * (x - ray.ztt) / (ray.z0 * ray.z0));
			for (int i = 0; i < 3; i++) {
				result[i] = ray.ov[i] + x * ray.rv[i] + k_wz * ray.nzs[i];
				// System.out.println(result[2]);
			}
		}
		return result;
	}

	public static double[] direction(Ray ray, double x) {
		final double[] result = new double[3];
		if (ray.kaustik == false) // Gerader Lichtstrahl
		{
			for (int i = 0; i < 3; i++) {
				result[i] = ray.rv[i];
			}
		}
		else // Strahl mit Kaustik
		{
			final double konst = (x - ray.ztt) * ray.k * ray.w0 / (ray.z0 * ray.z0 * Math.sqrt(1 + (x - ray.ztt) * (x - ray.ztt) / (ray.z0 * ray.z0)));
			for (int i = 0; i < 3; i++) {
				result[i] = ray.rv[i] + konst * ray.nzs[i];
			}
		}
		return result;
	}

	public static double[] perp_to_v1_in_v1v2_plane(double[] v1, double[] v2) {
		final double[] v3 = xproduct(v1, v2);
		return xproduct(v3, v1);
	}

	public static double nearest_param(double[] s, double[] t, double[] p) {
		final int N = p.length;
		double t2 = 0;
		double pmst = 0;
		for (int i = 0; i < N; i++) {
			t2 += t[i] * t[i];
			pmst += (p[i] - s[i]) * t[i];
		}
		final double x = pmst / t2;
		return x;
	}

	/**
	 * Methode berechnet das Produkt aus einer quadratischen Matrix und einem Vektor
	 */
	public static double[] product(double[][] m, double[] v) {
		final int N = v.length;
		final double[] res = new double[N];
		for (int zeile = 0; zeile < N; zeile++) {
			for (int spalte = 0; spalte < N; spalte++) {
				res[zeile] += m[zeile][spalte] * v[spalte];
			}
		}
		return res;
	}

	/**
	 * Methode berechnet das Produkt aus zwei quadratischen Matrizen
	 */
	public static double[][] product(double[][] m1, double[][] m2) {
		final int N = m1[0].length;
		final double[][] res = new double[N][N];
		for (int i = 0; i < N; i++) {
			for (int j = 0; j < N; j++) {
				for (int k = 0; k < N; k++) {
					res[i][j] += m1[i][k] * m2[k][j];
				}
			}
		}
		return res;
	}

	/**
	 * Methode berechnet das Kreuzprodukt zweier 3 komponentiger Vektoren
	 */
	public static double[] xproduct(double[] v1, double[] v2) {
		final double[] res = new double[3];
		res[0] = v1[1] * v2[2] - v1[2] * v2[1];
		res[1] = v1[2] * v2[0] - v1[0] * v2[2];
		res[2] = v1[0] * v2[1] - v1[1] * v2[0];
		return res;
	}

	/**
	 * Methode berechnet das Skalarprodukt zweier 3 komponentiger Vektoren
	 */
	public static double product(double[] v1, double[] v2) {
		return (v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2]);
	}

	/**
	 * Methode berechnet die Summe zweier 3 komponentiger Vektoren
	 */
	public static double[] summe(double[] v1, double[] v2) {
		return new double[] { v1[0] + v2[0], v1[1] + v2[1], v1[2] + v2[2] };
	}

	/**
	 * Methode berechnet die Differenz zweier 3 komponentiger Vektoren
	 */
	public static double[] differenz(double[] v1, double[] v2) {
		return new double[] { v1[0] - v2[0], v1[1] - v2[1], v1[2] - v2[2] };
	}

	/**
	 * Methode berechnet das Produkt zwischen Vektor und Skalar
	 */
	public static double[] product(double x, double[] v) {
		return new double[] { x * v[0], x * v[1], x * v[2] };
	}

	/**
	 * Methode berechnet die Laenge des uebergebenen Vektors
	 */
	public static double norm2(double[] v) {
		return Math.sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
	}

	/**
	 * Methode berechnet das Quadrat der Laenge des uebergebenen Vektors
	 */
	public static double norm2_square(double[] v) {
		return (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
	}

	/**
	 * Methode gibt Einheitsmatrix der Dimension n zurueck
	 */
	public static double[][] neutral(int n) {
		final double[][] e = new double[n][n];
		for (int i = 0; i < n; i++) {
			e[i][i] = 1.0;
		}
		return e;
	}

	/**
	 * Methode gibt die transponierte Matrix zurueck
	 */
	public static double[][] transpose(double[][] m) {
		final double[][] erg = new double[m[0].length][m.length];
		for (int i = 0; i < m.length; i++) {
			for (int j = 0; j < m[0].length; j++) {
				erg[j][i] = m[i][j];
			}
		}
		return erg;
	}

	/**
	 * Methode berechnet den zu einer Ebene mit dem Normalenvektor n gespiegelten Richtungsvektor t.
	 */
	public static double[] getMirrorDirection(double[] n, double[] t) {
		final double[] t_ = new double[3];
		final double k = 2 * (n[0] * t[0] + n[1] * t[1] + n[2] * t[2]) / (n[0] * n[0] + n[1] * n[1] + n[2] * n[2]);
		for (int i = 0; i < 3; i++) {
			t_[i] = -k * n[i] + t[i];
		}
		return t_;
	}

	/**
	 * Methode berechnet den zu einer Ebene mit dem Normalenvektor n gespiegelten Richtungsvektor t. Es wird vorausgesetzt, dass die Normalenvektor auf 1 normiert ist
	 */
	public static double[] getMirrorDirection_N_is_normalized(double[] n, double[] t) {
		final double[] t_ = new double[3];
		final double k = 2 * (n[0] * t[0] + n[1] * t[1] + n[2] * t[2]);
		for (int i = 0; i < 3; i++) {
			t_[i] = -k * n[i] + t[i];
		}
		return t_;
	}

	/**
	 * Methode zum uebergebenden Vektor den normierten auf die Laenge x zurueck
	 */
	public static double[] getNormTox(double[] v, double x) {
		double l = 0;
		final double[] erg = new double[v.length];
		for (int i = 0; i < v.length; i++) {
			l += v[i] * v[i];
		}
		l = Math.sqrt(l);
		for (int i = 0; i < v.length; i++) {
			erg[i] = x * v[i] / l;
		}
		return erg;
	}

	/**
	 * Methode gibt den Parameter x zurueck, welcher in eine Gerade sg + x*t eingesetzt werden muss um einen Punkt auf der Ebene se + r*v1 + k*v2 zu erreichen. Also Schnitt: Gerade - Ebene Es werden uebergeben:
	 *
	 * @se Ortsvektor der Ebene
	 * @nv Normalenvektor der Ebene (v1 x v2)
	 * @sg Ortsvektor der Gerade
	 * @tnv Skalarprodukt aus Richtungsvektor der Geraden und Normalenvektor der Ebene
	 */
	public static double paramx_intersection_plane_line(double[] se, double[] nv, double[] sg, double tnv) {
		double z = 0;
		for (int i = 0; i < 3; i++) {
			z += (se[i] - sg[i]) * nv[i];
		}
		return z / tnv;
	}

	/**
	 * Methode gibt die reellen Loesungen einer quadratischen Gleichung zurueck A * x**2 + B * x + C = 0 Das Array hat drei Werte: [Anzahl der Loesungen, Loesung 1, Loesung 2]
	 */
	public static double[] solvePolynom2(double A, double B, double C) {
		final double p = B / A;
		final double q = C / A;
		double f;
		double D;
		final double[] lsg = new double[3];
		if (Math.abs(p) > 1) {
			f = Math.abs(p);
			D = 0.25 - q / p / p;
		}
		else {
			f = 1;
			D = Math.pow(p / 2, 2) - q;
		}
		// Wenn es komplexe Loesungen gibt wuerden die durch -p/2 + i*f*sqrt(-D) und -p/2 - i*f*sqrt(-D) dargestellt
		// Dann ist D<0
		// System.out.println("D = "+D+"   A:"+A+" B:"+B+" C:"+C);
		D = D;
		if (D >= 0) {
			lsg[0] = 2;
			lsg[1] = Math.abs(p / 2) + f * Math.sqrt(D);
			if (p > 0) {
				lsg[1] = -lsg[1];
			}
			if (lsg[1] == 0) {
				lsg[2] = 0;
			}
			else {
				lsg[2] = q / lsg[1];

				// if (lsg[2]-lsg[1]<1e-10){
				// lsg[0]=1;
				// }

			}
		}
		else if (-D < 1e-15) { // Das wurde zugefuegt wegen Diskriminante manchmal nur wegen Numerik minimal kleiner Null anstatt Null
			lsg[0] = 1;
			lsg[1] = -p / 2;
		}
		return lsg;
	}

	/**
	 * Methode wandelt Array in String um
	 */
	public static String ArrayToString(double[][] matrix) {
		String s = "";
		for (int i = 0; i < matrix.length; i++) {
			for (int j = 0; j < matrix[0].length; j++) {
				s = s + matrix[i][j] + "\t";
			}
			s = s + "\n";
		}
		return s;
	}

	/**
	 * Methode wandelt Array in String um
	 */
	public static String ArrayToString(int[] vektor) {
		String s = "";
		for (int i = 0; i < vektor.length; i++) {
			s = s + vektor[i] + "\t";
		}
		s = s + "\n";

		return s;
	}

	/**
	 * Methode wandelt Array in String um
	 */
	public static String ArrayToString(double[] vektor) {
		String s = "";
		for (int i = 0; i < vektor.length; i++) {
			s = s + vektor[i] + "\t";
		}
		s = s + "\n";

		return s;
	}

	/*+
	 * Methode liefert zu einem Vektor einen senkrechten Vektor zurueck
	 */
	public static double[] getNormalVector(double[] v) {
		final double[] s = new double[3];
		if ((v[0] != 0) || (v[1] != 0)) {
			if (v[0] != 0) {
				s[0] = v[1];
				s[1] = -v[0];
			}
			else {
				s[0] = -v[1];
				s[1] = v[0];
			}
		}
		else {
			if ((v[0] != 0) || (v[2] != 0)) {
				if (v[0] != 0) {
					s[0] = v[2];
					s[2] = -v[0];
				}
				else {
					s[0] = -v[2];
					s[2] = v[0];
				}
			}
			else {
				if ((v[1] != 0) || (v[2] != 0)) {
					if (v[1] != 0) {
						s[1] = v[2];
						s[2] = -v[1];
					}
					else {
						s[1] = -v[2];
						s[2] = v[1];
					}
				}
			}
		}
		return s;
	}

	public static double[] paramx_intersection_plane_gaussianray(double[] ov, double[] ne, Ray ray) {

		final double k1 = product(ne, ray.rv);
		final double k2 = ray.k * ray.w0 * product(ne, ray.nzs);
		final double k3 = product(ne, differenz(ov, ray.ov));
		// System.out.println("k1:"+k1+"  k2:"+k2+"  k3:"+k3);
		final double z0_2 = ray.z0 * ray.z0;

		final double[] lsg = new double[3];

		if (Math.abs(((k1 * k1 * z0_2 + k1 * k1 * ray.ztt * ray.ztt - 2 * k1 * k3 * ray.ztt - k2 * k2 + k3 * k3)) / (k1 * k1 * z0_2 - k2 * k2)) < 1e-10) {
			lsg[0] = 1;
			lsg[1] = (k1 * k3 * z0_2 - k2 * k2 * ray.ztt) / (k1 * k1 * z0_2 - k2 * k2);
			lsg[2] = 0;
		}
		else {

			lsg[1] = (k1 * k3 * z0_2 - k2 * k2 * ray.ztt + k2 * ray.z0 * Math.sqrt((k1 * k1 * z0_2 + k1 * k1 * ray.ztt * ray.ztt - 2 * k1 * k3 * ray.ztt - k2 * k2 + k3 * k3))) / (k1 * k1 * z0_2 - k2 * k2);
			lsg[2] = -(k2 * k2 * ray.ztt - k1 * k3 * z0_2 + k2 * ray.z0 * Math.sqrt((k1 * k1 * z0_2 + k1 * k1 * ray.ztt * ray.ztt - 2 * k1 * k3 * ray.ztt - k2 * k2 + k3 * k3))) / (k1 * k1 * z0_2 - k2 * k2);
		}

		final double test1 = Math.abs(k2 * Math.sqrt(1 + Math.pow((lsg[1] - ray.ztt) / ray.z0, 2)) - k3 + lsg[1] * k1);
		final double test2 = Math.abs(k2 * Math.sqrt(1 + Math.pow((lsg[2] - ray.ztt) / ray.z0, 2)) - k3 + lsg[2] * k1);

		if ((test1 < 1e-8) && (test2 < 1e-8)) {
			lsg[0] = 2;
		}
		else if ((test1 < 1e-8) && (test2 > 1e-8)) {
			lsg[0] = 1;
			lsg[2] = 0;
		}
		else if ((test1 > 1e-8) && (test2 < 1e-8)) {
			lsg[0] = 1;
			lsg[1] = lsg[2];
			lsg[2] = 0;
		}
		else {
			lsg[0] = 0;
		}

		return lsg;
	}

	// // Move Triangle

	

	public static double get_current_laser_power(File beam_power_dat, double time) { // FFE, read the power_dat -->"timestep\tpower" and returns the power for the timestep closest to time
		String line;
		String[] line_split;
		final List<Double> times = new ArrayList<Double>();
		final List<Double> powers = new ArrayList<Double>();
		try {
			final BufferedReader br = new BufferedReader(new FileReader(beam_power_dat));

			while (((line = br.readLine()) != null) && (!line.trim().equals("END"))) { // Für jede line
				line = line.trim();
				line_split = line.split("\t");
				times.add(Double.valueOf(line_split[0]));
				powers.add(Double.valueOf(line_split[1]));
			}
		}
		catch (final IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		int i = 0;
		while (times.get(i) < time) {
			i++;
		}

		return powers.get(i);
	}

	public static String get_current_laser_position(File beam_osc_dat, double time) { // FFE, read the power_dat -->"timestep\tx, y, z" in µm and returns the power for the timestep closest to time
		String pos;
		String line;
		final List<Double> times = new ArrayList<Double>();
		final List<String> positions = new ArrayList<String>();
		try {
			final BufferedReader br = new BufferedReader(new FileReader(beam_osc_dat));

			while (((line = br.readLine()) != null) && (!line.trim().equals("END"))) { // Für jede line
				line = line.trim();
				pos = line.split("\t")[1];

				positions.add(pos);
			}
		}
		catch (final IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		int i = 0;
		while (times.get(i) < time) {
			i++;
		}

		return positions.get(i);
	}

	public static void set_current_laser_power(File raytracer_parameter_file, double power) { // FFE, sets the current laser power in the Parameter_Raytracer file - ungetested!!
		String complete_outfile = "";
		try (BufferedReader br = new BufferedReader(new FileReader(raytracer_parameter_file))) { // ASCII STL öffnen
			String line;
			String line_out = "";
			while (((line = br.readLine()) != null) && (!line.trim().equals("END"))) { // Für jede line
				line = line.trim();
				if (line.startsWith("beamPower")) {
					// System.out.println(line);

					line_out = line_out + "beamPower  = " + String.valueOf(power) + "\n";
				}
				else {
					line_out = line + "\n";
				}
				complete_outfile = complete_outfile + line_out;
			}
		}
		catch (final FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		catch (final IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		try (BufferedWriter br = new BufferedWriter(new FileWriter(raytracer_parameter_file))) {
			br.write(complete_outfile);
		}
		catch (final IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}

	public static void set_current_laser_position(File raytracer_parameter_file, String position) { // FFE, sets the current laser power in the Parameter_Raytracer file - ungetested!!
		String complete_outfile = "";
		try (BufferedReader br = new BufferedReader(new FileReader(raytracer_parameter_file))) { // ASCII STL öffnen
			String line;
			String line_out = "";
			while (((line = br.readLine()) != null) && (!line.trim().equals("END"))) { // Für jede line
				line = line.trim();
				if (line.startsWith("Incoming_beam_position_vector")) {
					// System.out.println(line);

					line_out = line_out + "Incoming_beam_position_vector  = " + position + "\n";
				}
				else {
					line_out = line + "\n";
				}
				complete_outfile = complete_outfile + line_out;
			}
		}
		catch (final FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		catch (final IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		try (BufferedWriter br = new BufferedWriter(new FileWriter(raytracer_parameter_file))) {
			br.write(complete_outfile);
		}
		catch (final IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}

	public static void make_dir(String name) {

		final File theDir = new File(name);
		// if the directory does not exist, create it
		if (!theDir.exists()) {
			System.out.println("creating directory: " + name);
			boolean result = false;

			try {
				theDir.mkdir();
				result = true;
			}
			catch (final SecurityException se) {
				// handle it
			}
			if (result) {
				System.out.println("DIR created");
			}
		}
	}

	public static void delete_file(File todelete) {
		try {
			System.out.println("Deleting: " + todelete.getName());
			Files.delete(Paths.get(todelete.getAbsolutePath()));
		}
		catch (final IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}

	public static double isNumeric(String str) // is that String parsable to a number ?
	{
		try {
			final double d = Double.parseDouble(str);
			return d;
		}
		catch (final NumberFormatException nfe) {
			return -1;
		}
	}

	public static double[] get_cap_dimensions(File cap_file) {
		/*
		 * liefert die maximale Ausdehnung einer stl-geometrie zurück
		 */
		double max_x = -9999999;
		double max_y = -9999999;
		double max_z = -9999999;
		double min_x = 9999999;
		double min_y = 9999999;
		double min_z = 9999999;
		final double[] cap_dims = { 0, 0, 0, 0, 0, 0 };
		try (BufferedReader br = new BufferedReader(new FileReader(cap_file))) { // ASCII STL öffnen
			String line;
			while (((line = br.readLine()) != null) && (!line.trim().equals("END"))) { // Für jede line
				line = line.trim();
				if (line.startsWith("vertex")) {
					line = line.split("vertex")[1];
					line = line.trim();
					//for (int jj = 0; jj < line.split(" ").length; jj++) {
					//	System.out.println(jj + " " + line.split(" ")[jj] + " ");
					//}
					//System.out.println();
					// System.out.println(line);
					if (Double.parseDouble(line.split(" ")[0]) > max_x) {
						max_x = Double.parseDouble(line.split(" ")[0]);
					}
					if (Double.parseDouble(line.split(" ")[1]) > max_y) {
						max_y = Double.parseDouble(line.split(" ")[1]);
					}
					if (Double.parseDouble(line.split(" ")[2]) > max_z) {
						max_z = Double.parseDouble(line.split(" ")[2]);
					}
					if (Double.parseDouble(line.split(" ")[0]) < min_x) {
						min_x = Double.parseDouble(line.split(" ")[0]);
					}
					if (Double.parseDouble(line.split(" ")[1]) < min_y) {
						min_y = Double.parseDouble(line.split(" ")[1]);
					}
					if (Double.parseDouble(line.split(" ")[2]) < min_z) {
						min_z = Double.parseDouble(line.split(" ")[2]);
					}
				}
			}
		}
		catch (final FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		catch (final IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		cap_dims[0] = max_x;
		cap_dims[1] = max_y;
		cap_dims[2] = max_z;
		cap_dims[3] = min_x;
		cap_dims[4] = min_y;
		cap_dims[5] = min_z;
		return cap_dims;
	}

	public static void generate_block_z(File blockmeshFile, double new_z) {
		String complete_out = "";
		boolean no_change = true;
		boolean in_vertices = false;
		int vert_counter = 0;
		try (BufferedReader br = new BufferedReader(new FileReader(blockmeshFile))) { // ASCII STL öffnen
			String line;
			while (((line = br.readLine()) != null) && (!line.trim().equals("END"))) {
				if (line.contains("vertices")) { // hier beginnt die block definition
					in_vertices = true;
					no_change = true;
					// System.out.println(line);
				}
				if (in_vertices == true) {
					no_change = true;
					vert_counter++; // Zeilen zählen
					// System.out.println(line + "  " + vert_counter);
				}
				if ((vert_counter > 6) && (vert_counter < 11)) { // this are the lines that are changed
					no_change = false;
					String firstline;
					// System.out.println(line);
					firstline = line.split("\\)")[0];
					for (int jj = 0; jj < firstline.split(" ").length - 1; jj++) {
						complete_out = complete_out + firstline.split(" ")[jj] + " ";
					}
					complete_out = complete_out + String.valueOf(new_z) + ")\n ";
					System.out.println();
					// complete_out = complete_out + line.split(" ")[0] + " " + line.split(" ")[1] + " " + String.valueOf(new_z) + ")\n";
					// System.out.println("edit " + line.split(" ")[0] + " " + line.split(" ")[1] + " " + String.valueOf(new_z) + ")\n");
				}
				if (vert_counter == 11) {
					in_vertices = false;
					no_change = true;
					// System.out.println("leaving block");
				}

				if (no_change) {
					complete_out = complete_out + line + "\n";
				}
			}
			br.close();
		}
		catch (final FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		catch (final IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		try (BufferedWriter br = new BufferedWriter(new FileWriter(blockmeshFile))) {
			br.write(complete_out);
		}
		catch (final IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}

	public static void write_logfile(File f, String new_line) {
		String complete_out = "";

		if (f.exists() && !f.isDirectory()) {

			try {
				final FileReader FR = new FileReader(f);
				final BufferedReader br = new BufferedReader(FR); // ASCII STL öffnen
				String line;
				while (((line = br.readLine()) != null) && (!line.trim().equals("END"))) {
					complete_out = complete_out + line + "\n";
				}
				br.close();
				System.out.println("logfile appended");
			}
			catch (final FileNotFoundException e) {
				System.out.println("logfile generated");
				e.printStackTrace();
			}
			catch (final IOException e) {
				e.printStackTrace();
			} //TODO Header für Logfile
	  } else {
		  complete_out="I.run"+"\t"+"\t"+"II.Cap_depth"+"\t"+"III.Cap_Surface"+"\t"+"IV.Infeed"+"\t"+"V.Scaling"+"\t"+"VI.Percentage"+"\n";
	  }
		try(BufferedWriter bw = new BufferedWriter(new FileWriter(f))){
			bw.write(complete_out+new_line+"\n");
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
	
	public static void change_melt(File runDir) {
		String t_melt;
		t_melt = Algorithmen.get_OF_parameter(new File (runDir, "constant/transportProperties"), "T_liquid");
		//Schmelzbadbild und csv datei
				
		PrintWriter writer = null;
		try {
			writer = new PrintWriter(new BufferedWriter(
					new FileWriter("D:\\OF_models\\trunk\\OpenFoam_Fetzer\\Bin\\FEM\\save_melt_data.temp")));
		} catch (IOException e1) {
			// TODO Auto-generated catch block
			e1.printStackTrace();
		}

		BufferedReader br = null;
		FileReader reader = null;
		try {
			reader = new FileReader("D:\\OF_models\\trunk\\OpenFoam_Fetzer\\Bin\\FEM\\save_melt_data.py");
			br = new BufferedReader(reader);
			String line;

			while ((line = br.readLine()) != null) {
				if (line.startsWith("Eval_Model = LegacyVTKReader")) {
					line = "Eval_Model = LegacyVTKReader( FileNames=['D:\\\\OF_models\\\\trunk\\\\OpenFoam_Fetzer\\\\Runs\\\\"
							+ runDir.getName() + "\\\\VTK\\\\" + runDir.getName() + "_4.vtk'] )";
				}
				if (line.startsWith("WriteImage")) {
					line = "WriteImage('D:\\\\OF_models\\\\trunk\\\\OpenFoam_Fetzer\\\\Runs\\\\" + runDir.getName()
							+ "\\\\melt.png')";
				}
				if (line.startsWith("writer = CreateWriter")) {
					line = "writer = CreateWriter(r\"D:\\OF_models\\trunk\\OpenFoam_Fetzer\\Runs\\"+runDir.getName()+"\\melt.csv\")";
				}
				if (line.startsWith("Contour4.Isosurfaces")) {
					line = "Contour4.Isosurfaces = ["+t_melt+"]";
				}
				writer.println(line);
			}
			writer.close();

		} catch (IOException e) {
			e.printStackTrace();
		} finally {
			try {
				br.close();
				reader.close();
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}
		File realName = new File("D:\\OF_models\\trunk\\OpenFoam_Fetzer\\Bin\\FEM\\save_melt_data.py");
		File f = new File("D:\\OF_models\\trunk\\OpenFoam_Fetzer\\Bin\\FEM\\save_melt_data.temp");
		Algorithmen.copy_file(f, realName);
		f.delete();
		}


	
	public static void change_scenerytext(File runDir) {

		
		// Szenenbild
		PrintWriter writer = null;
		try {
			writer = new PrintWriter(new BufferedWriter(
					new FileWriter("FEM\\Test_Scenery.temp")));
		} catch (IOException e1) {
			// TODO Auto-generated catch block
			e1.printStackTrace();
		}

		BufferedReader br = null;
		FileReader reader = null;
		try {
			reader = new FileReader("FEM\\Test_Scenery.py");
			br = new BufferedReader(reader);
			String line;

			while ((line = br.readLine()) != null) {
				if (line.startsWith("Eval_Model = LegacyVTKReader")) {
					line = "Eval_Model = LegacyVTKReader( FileNames=['D:\\\\OF_models\\\\trunk\\\\OpenFoam_Fetzer\\\\Runs\\\\"
							+ runDir.getName() + "\\\\VTK\\\\" + runDir.getName() + "_4.vtk'] )";
				}
				if (line.startsWith("WriteImage")) {
					line = "WriteImage('D:\\\\OF_models\\\\trunk\\\\OpenFoam_Fetzer\\\\Runs\\\\" + runDir.getName()
							+ "\\\\scene.png')";
				}
				writer.println(line);
			}
			writer.close();

		} catch (IOException e) {
			e.printStackTrace();
		} finally {
			try {
				br.close();
				reader.close();
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}
		File realName = new File("FEM\\Test_Scenery.py");
		File f = new File("FEM\\Test_Scenery.temp");
		Algorithmen.copy_file(f, realName);
		f.delete();
	
}
//Marcel 24.11.2015  
	  public static double[] get_csv_dimensions(File csv_file){ 
			double max_p0 = -9999999;
			double max_p1 = -9999999;
			double max_p2 = -9999999;
			double min_p0 = 9999999;
			double min_p1 = 9999999;
			double min_p2 = 9999999;
			double[] csv_dims = {0,0,0,0,0,0};
			try(BufferedReader br = new BufferedReader(new FileReader(csv_file))){ // ASCII STL öffnen
			      String line;	      
			      while ((line = br.readLine())!=null){  // Für jede line
			    	  if ((line.startsWith("\"T\"")) == false){

			    		  if (Double.parseDouble(line.split(",")[7]) > max_p0 ){
			    			  max_p0 = Double.parseDouble(line.split(",")[7]);
			    		  }
			    		  if (Double.parseDouble(line.split(",")[8]) > max_p1 ){
			    			  max_p1 = Double.parseDouble(line.split(",")[8]);
			    		  }
			    		  if (Double.parseDouble(line.split(",")[9]) > max_p2 ){
			    			  max_p2 = Double.parseDouble(line.split(",")[9]);
			    		  }
			    		  if (Double.parseDouble(line.split(",")[7]) < min_p0 ){
			    			  min_p0 = Double.parseDouble(line.split(",")[7]);
			    		  }
			    		  if (Double.parseDouble(line.split(",")[8]) < min_p1 ){
			    			  min_p1 = Double.parseDouble(line.split(",")[8]);
			    		  }
			    		  if (Double.parseDouble(line.split(",")[9]) < min_p2 ){
			    			  min_p2 = Double.parseDouble(line.split(",")[9]);
			    		  }
			    	  }
			      }
			      } catch (FileNotFoundException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				} catch (IOException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
			csv_dims[0] = max_p0; csv_dims[1]=max_p1; csv_dims[2] = max_p2; csv_dims[3]=min_p0; csv_dims[4] = min_p1; csv_dims[5]=min_p2;
			return csv_dims;
		}
		
 public static void controldict_set_tend_dt(File cdict,String end_time, String dt){
	 /*
	 String complete_out = "";
	  boolean changed = false;
	  if(cdict.exists() && !cdict.isDirectory()) { 
			try {
				final FileReader FR = new FileReader(cdict);
				final BufferedReader br = new BufferedReader(FR); // ASCII STL öffnen
				String line;
				while (((line = br.readLine()) != null) && (!line.trim().equals("END"))) {
					changed = false;
					if (line.startsWith("endTime")) {
						complete_out = complete_out + "endTime    " + end_time + ";\n";
						changed = true;
					}
					if (line.startsWith("deltaT")) {
						complete_out = complete_out + "deltaT    " + dt + ";\n";
						changed = true;
					}
					if (line.startsWith("writeControl")) {
						complete_out = complete_out + "writeControl    adjustableRunTime;\n";
						changed = true;
						System.out.println("write Interval changed");
					}
					if (line.startsWith("writeInterval")) {
						//complete_out = complete_out + "writeInterval    " + dt + ";\n";
						//complete_out = complete_out + "writeInterval    adjustableRunTime;\n";
						//changed = true;
						//System.out.println("write Interval changed");
					}
					if (changed == false) {
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
			try (BufferedWriter bw = new BufferedWriter(new FileWriter(cdict))) {
				bw.write(complete_out);
			}
			catch (final IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}
		*/

	}
 public static void set_ClosedLoop_rundirectory(File cdict,String rundir){
	 String complete_out = "";
	 boolean changed = false;
	  if(cdict.exists() && !cdict.isDirectory()) { 
			try {
				final FileReader FR = new FileReader(cdict);
				final BufferedReader br = new BufferedReader(FR); // ASCII STL öffnen
				String line;
				while (((line = br.readLine()) != null) && (!line.trim().equals("END"))) {
					changed = false;
					if (line.startsWith("caseDir")) {
						complete_out = complete_out + "caseDir    " + rundir + ";\n";
						changed = true;
					}
					if (changed == false) {
						complete_out = complete_out + line + "\n";
					}
				}
				br.close();
				System.out.println("ClosedLoop caseDir changed");
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
			try (BufferedWriter bw = new BufferedWriter(new FileWriter(cdict))) {
				bw.write(complete_out);
			}
			catch (final IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}

	}
 
 public static List<Double[]> read_colormap(String colormap_file) {

		List<Double[]> rgbs = new ArrayList<Double[]>();
		try (BufferedReader br = new BufferedReader(new FileReader(colormap_file))) { // ASCII STL öffnen
			String line;
			int counter = 0;
			while (((line = br.readLine()) != null) && (!line.trim().equals("END"))) { // Für jede line
				Double[] colors = new Double[3];
				line = line.trim().substring(1, line.length()-2);
				colors[0] = Double.valueOf(line.split(", ")[0]);
				colors[1] = Double.valueOf(line.split(", ")[1]);
				colors[2] = Double.valueOf(line.split(", ")[2]);
				//System.out.println(colors[2]);
				rgbs.add(colors);
			}
			counter++;
		}
		catch (final FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		catch (final IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		return rgbs;
}
 
 public static void triSet_to_stl(TriangleSet triSet, File datei){
	 
	 	
		FileWriter ausgabestrom;
		try {
			ausgabestrom = new FileWriter(datei);
		
		final PrintWriter ausgabe = new PrintWriter(ausgabestrom);

		ausgabe.println("solid corrected kapillare");

		for (int k = 0; k < triSet.getNumber(); k++) {
			ausgabe.println("facet normal " + triSet.triangles[k].nv[0] + " " + triSet.triangles[k].nv[1] + " " + triSet.triangles[k].nv[2]);
			ausgabe.println(" outer loop");
			ausgabe.println("  vertex " +triSet.triangles[k].p1[0] + " " + triSet.triangles[k].p1[1] + " " +triSet.triangles[k].p1[2]);
			ausgabe.println("  vertex " +triSet.triangles[k].p2[0] + " " + triSet.triangles[k].p2[1] + " " +triSet.triangles[k].p2[2]);
			ausgabe.println("  vertex " +triSet.triangles[k].p3[0] + " " + triSet.triangles[k].p3[1] + " " +triSet.triangles[k].p3[2]);
			ausgabe.println(" endloop");
			ausgabe.println("endfacet");
		}
		ausgabe.println("endsolid corrected kapillare");

		ausgabestrom.close();
		}
		catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
 }
 public static void append_doubles_to_file(File f, double[] d, double i){
	    String s="";
	    s = s+d[0] + ", " + d[1] + ", " + d[2] + ", "+i+"\n";
	    try {
	      Files.write(Paths.get(f.getAbsolutePath()), s.getBytes(), StandardOpenOption.APPEND);
	  }catch (IOException e) {
	      //exception handling left as an exercise for the reader
	  }
	  }
 
  public static String get_OF_last_writeStep(String OFrunDir){
	  String lastTime = "";
	  
	  File f = new File(OFrunDir); // current directory
	  List<Float> time_dir_list_f = new ArrayList<>();
	  List<String> time_dir_list_s = new ArrayList<>();
	  File[] files = f.listFiles();
	  for (int i = 0; i<files.length;i++){ // get all directories, whos name is a number
	    	if ( files[i].isDirectory() && NumberUtils.isNumber(files[i].getName()) ){
	    		System.out.print(files[i].getName() + "  ");
	    		System.out.println( Float.parseFloat(files[i].getName()));
	    		time_dir_list_f.add(Float.parseFloat(files[i].getName()));
	    		time_dir_list_s.add(files[i].getName());
	    	}
	  }
	  int max_i = 0; // find maximum time
	  float max = -999;
	  for (int i = 0; i < time_dir_list_f.size(); i++){
		  if (time_dir_list_f.get(i)>max) {
			  max = time_dir_list_f.get(i);
			  max_i = i;
		  }
	  }
	  System.out.println("Max time = " + time_dir_list_s.get(max_i));
	  return time_dir_list_s.get(max_i);
}
  
  public static void appendfiles(File out, File first, File second){
	  OutputStream out2;
	try {
		out2 = new FileOutputStream(out);
	    byte[] buf = new byte[12];
	        InputStream in = new FileInputStream(first);
	        int b = 0;
	        while ( (b = in.read(buf)) >= 0) {
	            out2.write(buf, 0, b);
	            out2.flush();
	        }
	        in.close();
	        in = new FileInputStream(second);
	        b = 0;
	        while ( (b = in.read(buf)) >= 0) {
	            out2.write(buf, 0, b);
	            out2.flush();
	        }
	        out2.close();
	}
	catch (IOException e) {
		// TODO Auto-generated catch block
		e.printStackTrace();
	}
	}
public static void write_beam_tracks_to_file(File out, List<double[]> beam_tracks){
	try{
		FileWriter fw = new FileWriter(out,true);
	    for (int i = 0;i<beam_tracks.size();i++){
	    	fw.write(String.valueOf(beam_tracks.get(i)[0]) +", " + String.valueOf(beam_tracks.get(i)[1]) + ", " + String.valueOf(beam_tracks.get(i)[2]) + ", " + String.valueOf(beam_tracks.get(i)[3]) + ", ");
	    }
	    fw.write("\n");
	    fw.close();
	} catch (IOException e) {
	   // do something
	}
}

public static void set_Raytracer_Param(File raytracer_parameter_file, String newline) { // FFE, sets the respective line to newline
	String complete_outfile = "";
	try (BufferedReader br = new BufferedReader(new FileReader(raytracer_parameter_file))) { // ASCII STL öffnen
		String line;
		while (((line = br.readLine()) != null) && (!line.trim().equals("END"))) { // Für jede line
			line = line.trim();
			String p_old = line.split("=")[0].trim();
			String p_new = newline.split("=")[0].trim();
			if (p_old.equals(p_new)) {
				line = newline;
			}
			complete_outfile = complete_outfile + line + "\n";
		}
	}
	catch (final FileNotFoundException e) {
		// TODO Auto-generated catch block
		e.printStackTrace();
	}
	catch (final IOException e) {
		// TODO Auto-generated catch block
		e.printStackTrace();
	}
	try (BufferedWriter br = new BufferedWriter(new FileWriter(raytracer_parameter_file))) {
		br.write(complete_outfile);
	}
	catch (final IOException e) {
		// TODO Auto-generated catch block
		e.printStackTrace();
	}
}
 
}