import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintWriter;
import java.io.UnsupportedEncodingException;

//Objekte dieser Klasse stellen einen Satz zusammengehörender Dreiecke dar
public class TriangleSet {
	Triangle[] triangles;
	String name;
	double[] bounds; // xmin, xmax, ymin, ymax, zmin, zmax
	int type;

	public TriangleSet(Triangle[] t, String s) {
		triangles = t;
		name = s;
		bounds = getMinMaxTriangleParameters();
	}

	public TriangleSet(Objectlist ol, String s) {
		final Object[] o = ol.getArray();
		triangles = new Triangle[o.length];
		for (int i = 0; i < o.length; i++) {
			triangles[i] = (Triangle) (o[i]);
		}
		name = s;
		bounds = getMinMaxTriangleParameters();

	}

	public int getNumber() {
		return triangles.length;
	}

	public void setType(int type) {
		this.type = type;
		for (int i = 0; i < triangles.length; i++) {
			triangles[i].type = type;
		}
	}

	public double getPower() {
		double p = 0;
		for (int i = 0; i < triangles.length; i++) {
			p += triangles[i].P;
		}
		return p;
	}

	public double maxTemperature() {
		double max = 0;

		for (int i = 0; i < triangles.length; i++) {
			if (triangles[i].getTemperature() > max) {
				max = triangles[i].getTemperature();
			}
		}

		return max;
	}

	public double getSurfaceArea() {
		double area = 0;
		for (int i = 0; i < triangles.length; i++) {
			area = area + (0.5 * Algorithmen.norm2(Algorithmen.xproduct(triangles[i].v1, triangles[i].v2)));
		}
		return area;
	}

	public double[] getMinMaxTriangleParameters() {
		double[] params = new double[8];
		if (triangles.length > 0) {
			double minL = Algorithmen.norm2(Algorithmen.differenz(triangles[0].p1, triangles[0].p2));
			double maxL = Algorithmen.norm2(Algorithmen.differenz(triangles[0].p1, triangles[0].p2));
			double minx = triangles[0].p1[0];
			double maxx = triangles[0].p1[0];
			double miny = triangles[0].p1[1];
			double maxy = triangles[0].p1[1];
			double minz = triangles[0].p1[2];
			double maxz = triangles[0].p1[2];

			for (int i = 1; i < triangles.length; i++) {
				final double[] p1 = triangles[i].p1;
				minx = p1[0] < minx ? p1[0] : minx;
				maxx = p1[0] > maxx ? p1[0] : maxx;
				miny = p1[1] < miny ? p1[1] : miny;
				maxy = p1[1] > maxy ? p1[1] : maxy;
				minz = p1[2] < minz ? p1[2] : minz;
				maxz = p1[2] > maxz ? p1[2] : maxz;
				final double[] p2 = triangles[i].p2;
				minx = p2[0] < minx ? p2[0] : minx;
				maxx = p2[0] > maxx ? p2[0] : maxx;
				miny = p2[1] < miny ? p2[1] : miny;
				maxy = p2[1] > maxy ? p2[1] : maxy;
				minz = p2[2] < minz ? p2[2] : minz;
				maxz = p2[2] > maxz ? p2[2] : maxz;
				final double[] p3 = triangles[i].p3;
				minx = p3[0] < minx ? p3[0] : minx;
				maxx = p3[0] > maxx ? p3[0] : maxx;
				miny = p3[1] < miny ? p3[1] : miny;
				maxy = p3[1] > maxy ? p3[1] : maxy;
				minz = p3[2] < minz ? p3[2] : minz;
				maxz = p3[2] > maxz ? p3[2] : maxz;
				final double l1 = Algorithmen.norm2(Algorithmen.differenz(p1, p2));
				final double l2 = Algorithmen.norm2(Algorithmen.differenz(p1, p3));
				final double l3 = Algorithmen.norm2(Algorithmen.differenz(p2, p3));
				minL = l1 < minL ? l1 : minL;
				minL = l2 < minL ? l2 : minL;
				minL = l3 < minL ? l3 : minL;
				maxL = l1 > maxL ? l1 : maxL;
				maxL = l2 > maxL ? l2 : maxL;
				maxL = l3 > maxL ? l3 : maxL;
			}
			params = new double[] { minL, maxL, minx, maxx, miny, maxy, minz, maxz };
		}
		return params;
	}
	
	
	public void write_tri_logfile(File f){
		/*
		 * prints everything into one ASCII file
		 */
		try{
			PrintWriter writer = new PrintWriter(f, "UTF-8");
			writer.println("Nr., spx, spy, spz, nvx, nvy, nvz, area, bend_radius, s_tension, dist_interface, inc_power, OF_power, delta_power, temperature, OFpressure, RecoilPressure, delta_pressure, verschiebung, material");
			for (int k = 0; k < triangles.length; k++){
				double A = 0.5 * Algorithmen.norm2(Algorithmen.xproduct(triangles[k].v1, triangles[k].v2)) / 1000000D;
				writer.print(k+ ", " + triangles[k].Schwerpunkt()[0] +", "+ triangles[k].Schwerpunkt()[1]+", "+ triangles[k].Schwerpunkt()[2]+ ", " + triangles[k].nv[0] + ", "+ triangles[k].nv[1] + ", "+ triangles[k].nv[2] + ", " + A);
				writer.print(", "+ triangles[k].Kruemmungsradius +", "+triangles[k].getSurfaceTension()+ ", " +triangles[k].DistanceInterface);
				writer.print(", "+ triangles[k].getPower() + ", " + triangles[k].Power_OF + ", "+ triangles[k].getDeltaP()+ ", "+ triangles[k].Temperature);
				writer.println(", " + triangles[k].OFPressure+", "+ triangles[k].RecoilPressure+ ", " + triangles[k].DeltaPressure+", " + triangles[k].Verschiebung+", " + triangles[k].material);
			}
			writer.close();
		}
		catch(IOException e){
			e.printStackTrace();
		}
	}
}
