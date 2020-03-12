import java.io.BufferedReader;
import java.io.FileReader;

/**
 * Objekte dieser Klasse beinhalten Parameter für eine Raytracingrechnung
 */

public class Parameters {
	String filename;
	int number_of_reflections = 30;
	int number_of_beams = 1000;
	String polarization = "linear_x";
	int number_of_sublevels = 18;
	double[] incoming_beam_position_vector = new double[] { 0, 0, 0 };
	double[] incoming_beam_direction_vector = new double[] { 0, 0, -1 };
	double beam_waist_radius = 1;
	double beam_wavelength = 1;
	double beam_M2 = 1;
	double beam_distance_to_focus = 0;
	double beam_source_radius = 1;
	double beam_source_width_x = 1;
	double beam_source_width_y = 1;
	String beam_source_type = "circle";
	double beam_min_distance = 1e-9;
	String beam_intensity_distribution = "gauss";
	String material_name = "Steel@1030nm 300K";
	double material_IOR_real = 2.59;
	double material_IOR_imag = 4.87;
	// 2. Material FFE
	String material_name2 = "Alul@1030nm 300K";
	double material_IOR_real2 = 3;
	double material_IOR_imag2 = 5;
	double rel_diff_scat;
	//
	boolean target_csv = true;
	boolean target_ply_lin = true;
	boolean target_ply_log = true;
	double bounding_box_spacing = 0.1;
	int number_of_threads = 1;
	public static double heatConductivity;
	public double beamPower = 3000.;
	public boolean multimaterial = false;
	public boolean enable_temperature_dependent_nk = false;

	public static double get_global_heatConductivity() {

		if (heatConductivity == 0 || Double.isNaN(heatConductivity)) { // Falls Parameter noch nicht initialisiert

			final Parameters hilfspara = new Parameters("Parameter_Raytracer.txt");
			heatConductivity = Parameters.heatConductivity;

		}

		return heatConductivity;

	}

	public Parameters(String filename) {
		this.filename = filename;
		try {
			final BufferedReader br = new BufferedReader(new FileReader(filename));
			String line;
			line = br.readLine();
			while ((line = br.readLine()) != null) {
				line = line.toLowerCase();
				line = line.trim();
				if (line.equals("")) {
					continue;
				}
				final char firstsymbol = line.charAt(0);
				if (firstsymbol == '#') {
					continue;
				}
				final String[] params = line.split("=");
				params[0] = params[0].trim();
				params[1] = params[1].trim();
				boolean validparam = false;
				if (params[0].equals("number_of_reflections")) {
					number_of_reflections = Integer.valueOf(params[1]).intValue();
					validparam = true;
				}
				if (params[0].equals("number_of_beams")) {
					number_of_beams = Integer.valueOf(params[1]).intValue();
					validparam = true;
				}
				if (params[0].equals("material_name")) {
					material_name = params[1];
					validparam = true;
				}
				if (params[0].equals("beam_intensity_distribution")) {
					beam_intensity_distribution = params[1];
					validparam = true;
				}
				if (params[0].equals("target_csv")) {
					target_csv = Boolean.valueOf(params[1]).booleanValue();
					validparam = true;
				}
				if (params[0].equals("rel_diff_scat")) {
					rel_diff_scat = Double.valueOf(params[1]).doubleValue();
					validparam = true;
				}
				if (params[0].equals("target_ply_lin")) {
					target_ply_lin = Boolean.valueOf(params[1]).booleanValue();
					validparam = true;
				}
				if (params[0].equals("target_ply_log")) {
					target_ply_log = Boolean.valueOf(params[1]).booleanValue();
					validparam = true;
				}
				if (params[0].equals("material_ior_real")) {
					material_IOR_real = Double.valueOf(params[1]).doubleValue();
					validparam = true;
				}
				if (params[0].equals("material_ior_imag")) {
					material_IOR_imag = Double.valueOf(params[1]).doubleValue();
					validparam = true;
				}
				// Material 2 FFE
				if (params[0].equals("material_ior_real2")) {
					material_IOR_real2 = Double.valueOf(params[1]).doubleValue();
					validparam = true;
				}
				if (params[0].equals("material_ior_imag2")) {
					material_IOR_imag2 = Double.valueOf(params[1]).doubleValue();
					validparam = true;
				}
				if (params[0].equals("material_name2")) {
					material_name2 = params[1];
					validparam = true;
				}
				//
				if (params[0].equals("bounding_box_spacing")) {
					bounding_box_spacing = Double.valueOf(params[1]).doubleValue();
					validparam = true;
				}
				if (params[0].equals("polarization")) {
					polarization = params[1];
					validparam = true;
				}
				if (params[0].equals("number_of_sublevels")) {
					number_of_sublevels = Integer.valueOf(params[1]).intValue();
					validparam = true;
				}
				if (params[0].equals("number_of_threads")) {
					number_of_threads = Integer.valueOf(params[1]).intValue();
					validparam = true;
				}
				if (params[0].equals("beam_waist_radius")) {
					beam_waist_radius = Double.valueOf(params[1]).doubleValue();
					validparam = true;
				}
				if (params[0].equals("beam_wavelength")) {
					beam_wavelength = Double.valueOf(params[1]).doubleValue();
					validparam = true;
				}
				if (params[0].equals("beam_min_distance")) {
					beam_min_distance = Double.valueOf(params[1]).doubleValue();
					validparam = true;
				}
				if (params[0].equals("beam_m2")) {
					beam_M2 = Double.valueOf(params[1]).doubleValue();
					if (beam_M2 < 1) {
						System.out.println("Warning: M-square smaller 1 senseless.");
					}
					validparam = true;
				}
				if (params[0].equals("beam_distance_to_focus")) {
					beam_distance_to_focus = Double.valueOf(params[1]).doubleValue();
					validparam = true;
				}
				if (params[0].equals("beam_source_radius")) {
					beam_source_radius = Double.valueOf(params[1]).doubleValue();
					validparam = true;
				}
				if (params[0].equals("beam_source_width_x")) {
					beam_source_width_x = Double.valueOf(params[1]).doubleValue();
					validparam = true;
				}
				if (params[0].equals("beam_source_width_y")) {
					beam_source_width_y = Double.valueOf(params[1]).doubleValue();
					validparam = true;
				}
				if (params[0].equals("beam_source_type")) {
					beam_source_type = params[1];
					validparam = true;
				}

				if (params[0].equals("heat_conductivity")) {
					heatConductivity = Double.valueOf(params[1]);
					validparam = true;
				}

				if (params[0].equals("incoming_beam_position_vector")) {
					final String[] xyz_p = params[1].split(",");
					incoming_beam_position_vector[0] = Double.valueOf(xyz_p[0]).doubleValue();
					incoming_beam_position_vector[1] = Double.valueOf(xyz_p[1]).doubleValue();
					incoming_beam_position_vector[2] = Double.valueOf(xyz_p[2]).doubleValue();
					validparam = true;
				}
				if (params[0].equals("incoming_beam_direction_vector")) {
					final String[] xyz_p = params[1].split(",");
					incoming_beam_direction_vector[0] = Double.valueOf(xyz_p[0]).doubleValue();
					incoming_beam_direction_vector[1] = Double.valueOf(xyz_p[1]).doubleValue();
					incoming_beam_direction_vector[2] = Double.valueOf(xyz_p[2]).doubleValue();
					validparam = true;
				}
				if (params[0].equals("beampower")) {
					beamPower = Double.valueOf(params[1]);
					validparam = true;
				}
				if (params[0].equals("multimaterial")) {
					multimaterial = Boolean.valueOf(params[1]).booleanValue();
					validparam = true;
				}
				if (params[0].equals("enable_temperature_dependent_nk")) {
					final String[] params_nk = params[1].split(" ");
					final String p = params_nk[0];
					enable_temperature_dependent_nk = Boolean.valueOf(p).booleanValue();
					validparam = true;
					if (params_nk.length > 1) {
						System.out.println("got, n,k list");
						final String ns = params_nk[1].substring(1, params_nk[1].length() - 2);
						System.out.println("ns: " + ns);
						final String ks = params_nk[2].substring(1, params_nk[2].length() - 2);
						System.out.println("ks: " + ks);

					}
				}
				if (validparam == false) {
					System.out.println("Warning: Line '" + line + "' not parsed.");
				}

			}
			br.close();
		}
		catch (final Exception e) {
			System.out.println("Probleme beim einlesen der Parameter");
			System.out.println("Einlesen " + e.toString());
		}

	}

	/**
	 * Methode gibt einen String mit wichtigen Eigenschaften des Dreiecks zurueck
	 */
	@Override
	public String toString() {
		String s = new String();
		s += "Parameters:\n";
		s += "Filename:                          =" + filename + "\n";
		s += "Number of reflections              =" + number_of_reflections + "\n";
		s += "Number of beams                    =" + number_of_beams + "\n";
		s += "Polarization                       =" + polarization + "\n";
		s += "Number of sublevels                =" + number_of_sublevels + "\n";
		s += "Incoming beam position vector      =" + incoming_beam_position_vector[0] + "," + incoming_beam_position_vector[1] + "," + incoming_beam_position_vector[2] + "\n";
		s += "Incoming beam direction vector     =" + incoming_beam_direction_vector[0] + "," + incoming_beam_direction_vector[1] + "," + incoming_beam_direction_vector[2] + "\n";
		s += "Beam waist radius                  =" + beam_waist_radius + "\n";
		s += "Beam wavelength                    =" + beam_wavelength + "\n";
		s += "Beam M-square                      =" + beam_M2 + "\n";
		s += "Beam distance to focus             =" + beam_distance_to_focus + "\n";
		s += "Beam source type                   =" + beam_source_type + "\n";
		s += "Beam source radius (if circle)     =" + beam_source_radius + "\n";
		s += "Beam source width x (if rectangle) =" + beam_source_width_x + "\n";
		s += "Beam source width y (if rectangle) =" + beam_source_width_y + "\n";
		s += "Beam intensity distribution        =" + beam_intensity_distribution + "\n";
		s += "Minimal propagation distance       =" + beam_min_distance + "\n";
		s += "Material name                      =" + material_name + "\n";
		s += "Material index of refraction r     =" + material_IOR_real + "\n";
		s += "Material index of refraction i     =" + material_IOR_imag + "\n";
		// MAterial 2 FF2
		s += "Material name2                      =" + material_name2 + "\n";
		s += "Material index of refraction r2     =" + material_IOR_real2 + "\n";
		s += "Material index of refraction i2     =" + material_IOR_imag2 + "\n";
		//
		s += "Bounding box spacing               =" + bounding_box_spacing + "\n";
		s += "Write CSV                          =" + target_csv + "\n";
		s += "Write PLY linear                   =" + target_ply_lin + "\n";
		s += "Write PLY logarithmic              =" + target_ply_log + "\n";
		s += "Number of Threads                  =" + number_of_threads + "\n";
		s += "beamPower			             =" + beamPower + "\n";
		s += "multimaterial			             =" + multimaterial + "\n";
		s += "enable_temperature_dependent_nk	 =" + enable_temperature_dependent_nk + "\n";
		return s;
	}

}
