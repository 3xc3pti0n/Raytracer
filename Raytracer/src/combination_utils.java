import java.io.File;
import java.io.FileWriter;
import java.io.PrintWriter;

public class combination_utils {
	/*
	 * holds the functionlity for the combination of several models.
	 * of which the last is the convMeltFoam.xxx model
	 */

	public static void run_continue_with_sclarTransportFOAM_melt(File rundir, String IDString, TriangleSet tset){
		/*
		 * in the rundir a completed case is expected. this case has a last, solved timestep, which is mapped to 0. then the new solver is applied.
		 * the system and the constant directory are taken from the FEM/models/scalarTransportFOAM_melt directory.
		 */
		// exchange the everything besides the constant/polyMesh directory
		Algorithmen.copy_file(new File("FEM/models/Combination/constant/transportProperties"), new File(rundir,"constant/transportProperties") );
		Algorithmen.copy_file(new File("FEM/models/Combination/constant/g"), new File(rundir,"constant/g") );
		Algorithmen.copy_dir(new File("FEM/models/Combination/system"), new File(rundir,"system") );
		
		// copy the missing field dictionaries 
		Algorithmen.copy_file(new File("FEM/models/Combination/0/p_rgh"), new File(rundir, "/0/p_rgh"));
		Algorithmen.copy_file(new File("FEM/models/Combination/0/alpha1"), new File(rundir, "/0/alpha1"));
		Algorithmen.copy_file(new File("FEM/models/Combination/0/alpha3"), new File(rundir, "/0/alpha3"));
		Algorithmen.copy_file(new File("FEM/models/Combination/0/h"), new File(rundir, "/0/h"));
		
		// rewriting them to match geometry
		write_boudaries_alpha1(rundir.getAbsolutePath() + "/0/alpha1", tset);
		write_boudaries_h(rundir.getAbsolutePath() + "/0/h", tset);
		write_boudaries_p(rundir.getAbsolutePath() + "/0/p", tset);
		write_boudaries_p_rgh(rundir.getAbsolutePath() + "/0/p_rgh", tset);
		
		// map the last timestep to 0
		Algorithmen.callbatchfile("mapfields_combination_after_scalarTransportFOAM", rundir.getAbsolutePath(), "", false);
		// change run-Settings to suit new solver
		Algorithmen.controldict_set_tend_dt(new File(rundir, "system/controlDict"), "0.3", "0.05");
		// run scalarTransportFoam_melt
		Algorithmen.callbatchfile("scalarTransportFOAM_melt.bat", rundir.getAbsolutePath(), "", false);
		
	}
	public static void run_continue_with_convMeltFoam(File rundir, String IDString, TriangleSet tset){
		/*
		 * in the rundir a completed case is expected. this case has a last, solved timestep, which is mapped to 0. then the new solver is applied.
		 * the system and the constant directory are taken from the FEM/models/scalarTransportFOAM_melt directory.
		 */
		// exchange the everything besides the constant/polyMesh directory
		Algorithmen.copy_file(new File("FEM/models/Combination/constant/transportProperties"), new File(rundir,"constant/transportProperties") );
		Algorithmen.copy_file(new File("FEM/models/Combination/constant/g"), new File(rundir,"constant/g") );
		Algorithmen.copy_dir(new File("FEM/models/Combination/system"), new File(rundir,"system") );
		
		
		// copy the missing field dictionaries 
		// ... all fields should be copied before in the scalarTransport_melt......
		
		// rewriting them to match geometry
		
		
		// map the last timestep to 0		
		Algorithmen.callbatchfile("mapfields_combination_after_scalarTransportFOAM_melt", rundir.getAbsolutePath(), "", false);
		
		// change run-Settings to suit new solver
		Algorithmen.controldict_set_tend_dt(new File(rundir, "system/controlDict"), "0.3", "0.005");
		
		// run scalarTransportFoam_melt
		Algorithmen.callbatchfile("convMeltFoam.bat", rundir.getAbsolutePath(), "", false);
		
	}
	
	public static void write_boudaries_alpha1( String outputfile, TriangleSet tri_set) {

		int bereich_zaehler;

		try {

			final File datei = new File(outputfile);
			final FileWriter ausgabestrom = new FileWriter(datei);
			final PrintWriter ausgabe = new PrintWriter(ausgabestrom);

			ausgabe.println("FoamFile");
			ausgabe.println("{");
			ausgabe.println("    " + "version" + "      " + "2.1;");
			ausgabe.println("    " + "format" + "       " + "ascii;");
			ausgabe.println("    " + "class" + "        " + "volScalarField;");
			ausgabe.println("    " + "location" + "     " + "\"0\";");
			ausgabe.println("    " + "object" + "       " + "U;");
			ausgabe.println("}");
			ausgabe.println("");
			ausgabe.println("dimensions    [ 0 0 0 0 0 0 0 ];    //K");
			ausgabe.println("internalField   uniform 0;"); // 0 0 0

			ausgabe.println("");
			ausgabe.println("boundaryField");
			ausgabe.println("{");
			ausgabe.println("");
			ausgabe.println("  ground");
			ausgabe.println("  {");
			ausgabe.println("      type calculated;");
			ausgabe.println("      value uniform 0;");
			ausgabe.println("  }");
			ausgabe.println("");
			ausgabe.println("  maxZ");
			ausgabe.println("  {");
			ausgabe.println("      type calculated;");
			ausgabe.println("      value uniform 0;");
			ausgabe.println("  }");
			ausgabe.println("");
			ausgabe.println("  minX");
			ausgabe.println("  {");
			ausgabe.println("      type calculated;");
			ausgabe.println("      value uniform 0;");
			ausgabe.println("  }");
			ausgabe.println("");
			ausgabe.println("  maxX");
			ausgabe.println("  {");
			ausgabe.println("      type calculated;");
			ausgabe.println("      value uniform 0;");
			ausgabe.println("  }");
			ausgabe.println("");
			ausgabe.println("  minY");
			ausgabe.println("  {");
			ausgabe.println("      type calculated;");
			ausgabe.println("      value uniform 0;");
			ausgabe.println("  }");
			ausgabe.println("");
			ausgabe.println("  maxY");
			ausgabe.println("  {");
			ausgabe.println("      type calculated;");
			ausgabe.println("      value uniform 0;");
			ausgabe.println("  }");
			ausgabe.println("");

			for (int i = 0; i < tri_set.getNumber(); i++) {
				bereich_zaehler = i + 1;

				ausgabe.println("  myfile_bereich" + bereich_zaehler);
				ausgabe.println("  {");
				ausgabe.println("      type calculated;");
				ausgabe.println("      value uniform 0;");

				ausgabe.println("  }");
				ausgabe.println("");

			}

			ausgabe.println("}");

			ausgabe.close();
			ausgabestrom.close();

		}
		catch (final Exception e) {
			System.out.println("Einlesen " + e.toString());
		}

	}
	
	public static void write_boudaries_h(String outputfile, TriangleSet tri_set) {

		int bereich_zaehler;

		try {

			final File datei = new File(outputfile);
			final FileWriter ausgabestrom = new FileWriter(datei);
			final PrintWriter ausgabe = new PrintWriter(ausgabestrom);

			ausgabe.println("FoamFile");
			ausgabe.println("{");
			ausgabe.println("    " + "version" + "      " + "2.1;");
			ausgabe.println("    " + "format" + "       " + "ascii;");
			ausgabe.println("    " + "class" + "        " + "volScalarField;");
			ausgabe.println("    " + "location" + "     " + "\"0\";");
			ausgabe.println("    " + "object" + "       " + "alpha1;");
			ausgabe.println("}");
			ausgabe.println("");
			ausgabe.println("dimensions    [0 2 -2 0 0 0 0];    //K");
			ausgabe.println("internalField   uniform -688275;"); // 0 0 0

			ausgabe.println("");
			ausgabe.println("boundaryField");
			ausgabe.println("{");
			ausgabe.println("");
			ausgabe.println("  ground");
			ausgabe.println("  {");
			ausgabe.println("      type calculated;");
			ausgabe.println("      value -688275;");
			ausgabe.println("  }");
			ausgabe.println("");
			ausgabe.println("  maxZ");
			ausgabe.println("  {");
			ausgabe.println("      type calculated;");
			ausgabe.println("      value -688275;");
			ausgabe.println("  }");
			ausgabe.println("");
			ausgabe.println("  minX");
			ausgabe.println("  {");
			ausgabe.println("      type calculated;");
			ausgabe.println("      value -688275;");
			ausgabe.println("  }");
			ausgabe.println("");
			ausgabe.println("  maxX");
			ausgabe.println("  {");
			ausgabe.println("      type calculated;");
			ausgabe.println("      value uniform 0;");
			ausgabe.println("  }");
			ausgabe.println("");
			ausgabe.println("  minY");
			ausgabe.println("  {");
			ausgabe.println("      type calculated;");
			ausgabe.println("      value uniform 0;");
			ausgabe.println("  }");
			ausgabe.println("");
			ausgabe.println("  maxY");
			ausgabe.println("  {");
			ausgabe.println("      type calculated;");
			ausgabe.println("      value -688275;");
			ausgabe.println("  }");
			ausgabe.println("");

			for (int i = 0; i < tri_set.getNumber(); i++) {
				bereich_zaehler = i + 1;

				ausgabe.println("  myfile_bereich" + bereich_zaehler);
				ausgabe.println("  {");
				ausgabe.println("      type calculated;");
				ausgabe.println("      value -688275;");

				ausgabe.println("  }");
				ausgabe.println("");

			}

			ausgabe.println("}");

			ausgabe.close();
			ausgabestrom.close();

		}
		catch (final Exception e) {
			System.out.println("Einlesen " + e.toString());
		}

	}
	
	public static void write_boudaries_p(String outputfile, TriangleSet tri_set) {

		int bereich_zaehler;

		try {

			final File datei = new File(outputfile);
			final FileWriter ausgabestrom = new FileWriter(datei);
			final PrintWriter ausgabe = new PrintWriter(ausgabestrom);

			ausgabe.println("FoamFile");
			ausgabe.println("{");
			ausgabe.println("    " + "version" + "      " + "2.1;");
			ausgabe.println("    " + "format" + "       " + "ascii;");
			ausgabe.println("    " + "class" + "        " + "volScalarField;");
			ausgabe.println("    " + "location" + "     " + "\"0\";");
			ausgabe.println("    " + "object" + "       " + "p;");
			ausgabe.println("}");
			ausgabe.println("");
			ausgabe.println("dimensions    [0 2 -2 0 0 0 0];    //K");
			ausgabe.println("internalField   uniform 0;"); // 0 0 0

			ausgabe.println("");
			ausgabe.println("boundaryField");
			ausgabe.println("{");
			ausgabe.println("");
			ausgabe.println("  ground");
			ausgabe.println("  {");
			ausgabe.println("      type zeroGradient;");
			ausgabe.println("  }");
			ausgabe.println("");
			ausgabe.println("  maxZ");
			ausgabe.println("  {");
			ausgabe.println("      type zeroGradient;");
			ausgabe.println("  }");
			ausgabe.println("");
			ausgabe.println("  minX");
			ausgabe.println("  {");
			ausgabe.println("      type zeroGradient;");
			ausgabe.println("  }");
			ausgabe.println("");
			ausgabe.println("  maxX");
			ausgabe.println("  {");
			ausgabe.println("      type fixedValue; value uniform 0;");
			ausgabe.println("  }");
			ausgabe.println("");
			ausgabe.println("  minY");
			ausgabe.println("  {");
			ausgabe.println("      type symmetryPlane;");
			ausgabe.println("  }");
			ausgabe.println("");
			ausgabe.println("  maxY");
			ausgabe.println("  {");
			ausgabe.println("      type zeroGradient;");
			ausgabe.println("  }");
			ausgabe.println("");

			for (int i = 0; i < tri_set.getNumber(); i++) {
				bereich_zaehler = i + 1;

				ausgabe.println("  myfile_bereich" + bereich_zaehler);
				ausgabe.println("  {");
				ausgabe.println("      type zeroGradient;");
				ausgabe.println("  }");
				ausgabe.println("");

			}

			ausgabe.println("}");

			ausgabe.close();
			ausgabestrom.close();

		}
		catch (final Exception e) {
			System.out.println("Einlesen " + e.toString());
		}

	}
	
	public static void write_boudaries_p_rgh(String outputfile, TriangleSet tri_set) {

		int bereich_zaehler;

		try {

			final File datei = new File(outputfile);
			final FileWriter ausgabestrom = new FileWriter(datei);
			final PrintWriter ausgabe = new PrintWriter(ausgabestrom);

			ausgabe.println("FoamFile");
			ausgabe.println("{");
			ausgabe.println("    " + "version" + "      " + "2.1;");
			ausgabe.println("    " + "format" + "       " + "ascii;");
			ausgabe.println("    " + "class" + "        " + "volScalarField;");
			ausgabe.println("    " + "location" + "     " + "\"0\";");
			ausgabe.println("    " + "object" + "       " + "p_rgh;");
			ausgabe.println("}");
			ausgabe.println("");
			ausgabe.println("dimensions    [0 2 -2 0 0 0 0];    //K");
			ausgabe.println("internalField   uniform 0;"); // 0 0 0

			ausgabe.println("");
			ausgabe.println("boundaryField");
			ausgabe.println("{");
			ausgabe.println("");
			ausgabe.println("  ground");
			ausgabe.println("  {");
			ausgabe.println("      type zeroGradient;");
			ausgabe.println("  }");
			ausgabe.println("");
			ausgabe.println("  maxZ");
			ausgabe.println("  {");
			ausgabe.println("      type zeroGradient;");
			ausgabe.println("  }");
			ausgabe.println("");
			ausgabe.println("  minX");
			ausgabe.println("  {");
			ausgabe.println("      type zeroGradient;");
			ausgabe.println("  }");
			ausgabe.println("");
			ausgabe.println("  maxX");
			ausgabe.println("  {");
			ausgabe.println("      type fixedValue; value uniform 0;");
			ausgabe.println("  }");
			ausgabe.println("");
			ausgabe.println("  minY");
			ausgabe.println("  {");
			ausgabe.println("      type symmetryPlane;");
			ausgabe.println("  }");
			ausgabe.println("");
			ausgabe.println("  maxY");
			ausgabe.println("  {");
			ausgabe.println("      type zeroGradient;");
			ausgabe.println("  }");
			ausgabe.println("");

			for (int i = 0; i < tri_set.getNumber(); i++) {
				bereich_zaehler = i + 1;

				ausgabe.println("  myfile_bereich" + bereich_zaehler);
				ausgabe.println("  {");
				ausgabe.println("      type zeroGradient;");
				ausgabe.println("  }");
				ausgabe.println("");

			}

			ausgabe.println("}");

			ausgabe.close();
			ausgabestrom.close();

		}
		catch (final Exception e) {
			System.out.println("Einlesen " + e.toString());
		}

	}
	
}
