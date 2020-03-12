import java.awt.BorderLayout;
import java.awt.EventQueue;
import java.awt.Graphics2D;
import java.awt.RenderingHints;
import java.awt.image.BufferedImage;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.text.DecimalFormat;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Calendar;
import java.util.List;

import javax.imageio.ImageIO;
import javax.swing.JFrame;
import javax.swing.JPanel;
import javax.swing.border.EmptyBorder;
import javax.swing.filechooser.FileNameExtensionFilter;

import javax.swing.ImageIcon;
import javax.swing.JButton;
import javax.swing.JFileChooser;
import javax.swing.JLabel;
import javax.swing.JTextField;
import javax.swing.JCheckBox;
import javax.swing.UIManager;
import javax.swing.UnsupportedLookAndFeelException;

import java.awt.event.MouseAdapter;
import java.awt.event.MouseEvent;
import java.awt.Font;
import java.awt.Color;

import javax.swing.JTabbedPane;

import java.awt.FlowLayout;

import javax.swing.JSeparator;

import java.awt.SystemColor;

import javax.swing.SwingConstants;

import java.awt.event.ActionListener;
import java.awt.event.ActionEvent;

import javax.swing.event.ChangeListener;
import javax.swing.event.ChangeEvent;
import javax.swing.JLayeredPane;
import java.awt.Component;
import javax.swing.Box;
import java.awt.Label;
import java.awt.TextField;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;

public class Main extends JFrame {

	private JPanel contentPane;

	public List<Thread> threadList = new ArrayList<Thread>();
	public JPanel panel_1;
	public JTextField text_number_beams;
	public JTextField text_number_reflections;
	public JTextField text_number_threads;
	public JTextField text_Polarisation;
	public JTextField text_beam_waist_radius;
	public JTextField text_beam_wavelength;
	public JTextField text_beam_m_square;
	public JTextField text_number_sublevel;
	public JTextField text_beam_source_width_y;
	public JTextField text_beam_source_width_x;
	public JTextField text_incoming_pos_vector;
	public JTextField text_beam_source_typ;
	public JTextField text_beam_source_radius;
	public JTextField text_incoming_direction_vector;
	public JTextField text_bounding_box_spacing;
	public JTextField text_beam_intensity_dist;
	public JTextField text_write_csv;
	public JTextField text_beam_distance_to_focus;
	public JTextField text_write_ply_log;
	public JTextField text_beam_min_dist;
	public JTextField text_rel_diff_scat;

	public File current_STL_file = new File("./Raytracer_local_input/TestEbene.stl");
	public File current_caseDir = new File("./");
	public File current_STL_detector = new File("./Raytracer_local_input/TestDetector.stl");
	public File current_outDir = new File("./Raytracer_local_output");
	public File[] listOfFiles;

	public JTextField txt_currentSTL;
	public JTextField txtOut;
	public JTextField text_write_ply_lin;
	public JTextField text_beam_material_name;
	public JTextField text_beam_IOR_real;
	public JTextField text_beam_IOR_imag;
	private JTextField nowAtTextField;
	private JTextField txt_currentDetector;
	private JTextField txtCraytrtheraytracerbatchtracepystls;
	private JTextField text_beamPower;

	public static BufferedImage resize(BufferedImage image, int width, int height) {
		BufferedImage bi = new BufferedImage(width, height, BufferedImage.TRANSLUCENT);
		Graphics2D g2d = (Graphics2D) bi.createGraphics();
		g2d.addRenderingHints(new RenderingHints(RenderingHints.KEY_RENDERING, RenderingHints.VALUE_RENDER_QUALITY));
		g2d.drawImage(image, 0, 0, width, height, null);
		g2d.dispose();
		return bi;
	}

	/**
	 * Launch the application.
	 */
	public static void main(String[] args) {
		try {
			UIManager.setLookAndFeel(UIManager.getSystemLookAndFeelClassName());
		}
		catch (ClassNotFoundException e1) {
			// TODO Auto-generated catch block
			e1.printStackTrace();
		}
		catch (InstantiationException e1) {
			// TODO Auto-generated catch block
			e1.printStackTrace();
		}
		catch (IllegalAccessException e1) {
			// TODO Auto-generated catch block
			e1.printStackTrace();
		}
		catch (UnsupportedLookAndFeelException e1) {
			// TODO Auto-generated catch block
			e1.printStackTrace();
		}

		EventQueue.invokeLater(new Runnable() {
			public void run() {
				try {
					Main frame = new Main();
					frame.setVisible(true);
				}
				catch (Exception e) {
					e.printStackTrace();
				}
			}
		});
	}

	/**
	 * Create the frame.
	 */
	public Main() {
		setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
		setBounds(100, 100, 1080, 769);
		contentPane = new JPanel();
		contentPane.setBorder(new EmptyBorder(5, 5, 5, 5));
		setContentPane(contentPane);
		contentPane.setLayout(null);

		final JTabbedPane tabbedPane = new JTabbedPane(JTabbedPane.TOP);

		tabbedPane.setToolTipText("AS");
		tabbedPane.setBounds(0, 0, 1043, 713);
		contentPane.add(tabbedPane);

		JPanel panel_4 = new JPanel();
		tabbedPane.addTab("Raytracing", null, panel_4, null);
		panel_4.setLayout(null);

		JLabel lblNewLabel = new JLabel("Number of Rays");
		lblNewLabel.setBounds(385, 18, 156, 14);
		panel_4.add(lblNewLabel);

		text_number_beams = new JTextField();
		text_number_beams.setBounds(588, 18, 73, 20);
		panel_4.add(text_number_beams);
		text_number_beams.setColumns(10);

		JLabel lblNumberOfReflections = new JLabel("Number of Reflections");
		lblNumberOfReflections.setBounds(385, 51, 184, 14);
		panel_4.add(lblNumberOfReflections);

		text_number_reflections = new JTextField();
		text_number_reflections.setColumns(10);
		text_number_reflections.setBounds(588, 51, 73, 20);
		panel_4.add(text_number_reflections);

		JLabel lblNOfThreads = new JLabel("N of threads");
		lblNOfThreads.setBounds(385, 85, 184, 14);
		panel_4.add(lblNOfThreads);

		text_number_threads = new JTextField();
		text_number_threads.setColumns(10);
		text_number_threads.setBounds(588, 82, 73, 20);
		panel_4.add(text_number_threads);

		text_Polarisation = new JTextField();
		text_Polarisation.setColumns(10);
		text_Polarisation.setBounds(388, 427, 73, 20);
		panel_4.add(text_Polarisation);

		JLabel lblPolarisationlinearxLineary = new JLabel("Polarisation: (linear_x, linear_y, azimutal, zirkular)");
		lblPolarisationlinearxLineary.setBounds(15, 423, 358, 29);
		panel_4.add(lblPolarisationlinearxLineary);

		JLabel lblBeamWaistRadius = new JLabel("Beam waist radius");
		lblBeamWaistRadius.setBounds(15, 214, 136, 14);
		panel_4.add(lblBeamWaistRadius);

		text_beam_waist_radius = new JTextField();
		text_beam_waist_radius.setColumns(10);
		text_beam_waist_radius.setBounds(244, 217, 73, 20);
		panel_4.add(text_beam_waist_radius);

		text_beam_wavelength = new JTextField();
		text_beam_wavelength.setColumns(10);
		text_beam_wavelength.setBounds(244, 250, 73, 20);
		panel_4.add(text_beam_wavelength);

		JLabel lblWavelength = new JLabel("Wavelength");
		lblWavelength.setBounds(15, 247, 90, 14);
		panel_4.add(lblWavelength);

		JLabel lblM = new JLabel("M^2");
		lblM.setBounds(15, 183, 90, 14);
		panel_4.add(lblM);

		text_beam_m_square = new JTextField();
		text_beam_m_square.setColumns(10);
		text_beam_m_square.setBounds(244, 186, 73, 20);
		panel_4.add(text_beam_m_square);

		text_number_sublevel = new JTextField();
		text_number_sublevel.setColumns(10);
		text_number_sublevel.setBounds(588, 365, 73, 20);
		panel_4.add(text_number_sublevel);

		JLabel lblSublevels = new JLabel("Sublevels");
		lblSublevels.setBounds(385, 371, 90, 14);
		panel_4.add(lblSublevels);

		text_beam_source_width_y = new JTextField();
		text_beam_source_width_y.setColumns(10);
		text_beam_source_width_y.setBounds(244, 354, 73, 20);
		panel_4.add(text_beam_source_width_y);

		JLabel lblBeamSourceWidth_1 = new JLabel("Beam source width (y)");
		lblBeamSourceWidth_1.setBounds(15, 357, 184, 14);
		panel_4.add(lblBeamSourceWidth_1);

		JLabel lblBeamSourceWidth = new JLabel("Beam source width (x)");
		lblBeamSourceWidth.setBounds(15, 324, 167, 14);
		panel_4.add(lblBeamSourceWidth);

		text_beam_source_width_x = new JTextField();
		text_beam_source_width_x.setColumns(10);
		text_beam_source_width_x.setBounds(244, 321, 73, 20);
		panel_4.add(text_beam_source_width_x);

		text_incoming_pos_vector = new JTextField();
		text_incoming_pos_vector.setColumns(10);
		text_incoming_pos_vector.setBounds(244, 49, 113, 20);
		panel_4.add(text_incoming_pos_vector);

		JLabel lblBeamSourceType = new JLabel("Beam source type");
		lblBeamSourceType.setBounds(15, 288, 136, 14);
		panel_4.add(lblBeamSourceType);

		text_beam_source_typ = new JTextField();
		text_beam_source_typ.setColumns(10);
		text_beam_source_typ.setBounds(244, 281, 73, 20);
		panel_4.add(text_beam_source_typ);

		JLabel lblBeamSourcePosition = new JLabel("Beam source position");
		lblBeamSourcePosition.setBounds(15, 49, 179, 14);
		panel_4.add(lblBeamSourcePosition);

		JLabel lblBeamSourceRadius = new JLabel("Beam source radius [\u00B5m]:\r\nbeacause not all generated \r\nrays are inside one radius");
		lblBeamSourceRadius.setBounds(385, 289, 529, 20);
		panel_4.add(lblBeamSourceRadius);

		text_beam_source_radius = new JTextField();
		text_beam_source_radius.setColumns(10);
		text_beam_source_radius.setBounds(945, 288, 55, 21);
		panel_4.add(text_beam_source_radius);

		text_incoming_direction_vector = new JTextField();
		text_incoming_direction_vector.setColumns(10);
		text_incoming_direction_vector.setBounds(244, 16, 113, 20);
		panel_4.add(text_incoming_direction_vector);

		JLabel lblIncomingBeamDirection = new JLabel("Incoming beam direction vector");
		lblIncomingBeamDirection.setBounds(15, 16, 225, 14);
		panel_4.add(lblIncomingBeamDirection);

		JLabel lblBoundingBoxSpacing = new JLabel("Bounding box spacing");
		lblBoundingBoxSpacing.setBounds(385, 326, 167, 29);
		panel_4.add(lblBoundingBoxSpacing);

		text_bounding_box_spacing = new JTextField();
		text_bounding_box_spacing.setColumns(10);
		text_bounding_box_spacing.setBounds(588, 325, 86, 20);
		panel_4.add(text_bounding_box_spacing);

		text_beam_intensity_dist = new JTextField();
		text_beam_intensity_dist.setColumns(10);
		text_beam_intensity_dist.setBounds(244, 113, 86, 20);
		panel_4.add(text_beam_intensity_dist);

		JLabel lblBeamPowerDistribution = new JLabel("Beam power distribution");
		lblBeamPowerDistribution.setBounds(15, 113, 167, 14);
		panel_4.add(lblBeamPowerDistribution);

		text_write_csv = new JTextField();
		text_write_csv.setColumns(10);
		text_write_csv.setBounds(588, 203, 86, 20);
		panel_4.add(text_write_csv);

		JLabel lblWritecsv = new JLabel("Write .csv");
		lblWritecsv.setBounds(385, 203, 90, 14);
		panel_4.add(lblWritecsv);

		JLabel lblWriteLinply = new JLabel("beam disance to focus");
		lblWriteLinply.setBounds(15, 83, 179, 14);
		panel_4.add(lblWriteLinply);

		text_beam_distance_to_focus = new JTextField();
		text_beam_distance_to_focus.setColumns(10);
		text_beam_distance_to_focus.setBounds(244, 83, 86, 20);
		panel_4.add(text_beam_distance_to_focus);

		text_write_ply_log = new JTextField();
		text_write_ply_log.setColumns(10);
		text_write_ply_log.setBounds(588, 172, 86, 20);
		panel_4.add(text_write_ply_log);

		JLabel lblWriteLogply = new JLabel("Write LOG .ply");
		lblWriteLogply.setBounds(385, 175, 90, 14);
		panel_4.add(lblWriteLogply);

		text_beam_min_dist = new JTextField();
		text_beam_min_dist.setColumns(10);
		text_beam_min_dist.setBounds(588, 113, 73, 20);
		panel_4.add(text_beam_min_dist);

		JLabel lblMinimumBeamDistance = new JLabel("Minimum beam distance");
		lblMinimumBeamDistance.setBounds(385, 116, 167, 14);
		panel_4.add(lblMinimumBeamDistance);

		JLabel lblDiffuseScattering = new JLabel("% diffuse scattering");
		lblDiffuseScattering.setBounds(15, 390, 172, 14);
		panel_4.add(lblDiffuseScattering);

		text_rel_diff_scat = new JTextField();
		text_rel_diff_scat.setColumns(10);
		text_rel_diff_scat.setBounds(244, 390, 73, 20);
		panel_4.add(text_rel_diff_scat);

		final JButton btn_trace_single = new JButton("Trace single file");
		btn_trace_single.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
			}
		});
		btn_trace_single.addMouseListener(new MouseAdapter() {
			@Override
			public void mouseClicked(MouseEvent e) {
				btn_trace_single.setText("Tracing for local"); // nur Gui stuff
				Color bg = new Color(1.f, 0.f, 0.f);
				btn_trace_single.setBackground(bg);

				Raytracer rayt; // rechnet Fresnel - orig. Michalowski
				Raytracer_diffus rayt_diff; // rechnet die streuung FFE, Lambert

				// first raytrace starting geometry
				String file_to_trace = current_STL_file.getAbsolutePath();
				String file_to_trace_detector = current_STL_detector.getAbsolutePath();
				System.out.println("Tracing locally: " + file_to_trace);
				String file_name = current_STL_file.getName().substring(0, current_STL_file.getName().length() - 4); // endung entfernen
				if (false) { // test zum detecten mit true, ohne detect: false
					Algorithmen.stl_to_tri(file_to_trace, "Raytracer_local_input/target_" + file_name + ".tri"); // convertieren in tri // in
					Algorithmen.stl_to_tri_detector(file_to_trace_detector, "Raytracer_local_input/detect_" + file_name + ".tri"); // convertieren in tri // in
					Algorithmen.appendfiles(new File("Raytracer_local_input/" + file_name + ".tri"), new File("Raytracer_local_input/target_" + file_name + ".tri"), new File("Raytracer_local_input/detect_" + file_name + ".tri"));
				}
				else { // ohne detetct
					Algorithmen.stl_to_tri(file_to_trace, "Raytracer_local_input/" + file_name + ".tri"); // convertieren in tri // in
				}
				String[] Steuer_Kommandos = new String[3];
				Steuer_Kommandos[0] = "Raytracer_local_output/" + file_name;// Name der ausgabe
				Steuer_Kommandos[1] = "Parameter_Raytracer.txt"; // Name der Parameter
				Steuer_Kommandos[2] = "Raytracer_local_input/" + file_name + ".tri";
				rayt = new Raytracer(Steuer_Kommandos); // ray-trace,
				// diffuse
				if (Double.valueOf(text_rel_diff_scat.getText()) > 0) {
					System.out.println("Now diffuse tracing");
					Steuer_Kommandos[0] = "Raytracer_local_output/" + "diff_" + file_name;// Name
					rayt_diff = new Raytracer_diffus(Steuer_Kommandos); // ray-trace
					System.out.println("Combining diffuse and fresnel");
					File datei = new File("Raytracer_local_output/comb_" + file_name + ".ply");
					FileWriter ausgabestrom;
					try { // combining
						ausgabestrom = new FileWriter(datei);
						PrintWriter combined_output = new PrintWriter(ausgabestrom);
						Algorithmen.write_combined_output_as_ply(rayt.get_global_triSet()[0], rayt_diff.get_global_triSet(), "FETZ", Double.valueOf(text_rel_diff_scat.getText()) * 0.01, combined_output);
						System.out.print("Combined - DONE");
					}
					catch (IOException e1) {
						// TODO Auto-generated catch block
						e1.printStackTrace();
					}
				}
				btn_trace_single.setText("Done Tracing");
				bg = new Color(0.f, 1.f, 0.f);
				btn_trace_single.setBackground(bg);
			}
		});

		btn_trace_single.setFont(new Font("Tahoma", Font.BOLD, 14));
		btn_trace_single.setBounds(882, 616, 145, 47);
		panel_4.add(btn_trace_single);

		JButton btnSelectCurrentStl = new JButton("Select STL file");
		btnSelectCurrentStl.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
			}
		});
		btnSelectCurrentStl.addMouseListener(new MouseAdapter() {
			@Override
			public void mouseClicked(MouseEvent e) {
				JFileChooser chooser = new JFileChooser();
				chooser.setCurrentDirectory(current_caseDir);
				FileNameExtensionFilter filter = new FileNameExtensionFilter("STL Files", "stl", "txt");
				chooser.setFileFilter(filter);
				int returnVal = chooser.showOpenDialog(null);
				if (returnVal == JFileChooser.APPROVE_OPTION) {
					// txt_currentSTL.setText(chooser.getSelectedFile().toString());
					current_STL_file = chooser.getSelectedFile();
					txt_currentSTL.setText(current_STL_file.getAbsolutePath());
				}
			}

		});
		btnSelectCurrentStl.setFont(new Font("Tahoma", Font.BOLD, 9));
		btnSelectCurrentStl.setBounds(766, 511, 103, 29);
		panel_4.add(btnSelectCurrentStl);

		txt_currentSTL = new JTextField();
		txt_currentSTL.setColumns(10);
		txt_currentSTL.setBounds(251, 513, 505, 27);
		txt_currentSTL.setText("Raytracer_local_input/" + "TestEbene.stl");
		panel_4.add(txt_currentSTL);

		JLabel lblToTrace = new JLabel("to trace");
		lblToTrace.setFont(new Font("Tahoma", Font.BOLD, 16));
		lblToTrace.setBounds(251, 483, 122, 20);
		panel_4.add(lblToTrace);

		JButton btn_save_raytracer_parameter = new JButton("Save Parameters");
		btn_save_raytracer_parameter.addMouseListener(new MouseAdapter() {
			@Override
			public void mouseClicked(MouseEvent e) {
				write_raytracer_parameter(); // saving to the parameter file
			}
		});
		btn_save_raytracer_parameter.setFont(new Font("Tahoma", Font.BOLD, 14));
		btn_save_raytracer_parameter.setBounds(15, 483, 179, 47);
		panel_4.add(btn_save_raytracer_parameter);

		txtOut = new JTextField();
		txtOut.setText("Success save ??");
		txtOut.setBounds(15, 549, 195, 35);
		panel_4.add(txtOut);
		txtOut.setColumns(10);

		JLabel lblWriteFetzply = new JLabel("Write custom .ply");
		lblWriteFetzply.setBounds(385, 234, 90, 14);
		panel_4.add(lblWriteFetzply);

		text_write_ply_lin = new JTextField();
		text_write_ply_lin.setColumns(10);
		text_write_ply_lin.setBounds(588, 234, 86, 20);
		panel_4.add(text_write_ply_lin);

		JLabel lblMaterialName = new JLabel("Material 1 Name");
		lblMaterialName.setBounds(757, 22, 129, 14);
		panel_4.add(lblMaterialName);

		text_beam_material_name = new JTextField();
		text_beam_material_name.setColumns(10);
		text_beam_material_name.setBounds(901, 18, 82, 20);
		panel_4.add(text_beam_material_name);

		JLabel lblMatN = new JLabel("Mat 1: n");
		lblMatN.setBounds(757, 55, 90, 14);
		panel_4.add(lblMatN);

		text_beam_IOR_real = new JTextField();
		text_beam_IOR_real.setColumns(10);
		text_beam_IOR_real.setBounds(901, 51, 82, 20);
		panel_4.add(text_beam_IOR_real);

		JLabel lblMatK = new JLabel("Mat 1: k");
		lblMatK.setBounds(757, 89, 90, 14);
		panel_4.add(lblMatK);

		text_beam_IOR_imag = new JTextField();
		text_beam_IOR_imag.setColumns(10);
		text_beam_IOR_imag.setBounds(901, 82, 82, 20);
		panel_4.add(text_beam_IOR_imag);

		txt_currentDetector = new JTextField();
		txt_currentDetector.setText("Raytracer_local_input/TestDetector.stl");
		txt_currentDetector.setBounds(251, 571, 505, 29);
		panel_4.add(txt_currentDetector);
		txt_currentDetector.setColumns(10);

		JLabel lblToDetect = new JLabel("Detector");
		lblToDetect.setFont(new Font("Tahoma", Font.BOLD, 16));
		lblToDetect.setBounds(251, 545, 82, 20);
		panel_4.add(lblToDetect);

		JButton btnSelectCurrentDetector = new JButton("Select STL file");
		btnSelectCurrentDetector.addMouseListener(new MouseAdapter() {
			@Override
			public void mouseClicked(MouseEvent e) {
				JFileChooser chooser = new JFileChooser();
				chooser.setCurrentDirectory(current_caseDir);
				FileNameExtensionFilter filter = new FileNameExtensionFilter("STL Files", "stl", "txt");
				chooser.setFileFilter(filter);
				int returnVal = chooser.showOpenDialog(null);
				if (returnVal == JFileChooser.APPROVE_OPTION) {
					// txt_currentSTL.setText(chooser.getSelectedFile().toString());
					current_STL_detector = chooser.getSelectedFile();
					txt_currentDetector.setText(current_STL_detector.getAbsolutePath());
				}
			}
		});
		btnSelectCurrentDetector.setFont(new Font("Tahoma", Font.BOLD, 9));
		btnSelectCurrentDetector.setBounds(766, 571, 103, 29);
		panel_4.add(btnSelectCurrentDetector);

		JLabel lblBeamPowerw = new JLabel("Beam Power [W]");
		lblBeamPowerw.setBounds(15, 149, 167, 14);
		panel_4.add(lblBeamPowerw);

		text_beamPower = new JTextField();
		text_beamPower.setColumns(10);
		text_beamPower.setBounds(244, 149, 86, 20);
		panel_4.add(text_beamPower);

		panel_1 = new JPanel();
		tabbedPane.addTab("BatchTracing", null, panel_1, null);
		panel_1.setLayout(null);

		JButton BatchTraceBtn = new JButton("Trace Batch");
		BatchTraceBtn.addMouseListener(new MouseAdapter() {
			@Override
			public void mouseClicked(MouseEvent e) {
				//String fparam = "C:\\Users\\ifsw-ffetzer\\Desktop\\STL_Rek\\parameter.txt"; // hoer liegen die Parameter in  textfile: Name\tLeistung\df\n
				String fparam = "G:\\jonas_durchschweissen\\Rek\\STL\\paras.txt"; // Jonas
				System.out.println("Started Tracing");
				// Here Batch Processing...
				Raytracer rayt;
				Raytracer_diffus rayt_diff; // rechnet die streuung FFE, Lambert
				for (int i = 0; i < listOfFiles.length; i++) {
					if ((listOfFiles[i].getName().endsWith(".stl")) || (listOfFiles[i].getName().endsWith(".STL"))) {
						System.out.print("Now at...\n" + listOfFiles[i].getPath());
						nowAtTextField.setText("Now at...\n" + listOfFiles[i].getPath());
						// get Raytr_params
						adapt_Ray_Params(listOfFiles[i].getName(), fparam, "Parameter_Raytracer.txt"); // STL file, "Labbook", ParasRaytr.
						// STL to Tri
						Algorithmen.stl_to_tri(listOfFiles[i].getPath(), listOfFiles[i].getPath().replace(".stl", ".tri")); // convertieren in tri // in
						// Raytracen - Fresnel
						String[] Steuer_Kommandos = new String[3];
						Steuer_Kommandos[0] = "./Raytracer_local_output/" + listOfFiles[i].getName().replace(".stl", "");// Name der ausgabe
						Steuer_Kommandos[1] = "Parameter_Raytracer.txt"; // Name der Parameter
						Steuer_Kommandos[2] = listOfFiles[i].getAbsolutePath().replace(".stl", ".tri");
						rayt = new Raytracer(Steuer_Kommandos); // ray-trace,
						// Raytracen - diffuse
						if (Double.valueOf(text_rel_diff_scat.getText()) > 0) {
							System.out.println("Now diffuse tracing");
							Steuer_Kommandos[0] = "./Raytracer_local_output/" + "diff_" + listOfFiles[i].getName().replace(".stl", "");// Name der ausgabe
							rayt_diff = new Raytracer_diffus(Steuer_Kommandos); // ray-trace
							System.out.println("Combining diffuse and fresnel");
							File datei = new File("Raytracer_local_output/comb_" + listOfFiles[i].getName().replace(".stl", ".ply"));
							FileWriter ausgabestrom;
							try { // combining
								ausgabestrom = new FileWriter(datei);
								PrintWriter combined_output = new PrintWriter(ausgabestrom);
								Algorithmen.write_combined_output_as_ply(rayt.get_global_triSet()[0], rayt_diff.get_global_triSet(), "FETZ", Double.valueOf(text_rel_diff_scat.getText()) * 0.01, combined_output);
								System.out.print("Combined - DONE");
							}
							catch (IOException e1) {
								// TODO Auto-generated catch block
								e1.printStackTrace();
							}
						}
						System.out.print("DONE");
					}
				}
			}
		});

		BatchTraceBtn.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
			}
		});
		BatchTraceBtn.setBounds(326, 219, 170, 75);
		BatchTraceBtn.setFont(new Font("Tahoma", Font.BOLD, 16));
		panel_1.add(BatchTraceBtn);

		JLabel lblSetPropertiesIn = new JLabel("Set Properties in \"Raytracing\" Panel");
		lblSetPropertiesIn.setBounds(10, 11, 282, 19);
		lblSetPropertiesIn.setFont(new Font("Tahoma", Font.BOLD, 15));
		panel_1.add(lblSetPropertiesIn);

		JPanel panel = new JPanel();
		panel.setBounds(693, 47, 1, 1);
		panel.setLayout(null);
		panel.setBackground(new Color(100, 149, 237));
		panel_1.add(panel);

		JSeparator separator_2 = new JSeparator();
		separator_2.setOrientation(SwingConstants.VERTICAL);
		separator_2.setForeground(SystemColor.textHighlight);
		separator_2.setBounds(582, 193, 134, 1);
		panel_1.add(separator_2);

		JButton btn_test_stuff = new JButton("test stuff");
		btn_test_stuff.addMouseListener(new MouseAdapter() {
			@Override
			public void mouseClicked(MouseEvent e) {
				test_func();
			}
		});

		btn_test_stuff.setFont(new Font("Tahoma", Font.BOLD, 15));
		btn_test_stuff.setBounds(42, 588, 170, 75);
		panel_1.add(btn_test_stuff);

		nowAtTextField = new JTextField();
		nowAtTextField.setColumns(10);
		nowAtTextField.setBounds(31, 346, 647, 32);
		panel_1.add(nowAtTextField);

		JLabel lblSelectFolderContaining = new JLabel("Select Folder containing the STLs");
		lblSelectFolderContaining.setFont(new Font("Tahoma", Font.BOLD, 15));
		lblSelectFolderContaining.setBounds(10, 41, 292, 19);
		panel_1.add(lblSelectFolderContaining);

		JLabel lblAllStlsWill = new JLabel("all STLs will be traced");
		lblAllStlsWill.setFont(new Font("Tahoma", Font.BOLD, 15));
		lblAllStlsWill.setBounds(10, 71, 292, 19);
		panel_1.add(lblAllStlsWill);

		txtCraytrtheraytracerbatchtracepystls = new JTextField();
		txtCraytrtheraytracerbatchtracepystls.setText(".");
		txtCraytrtheraytracerbatchtracepystls.setBounds(31, 123, 452, 25);
		panel_1.add(txtCraytrtheraytracerbatchtracepystls);
		txtCraytrtheraytracerbatchtracepystls.setColumns(10);

		JButton DirButtonBatch = new JButton("Select Dir");
		DirButtonBatch.addMouseListener(new MouseAdapter() {
			@Override
			public void mouseClicked(MouseEvent arg0) {
				JFileChooser chooser = new JFileChooser();
				chooser.setCurrentDirectory(current_caseDir);
				chooser.setFileSelectionMode(JFileChooser.DIRECTORIES_ONLY);
				int returnVal = chooser.showOpenDialog(null);
				if (returnVal == JFileChooser.APPROVE_OPTION) {
					File current_Dir = chooser.getSelectedFile();
					System.out.println(current_Dir.getAbsolutePath());
					listOfFiles = current_Dir.listFiles();
					for (int i = 0; i < listOfFiles.length; i++) {
						// STL to Tri
						System.out.println(i + listOfFiles[i].getPath());
					}
				}
			}
		});
		DirButtonBatch.setBounds(532, 124, 123, 24);
		panel_1.add(DirButtonBatch);

		JLabel lblNowAt = new JLabel("Now at:");
		lblNowAt.setFont(new Font("Tahoma", Font.BOLD, 15));
		lblNowAt.setBounds(31, 319, 292, 19);
		panel_1.add(lblNowAt);

		tabbedPane.addChangeListener(new ChangeListener() {
			public void stateChanged(ChangeEvent e) {
				if (tabbedPane.getSelectedIndex() == 0) { // reads the Raytracer
															// Parameters if the
															// tab is selected
					System.out.println("Reading Raytracer Parameters");

					Parameters parameter = new Parameters("Parameter_Raytracer.txt");

					text_beam_distance_to_focus.setText("" + parameter.beam_distance_to_focus);
					text_beam_IOR_imag.setText("" + parameter.material_IOR_imag);
					text_beam_IOR_real.setText("" + parameter.material_IOR_real);
					text_beamPower.setText("" + parameter.beamPower);
					text_beam_m_square.setText("" + parameter.beam_M2);
					text_beam_material_name.setText("" + parameter.material_name);
					text_beam_min_dist.setText("" + parameter.beam_min_distance);
					text_beam_source_radius.setText("" + parameter.beam_source_radius);
					text_beam_source_typ.setText("" + parameter.beam_source_type);
					text_beam_source_width_x.setText("" + parameter.beam_source_width_x);
					text_beam_source_width_y.setText("" + parameter.beam_source_width_y);
					text_beam_waist_radius.setText("" + parameter.beam_waist_radius);
					text_beam_wavelength.setText("" + parameter.beam_wavelength);
					text_bounding_box_spacing.setText("" + parameter.bounding_box_spacing);
					text_incoming_direction_vector.setText("" + parameter.incoming_beam_direction_vector[0] + ", " + parameter.incoming_beam_direction_vector[1] + ", " + parameter.incoming_beam_direction_vector[2]);
					text_incoming_pos_vector.setText("" + parameter.incoming_beam_position_vector[0] + ", " + parameter.incoming_beam_position_vector[1] + ", " + parameter.incoming_beam_position_vector[2]);
					text_number_beams.setText("" + parameter.number_of_beams);
					text_number_reflections.setText("" + parameter.number_of_reflections);
					text_number_sublevel.setText("" + parameter.number_of_sublevels);
					text_number_threads.setText("" + parameter.number_of_threads);
					text_Polarisation.setText("" + parameter.polarization);
					text_write_csv.setText("" + parameter.target_csv);
					text_write_ply_lin.setText("" + parameter.target_ply_lin);
					text_write_ply_log.setText("" + parameter.target_ply_log);
					text_beam_intensity_dist.setText("" + parameter.beam_intensity_distribution);
					text_rel_diff_scat.setText("" + parameter.rel_diff_scat);

				}
			}
		});

	}

	public void write_raytracer_parameter() // writes the settings of the form
											// to the Parameters.txt file
	{

		try {

			File datei = new File("Parameter_Raytracer.txt");
			FileWriter ausgabestrom = new FileWriter(datei);
			PrintWriter ausgabe = new PrintWriter(ausgabestrom);

			ausgabe.println("#Parameters");
			ausgabe.println("Number_of_threads" + "  = " + text_number_threads.getText());
			ausgabe.println("Number_of_reflections" + "  = " + text_number_reflections.getText());
			ausgabe.println("Number_of_beams" + "  = " + text_number_beams.getText());
			ausgabe.println("#Polarisations: linear_x or linear_y or zirkular or radial or azimutal");
			ausgabe.println("polarization" + "  = " + text_Polarisation.getText());
			ausgabe.println("Number_of_sublevels" + "  = " + text_number_sublevel.getText());
			ausgabe.println("Beam_waist_radius" + "  = " + text_beam_waist_radius.getText());
			ausgabe.println("Beam_wavelength" + "  = " + text_beam_wavelength.getText());
			ausgabe.println("Beam_M2" + "  = " + text_beam_m_square.getText());
			ausgabe.println("Beam_distance_to_focus" + "  = " + text_beam_distance_to_focus.getText());
			ausgabe.println("#Geometry of beam source: circle or rectangle");
			ausgabe.println("Beam_source_type" + "  = " + text_beam_source_typ.getText());
			ausgabe.println("#If beam source type: rectangle");
			ausgabe.println("Beam_source_width_x" + "  = " + text_beam_source_width_x.getText());
			ausgabe.println("Beam_source_width_y" + "  = " + text_beam_source_width_y.getText());
			ausgabe.println("#If beam source type: circle");
			ausgabe.println("Beam_source_radius" + "  = " + text_beam_source_radius.getText());
			ausgabe.println("Beam_min_distance" + "  = " + text_beam_min_dist.getText());
			ausgabe.println("Material_name" + "  = " + text_beam_material_name.getText());
			ausgabe.println("Material_IOR_real" + "  = " + text_beam_IOR_real.getText());
			ausgabe.println("Material_IOR_imag" + "  = " + text_beam_IOR_imag.getText());
			ausgabe.println("Target_csv" + "  = " + text_write_csv.getText());
			ausgabe.println("Target_ply_lin" + "  = " + text_write_ply_lin.getText());
			ausgabe.println("Target_ply_log" + "  = " + text_write_ply_log.getText());
			ausgabe.println("bounding_box_spacing" + "  = " + text_bounding_box_spacing.getText());
			ausgabe.println("Incoming_beam_position_vector" + "  = " + text_incoming_pos_vector.getText());
			ausgabe.println("#Beam direction vector must be parallel to z-axis (other direction not yet implemented, coming soon.)");
			ausgabe.println("Incoming_beam_direction_vector" + "  = " + text_incoming_direction_vector.getText());
			ausgabe.println("#Beam intensity distributions: ring or  gauss or tophat or planewave");
			ausgabe.println("Beam_intensity_distribution" + "  = " + text_beam_intensity_dist.getText());
			ausgabe.println("beamPower" + "  = " + text_beamPower.getText());
			ausgabe.println("rel_diff_scat" + "  = " + text_rel_diff_scat.getText());

			txtOut.setText("Saved");

			ausgabestrom.close();
			ausgabe.close();

		}
		catch (Exception e) {
			System.out.println("Fehler:" + e.toString());
		}

	}
	
	
	public static String get_ProcParams(String fname, String fparam_name) { // STL, labbook
		//String fname = "cut_Cu_20.tif_intens_p4_f_rot.stl";
		String id = fname.split("\\.")[0].replace("cut_", "").substring(1);
		System.out.println(id);
		File fparam = new File(fparam_name);
			try (BufferedReader br = new BufferedReader(new FileReader(fparam))) { // ASCII STL öffnen
				String line;
				String n = "0";
				String k = "0";
				String beamPower = "0";
				String df = "680";
				String sr = "900";
				String dx = "340";
				String nrefl = "12";
				
				String pname; // name of process
				String ppower; // power in exp
				String pdf; // df in exp
				
				line = br.readLine();
				while (((line = br.readLine()) != null) && (!line.trim().equals("END"))) { // Für jede line
					// System.out.println(line);
					line = line.trim();
					System.out.println(line);
					String[] parts = line.split("\t");
					for (int i =0;i<parts.length;i++) {System.out.print(parts[i] + ", ");}
					System.out.println("");
					if (parts.length >2) {
						pname = parts[0];
						ppower = parts[1];
						pdf = parts[2];
						df = Double.toString((Float.parseFloat(pdf)/2.0));
						dx = Double.toString((Float.parseFloat(pdf)/2.0));
						sr = Double.toString((Float.parseFloat(pdf)/2.0*1.)); // fakt. 1 bei tophat
						System.out.println("id: " + id+ "   pname: "+ pname);
						if (pname.contains(id)) { // in der Zeile stehen die Parameter
							
							System.out.println("IDIDIDIDI");
							beamPower = ppower;
							if (pname.contains("82_") || pname.contains("Mg3_")) { // alu
								n = "4";
								k = "10";
								nrefl = "25";
							}
							if (pname.contains("S235_")) { // aSt235
								n = "3.6";
								k = "5";
								nrefl = "12";
							}
							if (pname.contains("Cu_")) { // Cu
								n = "0.54";
								k = "6.74";
								nrefl = "25";
							}
						System.out.println(line);
						System.out.println(pname + ", n=" +n + ", k=" + k + ", nrefl=" + nrefl + ", power=" + beamPower + ", df=" + df + ", dx=" + dx);
						String Pstring = pname + "; Material_IOR_real  = " +n + "; Material_IOR_imag  = " + k + "; Number_of_reflections  = " + nrefl + "; beamPower  = " + beamPower + "; Beam_waist_radius  = " + df + "; Beam_source_radius  = " + sr + "; Incoming_beam_position_vector  = " + dx + ", 0.0, -50.0";
					    return(Pstring);
						}
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
			return "";
		}
	
	public static String get_ProcParams_jonas(String fname, String fparam_name) { // STL file, labbook
		//String fname = "cut_Cu_20.tif_intens_p4_f_rot.stl";
		String id = fname.split("_")[2];
		System.out.println(id);
		File fparam = new File(fparam_name);
			try (BufferedReader br = new BufferedReader(new FileReader(fparam))) { // ASCII STL öffnen
				String line;
				String n = "4";
				String k = "10";
				String beamPower = "0";
				String df = "680";
				String sr = "900";
				String dx = "340";
				String nrefl = "25";
				
				String pname; // name of process
				String ppower; // power in exp
				String pdf; // df in exp
				
				line = br.readLine();
				while (((line = br.readLine()) != null) && (!line.trim().equals("END"))) { // Für jede line
					// System.out.println(line);
					line = line.trim();
					System.out.println(line);
					String[] parts = line.split("\t");
					for (int i =0;i<parts.length;i++) {System.out.print(parts[i] + ", ");}
					System.out.println("");
					if (parts.length >2) {
						pname = parts[0];
						ppower = parts[1];
						pdf = parts[2];
						df = Double.toString((Float.parseFloat(pdf)/2.0));
						dx = Double.toString((Float.parseFloat(pdf)/2.0));
						sr = Double.toString((Float.parseFloat(pdf)/2.0*1.)); // fakt. 1 bei tophat
						System.out.println("id: " + id+ "   pname: "+ pname);
						if (pname.contains(id)) { // in der Zeile stehen die Parameter
							
							System.out.println("IDIDIDIDI");
							beamPower = ppower;
							System.out.println(line);
							System.out.println(pname + ", n=" +n + ", k=" + k + ", nrefl=" + nrefl + ", power=" + beamPower + ", df=" + df + ", dx=" + dx);
							String Pstring = pname + "; Material_IOR_real  = " +n + "; Material_IOR_imag  = " + k + "; Number_of_reflections  = " + nrefl + "; beamPower  = " + beamPower + "; Beam_waist_radius  = " + df + "; Beam_source_radius  = " + sr + "; Incoming_beam_position_vector  = " + dx + ", 0.0, -50.0";
							return(Pstring);
						}
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
			return "";
		}

	public void adapt_Ray_Params(String fname, String fparam, String rayname) { // STL file, Labbook, ...
		File rayfile = new File(rayname);// "C:\\Users\\ifsw-ffetzer\\Desktop\\SVN\\TheRaytracer\\Bin\\Parameter_Raytracer.txt"
		//String params = get_ProcParams(fname, fparam); // STL, labbook
		String params = get_ProcParams_jonas(fname, fparam); // STL, labbook --> Jonas
		System.out.println(params);
		String[] para = params.split("; ");
		for (int i =1;i<para.length;i++) {
			System.out.println(para[i]);
			Algorithmen.set_Raytracer_Param(rayfile, para[i]);
		}
		
	}
	
	public void test_func() {

		int st_depth = Math.max(1, (int) Math.ceil(Math.log(0.0) / Math.log(2)) - 1);
		//if (st_depth >20) {
			//st_depth = 16;
		//}
		System.out.println("Setting " + st_depth + " as ST_depth");
		System.out.println(String.valueOf(String.valueOf(st_depth)));
		/*
		String fname = "AVG_SK_032_C1.tif_intens_p4_f_rot.stl";
		String fparam = "G:\\jonas_durchschweissen\\Rek\\STL\\paras.txt";
		String rayname = "C:\\Users\\ifsw-ffetzer\\Desktop\\SVN\\TheRaytracer\\Bin\\Parameter_Raytracer.txt";
		adapt_Ray_Params(fname, fparam, rayname);
		*/
	}
	
	}

