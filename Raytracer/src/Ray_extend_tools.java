import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Scanner;

public class Ray_extend_tools {
    public static void write_array_to_file(String fname, double[][] arr) {
    	FileWriter writer;
		try {
			writer = new FileWriter(fname); 
			for(int i =0;i< arr.length;i++) {
				for(int j = 0;j<arr[0].length;j++) {
					writer.write(String.valueOf(arr[i][j]) + "\t");
					
				}
				writer.write(System.lineSeparator());
			}
			writer.close();   
		} catch (IOException e) {
		// TODO Auto-generated catch block
		e.printStackTrace();
		}
    }
    
    public static double[][] read_array_from_file(String fname){ // 2 d array, columns x, rows y returns list of arrays, x, y, z
    	Scanner sc;
    	double [][] myArray = null;
    	//System.out.println("1"); 
		try {
			sc = new Scanner(new BufferedReader(new FileReader(fname)));
			//System.out.println("2"); 
			String[] linelist = sc.nextLine().trim().split("\t");
			sc.nextLine();
			int rows = 2;
			int columns = linelist.length+1;
			//System.out.println(columns);
			while(sc.hasNextLine()) {
				sc.nextLine();
				rows = rows+1;
				//System.out.println(rows); 
			}
			//System.out.println(rows); 
			sc.close();
			myArray = new double[rows][columns];
			// open file for read
			try {
			File file=new File(fname);
	        BufferedReader br=new BufferedReader(new FileReader(file));
	        String line;
	        line=br.readLine();
			// first line
			linelist = line.trim().split("\t");
			//System.out.print("first : " + linelist[0] + "\n");
			myArray[0][0] = -999;
            for (int j=0; j<linelist.length; j++) {
               myArray[0][j+1] = Double.parseDouble(linelist[j]);
            }
            int i = 1;
            while((line=br.readLine())!=null){
	            linelist = line.trim().split("\t");
	            //System.out.print(linelist[0]);
	            for (int j=0; j<linelist.length; j++) {
	               myArray[i][j] = Double.parseDouble(linelist[j]);
	            }
            	i = i+1;
	        }
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
	        return myArray;
		} 
		catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		 return myArray;
	    }
    
    public static double lin_interpolate2d(double x, double y, double[][] z) { // erster index :x
    	double iz = 0;
    	double[] xx = new double[z.length-1];
    	//System.out.println(xx.length);
    	double[] yy = new double[z[0].length-1];
    	//System.out.println(yy.length);
    	for (int ii = 0; ii < xx.length;ii++) {
    		xx[ii]=z[ii+1][0];
    	}
    	for (int jj = 0; jj < yy.length;jj++) {
    		yy[jj]=z[0][jj];
    	}
    	// search in x,y for neighbour
    	int i = 0;
    	int j = 0;
    	while (x>xx[i] && i<xx.length-i) {i++;}
    	while (y>yy[j] && j< yy.length-1) {j++;}
    	if (i>0 && j>0) {
    		i = i+1; j = j+1;
    		//System.out.println("i: "+ i + ", j: " + j);
    		double i1 = z[i-1][j-1] + (z[i-1][j]-z[i-1][j-1])/(yy[i]-yy[i-1])*(y-yy[i-1]);
    		double i2 = i1 + (z[i][j]-z[i-1][j])/(xx[i]-xx[i-1])*(x-xx[i-1]);
    		
    		double i3 = z[i-1][j-1] + (z[i][j-1]-z[i-1][j-1])/(xx[i]-xx[i-1])*(x-xx[i-1]);
    		double i4 = i3 + (z[i][j]-z[i][j-1])/(yy[i]-yy[i-1])*(y-yy[i-1]);
    		iz = (i2+i4)/2.0;
    		//iz = (z[i-1][j-1] + z[i][j] +z[i-1][j]+z[i][j-1])/4.;
    	}
    	else {
    		iz = 0;
    	}
    	return iz;
    }
}
