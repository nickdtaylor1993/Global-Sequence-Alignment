import java.io.*; 
import java.util.*;
import java.lang.*; 
import Jama.*; 
//Author: Nicholas Taylor

/*
 * Class used for parsing the multiple files for the appropriate inputs for the sequence comparison algorithm. 
 */
public class Parser{
	private static Map<String,String> table = new HashMap<String,String>(); //Hash map for the appropriate values in the config file.  
	private static final int NUM_GAP_PENALTIES = 4; //Constant for the total number of gap penalties for the DNA sequences.  
	private static char[] sequenceA; //Char array of the sequence A. 
	private static char[] sequenceB; //Char array of the sequence B. 

	/*
	 * Helper method for making a HashMap table of possible options the config file can take. 
	 */
	private static void initTable(){
		table.put("sequenceInputFile","");
		table.put("sequenceA",""); 
		table.put("sequenceB",""); 
		table.put("weightMatrixFile",""); 
		table.put("gapPenaltyFile","");
		table.put("maxColumns",""); 
		table.put("writeToFile","");
		table.put("outputFile","");
		table.put("debugMode",""); 
	}
			
	/*
	 * Method used to scan the Configuration file properties and sets them into the System properties. 
	 */
	public static void parseConfigFile(File file) throws IOException{
		initTable(); //inits the hashmap table. 

		boolean writeToFile = false; //Used to determine if config file should be written to a file or not. 

		int lineCount = 0; //Keeps track of the line numbers in the config file.

		System.setProperty("maxColumns","80"); //Sets the default maxColumns to 80. 
		
		Scanner sc = new Scanner(file); //Inits the scanner with the conf file.  
			
		while(sc.hasNextLine()){
			lineCount++; 

			String token = sc.nextLine().replaceAll("\\s","");  //Removes all the white space from the string. 
				
			if(token.length() < 1 || token.charAt(0) == '#') //Checks if the line is blank or is a comment.  
				continue; 

			int endIndex = token.indexOf('='); //Gets the end index of the config key option value.
		
			if(endIndex == -1)  //If there is no = sign means invalid line is in the file throw exception clause. 
				throw new IOException("Error invalid line at: " + lineCount);

			String keyToken = token.substring(0,endIndex); //Gets the configuration key option value.
			
			if(table.containsKey(keyToken)){ //Checks if the line in the config file key value is appropriate by looking in the hash table. 
				if(keyToken.equals("writeToFile")){ //Checks the writeToFile option if it is true sets the boolean writeToFile as true and vice versa.  
					String valueToken = token.substring(endIndex+1);
					if(valueToken.equals("true"))
					  writeToFile = true;
					else
					 writeToFile = false; 
				}
				
				String valueToken = token.substring(endIndex+1); //Gets the configuration value for the key option value.  
				
				if(writeToFile == true && keyToken.equals("outputFile")){ //If writeToFile is true checks the outputFile value to confirm if its appropriate or not. 
					if(valueToken.length() < 1)
						throw new IOException("Error Invalid line at: " + lineCount); //No output file even though was specified. 
					System.setProperty(keyToken, valueToken); //Sets the config option key with its value from the config with System's property. 
				}

				else if(keyToken.equals("outputFile")){	 //If writeToFile is false no needed valueToken needed.  
					System.setProperty(keyToken, valueToken); //Sets the config option key with its value from the config with System's property. 
				}

				else{

					if(valueToken.length() < 1)
						throw new IOException("Error incomplete line of option at: " + lineCount); 

					System.setProperty(keyToken, valueToken); //Sets the config option key with its value from the config with System's property. 
				}
			}
			else
			   throw new IOException("Error invalid line at: " + lineCount); //If the keyToken is not in the table is a invalid option.  
		}
	}

	/*
	 * Method used to parse the Weight Matrix file and returns the Matrix representing the file's inputs. 
	 */
	public static Matrix parseWeightMatrixFile(String str)throws IOException{  
		Scanner sc = new Scanner(new File(str)); //Creates a scanner object to scan the Weight Matrix file. 
	
		List<Double> list = new ArrayList<Double>(); //Used to store all the values of the Weight Matrix file to enter into the Matrix. 
		
		int row = 0; //Row and Col used for defining the size of the array before creating an Matrix object.
		int col = 0; 
		
		String line = ""; //Default empty string before scanning the file. 
		
		while(sc.hasNextLine()){ //Scanning the Weight Matrix file. 
			line = sc.nextLine();  //Grabs the line of weights doubles. 

			if(line.length() < 1 || line.charAt(0) == '#') //Checks if the line is empty or is a comment. 
				continue; 

			row++; //increments the number of rows for the Weight Matrix. 
			
			String[] arr = line.split("\\s");  //Splits the weight matrix at the particular row into an array of double Strings. 
			
			col = arr.length; //Gets the array length which is the number of weights in the line of the file. 
			for(int i = 0; i < arr.length; i++){ //converts all the values in the current line into doubles and puts them into a list. 
				list.add(Double.valueOf(arr[i])); 
			}
		}

		double[][] array = new double[row][col]; //Inits an array with the row and col figured after scanning the file.  
		int k = 0;
		for(int i = 0; i < row; i++){ //Iterates and puts into the array all of the values from the file into the 2d Array. 
			for(int j = 0; j < col; j++){
				array[i][j] = (double) list.get(k);
				k++; 
			}
		}
		return new Matrix(array); 
	}
	
	/*
	 * Method parses the gap penalty file with each gap penalty associated with A,C,G,and T respectively returned with an double single array. 
	 */
	public static double[] parseGapPenaltyFile(String str)throws IOException {
		Scanner sc = new Scanner(new File(str)); //Creates a scanner object to scan the Gap Penalty file.
		
		double[] array = new double[NUM_GAP_PENALTIES] ; //Adds each of the items in the list into a double array. 
		
		String line = ""; 
		
		while(sc.hasNextLine()){ 
			line = sc.nextLine(); 
			
			if(line.length() < 1 || line.charAt(0) == '#') //Checks if the line is empty or is a comment. 
				continue; 

			String[] penalties = line.split("\\s");  //Splits the weight matrix at the particular row into an array of double Strings. 
			
			for(int i = 0; i < array.length; i++){ //converts all the values in the current line into doubles and puts them into a list. 
				array[i] = (double) Double.valueOf(penalties[i]); 
			}
		}
		return array; 
	}

	/*
	 * Method used to parse the Sequence Input file and initializes the two sequences A and B. 
	 */
	public static void parseSequenceInputFile(String str)throws IOException{
		Scanner sc = new Scanner(new File(str)); //Creates a scanner object to scan the Gap Penalty file.
		
		String line = ""; 
		while(sc.hasNextLine()){ //Scans the sequence file. 
			line = sc.nextLine().toUpperCase(); 
			
			if(line.length() < 1 || line.charAt(0) == '#') //Checks if the line is empty or is a comment. 
			  continue;
			int isSequenceA = line.indexOf(System.getProperty("sequenceA")); //Use to determine if the sequence is for sequenceA 
				
			int isSequenceB = line.indexOf(System.getProperty("sequenceB")); //Use to determine if the sequence is for sequenceB
			
			if(isSequenceA != -1){ //If the current line is the sequenceA proceed. 
				line = line.substring(line.lastIndexOf(":")+1).toUpperCase(); //Gets the entire line from the ":" index. 
				while(sc.hasNextLine()){ //Used to scan the next lines if the entire sequence is not on the entire line. 
			 	   if(line.indexOf("X") != -1) //Checks if the line has the terminator for the end of the Sequence if so breaks. 
					  break; 
				line += sc.nextLine().toUpperCase(); //Otherwise sequence must break to the next line so add the next line into the sequence. 
			  }
			  	line = line.substring(0, line.indexOf("X")); //Grabs the entire sequence till the Terminator
			
				sequenceA = line.toCharArray(); //converts the string into a array of Chars'. 
			}

			else if(isSequenceB != -1){ //Checks the current to see if is sequenceB then proceeds. 
			 	 line = line.substring(line.lastIndexOf(":")+1); //Grabs the entire line from the ":" index. 
			 
			 	 while(sc.hasNext()){ 
			 		if(line.indexOf("X") != -1) //Checks if the line has termiantor for the end of the sequence. 
					  break; 
				
				line += sc.nextLine().toUpperCase(); //If no terminator add the next line to the sequence and repeat. 
			 }
			 	line = line.substring(0, line.indexOf("X")); //Grabs the entire sequence until the terminator indicator.  
				
				sequenceB = line.toCharArray(); //Converts the string into a char array 
			}

			else
			  continue;
		}
		if(sequenceA == null || sequenceB == null)
			throw new IOException("Invalid sequence number specification"); 
	}
	
	/*
	 * Returns sequence A as char array. 
	 */
	public static char[] getSequenceA(){
		return sequenceA; 
	}
	
	/*
	 * Returns sequence B as char array. 
	 */
	public static char[] getSequenceB(){
		return sequenceB; 
	}
}
