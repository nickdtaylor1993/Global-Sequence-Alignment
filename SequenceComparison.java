import java.io.*; 
import java.util.*;
import java.lang.*; 
import Jama.*; 
//Author: Nicholas Taylor 

/*
 * Class used for generating the optimal alignment between two sequences by performing the Smith-Waterman Global Sequence Comparison algorithm. 
 */
public class SequenceComparison{
	private Matrix weightMatrix; //The weight matrix for determining the optimal alignment between the sequences.  
	private char[] sequenceA; //SequenceA which represents a DNA sequence. 
	private char[] sequenceB; //SequenceB which represents a DNA sequence. 
	private double[] gapPenalties; //gapPenalties for each particular nucleotide A,C,G and T. 
	
	private String optimalAlignmentA = ""; //The string representation of the optimalAlignment of sequenceA with respect to 
	private String optimalAlignmentB = "";  //The other optimal alignment of sequenceB. 
	private double optimalSimilarityScore; //Gets the similarityScore between the two optimal Alignments of sequences  A and B.  
	private double[][] sequenceAnalysisArray;
	private char[][] directionOfPreviousMax; 
	private static final int A = 0; //Each of these nucleotide DNA are used to index into the gapPenality as well as the Weight Matrix. 
	private static final int C = 1; 
	private static final int G = 2; 
	private static final int T = 3; 
	
	/*
	 * Constuctor for init the SequenceComparison instance with its weightMatrix, gapPenalties, and two sequences. 
	 */
	public SequenceComparison(Matrix weightMatrix, double[] gapPenalties,
		char[] sequenceOne, char[] sequenceTwo){ 
		
		this.weightMatrix = weightMatrix; //WeightMatrix. 
		
		this.gapPenalties = new double[gapPenalties.length]; //gapPenalties. 

		for(int i = 0; i < gapPenalties.length; i++){ //Deep copy of the gap Penalities. 
			this.gapPenalties[i] = gapPenalties[i]; 
		}
		
		this.sequenceA = new char[sequenceOne.length]; 
		
		for(int j = 0; j < sequenceOne.length; j++){ //Deep copy of the sequenceA. 
			this.sequenceA[j] = sequenceOne[j];
		}
		
		this.sequenceB = new char[sequenceTwo.length]; 
		for(int k = 0; k < sequenceTwo.length; k++){ //Deep copy of the sequenceB. 
			this.sequenceB[k] = sequenceTwo[k]; 
		}
		this.sequenceAnalysisArray = new double[this.sequenceA.length+1][this.sequenceB.length+1]; //Creates the sequenceAnalysisArray.
		
		
		this.directionOfPreviousMax = new char[sequenceAnalysisArray.length][sequenceAnalysisArray[0].length]; //Used to indicate direction of the previous max direction with 'D' being diagonal, 'U' being up, and 'L' being left.
	}
	
	/*
	 * Helper function that initializes the base cases for the Sequences set such as S(0,0) = 0, S(i,0) = g(uk), and S(0,j) = g(vk). 
	 */
	 private void initBaseCases(){
		sequenceAnalysisArray[0][0] = 0; //Sets the S(0,0) as 0
		double gap = 0; //Used as the gap.  
		
		for(int i = 0; i < this.sequenceA.length; i++){ //Sets the S(i,0) of the sequenceAnalaysisArray. 		
			switch(sequenceA[i]){ //Using a switch case for each nucleotide case then adding the gap penalty for that particular nucleotide. 
				case 'A': gap+=gapPenalties[this.A];
					   break; 
				case 'C': gap+=gapPenalties[this.C]; 
					  break; 
				case 'G': gap+=gapPenalties[this.G]; 
					  break;
				case 'T': gap+=gapPenalties[this.T]; 
					  break; 
				default:  
					  break; 
			 }
			 sequenceAnalysisArray[i+1][0] = gap;
			 directionOfPreviousMax[i+1][0] = 'U'; 
		}
		gap = 0; //Reset the gap to 0 for doing the next sequence. 
		
		for(int j = 0; j < sequenceB.length; j++){ //Sets the S(0,j) of the sequenceAnalysisArray.
			switch(sequenceB[j]){
				case 'A': gap+=gapPenalties[this.A];
					  break;
				case 'C': gap+=gapPenalties[this.C]; 
					  break;
				case 'G': gap+=gapPenalties[this.G]; 
					  break;
			 	case 'T': gap+=gapPenalties[this.T]; 
					  break; 
				default:
					  break; 
			}
			sequenceAnalysisArray[0][j+1] = gap; 
			directionOfPreviousMax[0][j+1] = 'B'; 
		}
	 }
	
	/*
	 * Helper method for calculating the Gap Penalty of the particular DNA Nucleotide. 
	 */
	private double getGapPenaltyValue(char c){
		double gap = 0; 
		switch(c){ //Using a switch statement gets the appropriate gap penalty with the DNA Nucleotide.  
			case 'A': gap=gapPenalties[this.A];
				   break;
			case 'C': gap=gapPenalties[this.C]; 
				   break;
			case 'G': gap=gapPenalties[this.G]; 
				   break;
			case 'T': gap=gapPenalties[this.T]; 
				   break; 
			default:
				   break; 
			}
		return gap;  //Returns the gap penalty. 
	}

	/*
	 * Helper method to calculate the similarity score between two DNA nucleotides with the weight matrix. 
	 */
	private double similarityScore(char ui, char vi){
		int i = 0; //Holders for the value of A,C,G,or T
		int j = 0; //same as above. 

		switch(ui){ //Using a switch statement gets the appropriate DNA Nucleotide.  
			case 'A': i=this.A;
				   break;
			case 'C': i=this.C; 
				   break;
			case 'G': i=this.G; 
				   break;
			case 'T': i=this.T; 
				   break; 
			default:
				   break; 
		}
		
		switch(vi){ //Using a switch statement gets the appropriate DNA Nucleotide.  
			case 'A': j=this.A;
				   break;
			case 'C': j=this.C; 
				   break;
			case 'G': j=this.G; 
				   break;
			case 'T': j=this.T; 
				   break; 
			default:
				   break; 
		}
		return weightMatrix.get(i,j); //Returns the weight of the two DNA nucleotides.  
	}
	

	/*
	 * Helper method to get the maximum 
	 */
	private char getMax(double up, double diagonal, double back){
		double epsilon = .01;

		if((up >= diagonal) && (up >= back)){ //If all three are equal up has highest priority. 
			return 'U'; 
		}

		else if((diagonal >= up) && (diagonal >= back)){ //If up and diagonal equal for max then up is max. 
			return 'D'; 
		}
		
		else if((back >= diagonal) && (back >= up)){ //If up and back are equal and up is bigger than diagonal, up is max. 
			return 'B'; 
		}

		else{ //Otherwise if diagonal and up is less than Back is max. 
			System.out.println("Something bad happend");
			System.exit(1);
			return 0; 
		}
	}

	/*
	 * Method computes the Smith-Waterman Global Sequence Comparison Algorithm. 
	 */
	public void computeOptimalAlignments(){
		initBaseCases();  //Inits the sequenceAnalaysisArray with the base settings. 
		
		double up = 0; //Used for representing the value from the left. 
		
		double diagonal = 0; //Used for representing the value from the diagonal. 
		
		double back = 0; //Used for representing the value of the left. 

		for(int i = 1; i < sequenceA.length+1; i++){ //Iterates through the sequenceAnalysisArray till we reach the furthest bottom right corner. 
			for(int j = 1; j < sequenceB.length+1; j++){
				up = sequenceAnalysisArray[i-1][j] + getGapPenaltyValue(sequenceA[i-1]); //Gets the up value.
			
				diagonal = sequenceAnalysisArray[i-1][j-1] + similarityScore(sequenceA[i-1],sequenceB[j-1]); //Gets the diagonal value. i			   
				back = sequenceAnalysisArray[i][j-1] + getGapPenaltyValue(sequenceB[j-1]); //Gets the back value. 
				
				directionOfPreviousMax[i][j] = getMax(up,diagonal,back); //returns which direction was the max into the char array indicating the direction of the previous max.
				
				switch(directionOfPreviousMax[i][j]){ //Inserts the max value of either the Up, Diagonal, and Back into the current position of the sequenceAnalayisArray. 
					case 'U': sequenceAnalysisArray[i][j] = up; 
						  break; 
				    	case 'D': sequenceAnalysisArray[i][j] = diagonal; 
						  break; 
					case 'B': sequenceAnalysisArray[i][j] = back; 
						  break; 
					default: 
						  break; 
				}
			}
		}
		
		int i = sequenceAnalysisArray.length-1; 
		
		int j = sequenceAnalysisArray[0].length-1; 
		optimalSimilarityScore = sequenceAnalysisArray[i][j]; //(sequenceA[i-1], sequenceB[j-1]); //Gets the optimal similarity score between the two optimal alignments of sequences A and B. 
		while(i+j != 0){ //Tracing the optimal Alignment back to the origin of the sequenceAnalysisArray. 
			switch(directionOfPreviousMax[i][j]){ //Following the previous max
				case 'U': optimalAlignmentA = sequenceA[i-1] + optimalAlignmentA; //Case if the previous max direction was up we set the sequenceA to its nucleoTide but make a gap in sequenceB
					  optimalAlignmentB = "-" + optimalAlignmentB; 
					  i--; //We then move up. 
				          break; 
				case 'D': optimalAlignmentA = sequenceA[i-1] + optimalAlignmentA; //Case if the previous max direction was diagonal we set both sequenceA and sequenceB to their nucleoTides. 
					  optimalAlignmentB = sequenceB[j-1] + optimalAlignmentB; 
					  i--;
					  j--; //Then move diagonal and continue till we reach s(0,0).     
					  break;

				case 'B': optimalAlignmentA = "-" + optimalAlignmentA; //Case where the previous max direction was back we set just the sequenceA as a gap and set the sequenceB as its nucleotide DNA sequence.  
					  optimalAlignmentB = sequenceB[j-1] + optimalAlignmentB; 
					  j--; //Move back and continue. 
					  break; 
				default: 
					  break; 
				}
		}
	}

	
	/*
	 * Method either prints to standard output or to a given file. 
	 */
	public void displayOutput()throws IOException{
			String cols = System.getProperty("maxColumns"); //Gets the max column. 
			int num_of_cols = (int) Integer.valueOf(cols); //same as above. 
			
			String lineA = ""; 
			String lineB = "";
			
			if(System.getProperty("writeToFile").equals("true")){ //For the case of writing to a file. 
				String absolutePath = System.getProperty("user.dir"); //Get the absolutePath of the directory to write the file. 
				
				File fileToWrite = new File(absolutePath +File.separator+System.getProperty("outputFile")); //Create the file to write. 
			
				if(!fileToWrite.exists()) //If the file does not yet exist we create it. 
					fileToWrite.createNewFile(); 
				
				PrintWriter pw = new PrintWriter(fileToWrite);  //Create instance of a PrintWriter object from the File object to write into the file. 
				
				pw.print("Optimal Similarity Score: " + optimalSimilarityScore); //Write the optimal similarity score. 
				pw.println(); 
				pw.print("Optimal Sequence Alignments A and B");
				pw.println();
				
				int count = 0; 
				for(int i = 0; i < optimalAlignmentA.length(); i++){ //Write the Optimal Alignment Sequence for Sequence A with the desired maximum Columns. 
					if(count < num_of_cols){
						lineA +=optimalAlignmentA.charAt(i); 
						lineB +=optimalAlignmentB.charAt(i); 
					}
					else{
							pw.println(lineA);
							pw.println(lineB);
							pw.println(); 
							count = 0; 
							lineA = ""; 
							lineB = ""; 
							lineA +=optimalAlignmentA.charAt(i); 
							lineB +=optimalAlignmentB.charAt(i); 

					}
					count++; 
				}

				if(lineA.length() >=1 && lineB.length() >=1){
						pw.println(lineA); 
						pw.println(lineB);
				}

				pw.close(); 
			}

			else{ //The case we do not specify to write to a file and instead which to print to standard out. 	
				System.out.println("Optimal Similarity Score: " + optimalSimilarityScore);
				System.out.println("Optimal Sequence Alignments A and B");
				int count = 0; 
				for(int i = 0; i < optimalAlignmentA.length(); i++){
					if(count < num_of_cols){
						lineA +=optimalAlignmentA.charAt(i);
						lineB +=optimalAlignmentB.charAt(i);
					}

					else{
						System.out.println(lineA);
						System.out.println(lineB); 
						System.out.println(); 
						count = 0;
						lineA = ""; 
						lineB = "";
						lineA +=optimalAlignmentA.charAt(i);
						lineB +=optimalAlignmentB.charAt(i);

					}
					count++; 
				}
				
				if(lineA.length() >=1 && lineB.length() >=1){
						System.out.println(lineA); 
						System.out.println(lineB);
				}
			}
	}

	public static void main(String[] args){
		try{
			if(args.length < 1){ //Checks if includes a command line argument for the configuration file. 
				System.out.println("Error: must include configuration file"); 
				System.exit(0); 
			}
		
			Parser.parseConfigFile(new File(args[0])); //Parses the config file. 
			
			Matrix weightMatrix = Parser.parseWeightMatrixFile(System.getProperty("weightMatrixFile")); //parses the Weight Matrix file and obtains the actual weightMatrix. 
			
			double[] penalties = Parser.parseGapPenaltyFile(System.getProperty("gapPenaltyFile")); //parses the Gap Penalty file and creates a double array of the penalties from the file. 
			
			Parser.parseSequenceInputFile(System.getProperty("sequenceInputFile")); //parses the sequenceInput files to obtain each sequence A and B. 
			SequenceComparison sc = new SequenceComparison(weightMatrix,penalties,Parser.getSequenceA(),Parser.getSequenceB()); //Creates the SequenceComparison object. 
			sc.computeOptimalAlignments(); //Performs the Smith Watermon Global Alignment algorithm. 
			
			sc.displayOutput(); //Displays the output to either standard out or to specified output file.  

		}
		catch(IOException e){
			System.out.println(e.getMessage());
			System.exit(0); 
		}
	}
}
