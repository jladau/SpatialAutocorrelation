package edu.ucsf.SpatialAutocorrelation;

import java.util.ArrayList;

import com.google.common.base.Strings;
import com.google.common.collect.Range;

import edu.ucsf.base.ClusterIterator;
import edu.ucsf.base.Permutation;
import edu.ucsf.geospatial.SpatialAutocorrelation;
import edu.ucsf.geospatial.SpatialWeightsMatrix;
import edu.ucsf.io.ArgumentIO;
import edu.ucsf.io.BiomIO;
import edu.ucsf.io.DataIO;

/**
 * This class implements spatial autocorrelation analyses for OTU tables.
 * @author Joshua Ladau, jladau@gmail.com
 */
public class SpatialAutocorrelationLauncher {

	/**
	 * Runs spatial autocorrelation analysis. Writes table with Moran's I and significance values for each OTU to specified output path. Can also optionally write a table of data that was analyzed and temporal and spatial distances between pairs of samples.
	 * @param rgsArgs Arguments. Arguments should be passed as --{argument name}={argument value}. Name-value pairs are:
	 * 				
	 * 				<p>
	 * 				<h4 class="list-heading">Required arguments</h4>
	 * 				<ul> 
	 * 				<li>sDataPath [string] = Absolute path to a file containing an OTU table. File should be in BIOM format (HDF5).
	 * 				<p>
	 * 				<li>sOutputPath [string] = Absolute path for output file. Should have a "csv" suffix.
	 * 				<p>
	 * 				<li>iMCMCChains [integer] = Number of independent MCMC chains to run for calculating p-values.
	 * 				<p>
	 * 				<li>iMCMCIterations [integer] = Total number of MCMC iterations. Each chain will have iMCMCIterations/iMCMCChains iterations.
	 *				<p>
	 * 				<li>rgsDistanceNeighborhoods [list of strings] = List of distance neighborhoods. Minimum and maximum distances should be separated by dashes; for instance 1-10,10-100,100-1000. If locations are given by latitude and longitude, then distances should be in kilometers. If locations are given by arbitrary x,y coordinates, then the distance units should be those of the coordinate system.
	 *				</ul>
	 *
	 *				<p>				
	 *				<h4 class="list-heading">Optional arguments: autocorrelation</h4>
	 * 				<ul>				
	 *				<li>rgsDirectionNeighborhoods [list of strings] [OPTIONAL] = List of direction neighborhoods. Minimum and maximum directions should be separated by dashes; for instance 0-180,180-270,270-360.
	 *				<p>
	 *				<li>rgsTimeNeighborhoods [list of strings] [OPTIONAL] = List of time neighborhoods. Minimum and maximum time differences should be separated by dashes; for instance 1-2,2-5,5-12.
	 *				<p>
	 *				<li>bOutputData [boolean] [OPTIONAL] = Flag for whether to output table of data (after filtering for prevalence, etc) that was analyzed. Defaults to true. 
	 *				<p>
	 *				<li>bOutputDistances [boolean] [OPTIONAL] = Flag for whether to output table time differences and geographic distances for samples. Defaults to false. 				
	 *				<p>
	 * 				<li>sAnalysis [string] [OPTIONAL] = Type of analysis to run. Accepted values are "MoranScatter" for Moran's scatter plots and "MoransI" for Moran's I and significance test. Defaults to MoransI. 
	 *				<p>
	 *				<li>sWeight [string]  [OPTIONAL] = Type of weighting to use in spatial weights matrix. Accepted values are "inverse" or "binary". Defaults to "binary".
	 *				</ul>
	 *
	 *				<p>				
	 *				<h4 class="list-heading">Optional arguments: parallelization</h4>
	 *				<ul>
	 *				<li>iTotalTasks [integer] [OPTIONAL] = Total number of tasks for parallelizing. Omitting this argument results in no parallelization.
	 *				<p>
	 *				<li>iTaskID [integer] [OPTIONAL] = Task to be run if parallelized. Omitting this argument results in no parallelization.
	 * 				</ul>
	 * 
	 * 				<p>
	 *				<h4 class="list-heading">Optional arguments: data table</h4>
	 *				<ul>
	 *				<li>sTaxonRank [string] [OPTIONAL] = Taxonomic units on which to collapse table. Accepted values are "kingdom", "phylum", "class", "order", "family", "genus", "species", or "otu." The value of "otu" will cause table to not be modified.
	 *              <p>
	 *              <li>sSampleMetadataPath [string] [OPTIONAL] = Path to text file containing sample metadata formatted according to http://biom-format.org/documentation/adding_metadata.html. For use if BIOM file does not contain metadata. Must include "id", "datetime" and "latitude", "longitude" or "x","y" fields
	 *				<p>
	 *              <li>sSamplesToKeepPath [string] [OPTIONAL] = Path to file with list of samples to keep. File should contain a list of sample names.
	 *              <p>
	 *              <li>sObservationsToKeepPath [string] [OPTIONAL] = Path to file with list of observations to keep. File should contain a list of observation names.
	 *				<p>
	 *              <li>iRandomSampleSubsetSize [integer] [OPTIONAL] = Number of randomly chosen samples to use. Useful for analyzing large data tables quickly.
	 *              <p>
	 *              <li>iRandomObservationSubsetSize [integer] [OPTIONAL] = Number of randomly chosen observations to use. Useful for analyzing large data tables quickly. 
	 *              <p>
	 *              <li>bCheckRarefied [boolean] [OPTIONAL] = Flag for whether to check for rarefaction. If enabled and table is not rarefied, error will be thrown. Default is true.
	 *              <p>
	 *              <li>bNormalize [boolean] [OPTIONAL] = Flag for whether to normalize within each sample so that entries total to 1. Default is true.
	 *              <p>
	 *              <li>iPrevalenceMinimum [integer] [OPTIONAL] = Minimum prevalence: observations that occur in fewer samples will be omitted from analysis. Default is 10.
	 *              <p>
	 *              <li>bPresenceAbsence [boolean] [OPTIONAL] = Flag for whether data should be reduced to presence-absence data. Default is false.
	 *              <p>
	 *              <li>iRarefactionTotal [integer] [OPTIONAL] = Total count to which to rarefy samples.
	 *				</ul>
	 **/
	public static void main(String rgsArgs[]){
		
		//iCounter = counter
		//sMode = latitude-longitude or euclidean
		//lstOut = output
		//rgp1 = random permutations
		//spw1 = spatial weights matrix
		//arg1 = arguments
		//spa1 = geospatial statistics object
		//cit1 = cluster iterator
		//bio1 = biom object
		
		int iCounter;
		String sMode;
		ArrayList<String> lstOut=null;
		Permutation<Integer> rgp1[];
		SpatialWeightsMatrix spw1 = null;
		ArgumentIO arg1;
		SpatialAutocorrelation spa1;
		ClusterIterator cit1;
		BiomIO bio1;
		
		//loading arguments
		arg1 = new ArgumentIO(rgsArgs);
		arg1.setErrorReporting(true);
		
		//loading defaults
		if(!arg1.containsArgument("sAnalysis")){
			arg1.updateArgument("sAnalysis", "MoransI");
		}
		if(!arg1.containsArgument("sWeight")){
			arg1.updateArgument("sWeight", "binary");
		}
		if(!arg1.containsArgument("iPrevalenceMinimum")){
			arg1.updateArgument("iPrevalenceMinimum", 10);
		}
		if(!arg1.containsArgument("bNormalize")){
			arg1.updateArgument("bNormalize", true);
		}
		if(!arg1.containsArgument("bOutputData")){
			arg1.updateArgument("bOutputData", true);
		}
		if(!arg1.containsArgument("bOutputDistances")){
			arg1.updateArgument("bOutputDistances", false);
		}
		if(!arg1.containsArgument("bPresenceAbsence")){
			arg1.updateArgument("bPresenceAbsence", false);
		}
		if(!arg1.containsArgument("bCheckRarefied")){
			arg1.updateArgument("bCheckRarefied", false);
		}
		if(!arg1.containsArgument("rgsDirectionNeighborhoods")){
			arg1.updateArgument("rgsDirectionNeighborhoods", (Strings.repeat("0-360,",arg1.getValueStringArray("rgsDistanceNeighborhoods").length-1) + "0-360").split(","));
		}
		if(!arg1.containsArgument("rgsTimeNeighborhoods")){
			arg1.updateArgument("rgsTimeNeighborhoods", (Strings.repeat("0-12,",arg1.getValueStringArray("rgsDistanceNeighborhoods").length-1) + "0-12").split(","));
		}
		
		arg1.updateArgument("rgsRequiredSampleMetadata", "datetime");
		
		//loading data
		System.out.println("Loading data...");
		try {
			bio1 = new BiomIO(arg1.getValueString("sDataPath"),arg1.getAllArguments());
		} catch (Exception e) {
			e.printStackTrace();
			return;
		}
		
		//loading mode and filtering samples without metadata
		System.out.println("Filtering samples without location metadata...");
		if(bio1.axsSample.hasMetadataField("latitude") && bio1.axsSample.hasMetadataField("longitude")){
			sMode="latitude-longitude";
			try{
				bio1.filterByNoMetadata(new String[]{"latitude","longitude"}, bio1.axsSample);
			}catch(Exception e){
				e.printStackTrace();
			}
		}else if(bio1.axsSample.hasMetadataField("x") && bio1.axsSample.hasMetadataField("y")){
			sMode="euclidean";
			try{
				bio1.filterByNoMetadata(new String[]{"x","y"}, bio1.axsSample);
			}catch(Exception e){
				e.printStackTrace();
			}
		}else{	
			System.out.println("Error: spatial metadata not found. Exiting.");
			return;
		}
		
		//initializing output
		lstOut = new ArrayList<String>();
		lstOut.add("Taxon" +
				",Morans I" +
				",Randomized p-value" +
				",Z-score" +
				",Distance neighborhood (km)" +
				",Direction neighborhood (degrees)" +
				",Time neighborhood (months)" +
				",P-value iterations" +
				",Number non-zero weights" +
				",Fraction observations with neighbors" +
				",Taxon prevalence (number of samples)" +
				",Taxon mean relative abundance");
		
		//initializing cluster iterator
		cit1 = new ClusterIterator(arg1.getValueInt("iTaskID"),arg1.getValueInt("iTotalTasks"));
		
		//looping through distances
		for(int k=0;k<arg1.getValueStringArray("rgsDistanceNeighborhoods").length;k++){
			
			//initializing spatial weights matrix
			try {
				spw1 = 	new SpatialWeightsMatrix(bio1,
						Range.closed(Double.parseDouble(arg1.getValueStringArray("rgsDistanceNeighborhoods")[k].split("-")[0]),Double.parseDouble(arg1.getValueStringArray("rgsDistanceNeighborhoods")[k].split("-")[1])),
						Range.closed(Double.parseDouble(arg1.getValueStringArray("rgsDirectionNeighborhoods")[k].split("-")[0]),Double.parseDouble(arg1.getValueStringArray("rgsDirectionNeighborhoods")[k].split("-")[1])),
						Range.closed(Double.parseDouble(arg1.getValueStringArray("rgsTimeNeighborhoods")[k].split("-")[0]),Double.parseDouble(arg1.getValueStringArray("rgsTimeNeighborhoods")[k].split("-")[1])),
						arg1.getValueString("sWeight"), 
						sMode,
						false);
			}catch(Exception e){
				e.printStackTrace();
				return;
			}
			
			//loading permutations
			rgp1 = new Permutation[arg1.getValueInt("iMCMCChains")];
			for(int i=0;i<arg1.getValueInt("iMCMCChains");i++){
				rgp1[i]=new Permutation<Integer>(spw1.getVertexIDs());
				rgp1[i].loadRandomPermutation();	
			}
			
			//looping through values
			iCounter=0;
			for(String s:bio1.axsObservation.getIDs()){
			
				//updating progress
				iCounter++;
				System.out.println("Analyzing " + s + ", taxon " + iCounter + " of " + bio1.axsObservation.size() + ", iteration " + (k+1) + " of " + arg1.getValueStringArray("rgsDistanceNeighborhoods").length + "...");
			
				//updating cluster iterator
				cit1.next();
				if(cit1.bInclude==false){
					continue;
				}
				
				//loading values
				for(String t:bio1.axsSample.getIDs()){
					spw1.getVertex(bio1.axsSample.getIndex(t)).put("dValue", bio1.getValueByIDs(s, t));
				}
				
				//loading geospatial statistics object
				spa1 = new SpatialAutocorrelation(spw1);
				
				//running analysis
				if(arg1.getValueString("sAnalysis").startsWith("MoransI")){
					spa1.calculateMoransIMCMC(rgp1,arg1.getValueInt("iMCMCIterations"),true);
					if(Double.isInfinite(spa1.mrn1.dObsMoransI)){
						spa1.mrn1.dObsMoransI=Double.NaN;
						spa1.mrn1.dZMoransI=Double.NaN;
						spa1.mrn1.dPValueMoransI=Double.NaN;
						spa1.mrn1.iPValueIterations=0;
					}
					
					lstOut.add(s + "," + spa1.mrn1.dObsMoransI + 
							"," + spa1.mrn1.dPValueMoransI + 
							"," + spa1.mrn1.dZMoransI + 
							"," + spw1.rngDistances.lowerEndpoint() + "-" + spw1.rngDistances.upperEndpoint() +  
							"," + spw1.rngDirections.lowerEndpoint() + "-" + spw1.rngDirections.upperEndpoint() +
							"," + spw1.rngTimeDifferences.lowerEndpoint() + "-" + spw1.rngTimeDifferences.upperEndpoint() +
							"," + spa1.mrn1.iPValueIterations + 
							"," + spw1.iNonzeroWeights + 
							"," + spw1.dFractionWithNeighbors + 
							"," + bio1.getNonzeroCount(bio1.axsObservation, s) + 
							"," + bio1.getMean(bio1.axsObservation, s));
				}else if(arg1.getValueString("sAnalysis").equals("MoranScatter")){
					lstOut=spa1.calculateMoranScatterPlot();
					break;
				}	
			}
		}
		 
		//printing results and writing completion file
		if(arg1.getValueInt("iTaskID")==-9999 || arg1.getValueInt("iTotalTasks")==-9999){	
			if(arg1.getValueBoolean("bOutputData")){
				DataIO.writeToFile(outputDataTable(bio1), arg1.getValueString("sOutputPath").replace(".csv",".data.csv"));
			}
			if(arg1.getValueBoolean("bOutputDistances")){
				//TODO write code here to output all distances
				//gph1 = loadGraph(bio1, sMode, spw1, true);
				//DataIO.writeToFile(outputDistances(gph1), arg1.getValueString("sOutputPath").replace(".csv",".distances.csv"));
			}
			DataIO.writeToFile(lstOut, arg1.getValueString("sOutputPath"));
		}else{
			if(arg1.getValueBoolean("bOutputData") && arg1.getValueInt("iTaskID")==arg1.getValueInt("iTotalTasks")){	
				DataIO.writeToFile(outputDataTable(bio1), arg1.getValueString("sOutputPath").replace(".csv",".data.csv"));
			}
			if(arg1.getValueBoolean("bOutputDistances")){
				//TODO write code here to output all distances
				//gph1 = loadGraph(bio1, sMode, spw1, true);
				//DataIO.writeToFileWithCompletionFile(outputDistances(gph1), arg1.getValueString("sOutputPath").replace(".csv",".distances.csv"), arg1.getValueInt("iTaskID"));
			}
			DataIO.writeToFileWithCompletionFile(lstOut, arg1.getValueString("sOutputPath"), arg1.getValueInt("iTaskID"));
		}
		
		//terminating
		System.out.println("Done.");
	}
	
	/**
	 * Outputs data time differences and geographic distances
	 * @param gph1 Graph with distances
	 * @return List of distances and time differences
	 */
	private static ArrayList<String> outputDistances(SpatialWeightsMatrix spw1){
		
		//lstOut = output
		//lstTime = list of time differences
		//lstDist = list of differences
		
		ArrayList<String> lstOut;
		ArrayList<Double> lstTime;
		ArrayList<Double> lstDist;
		
		lstTime = spw1.getEdgeProperties("dTimeDifference");
		lstDist = spw1.getEdgeProperties("dLength");
		lstOut = new ArrayList<String>(lstTime.size()+1);
		lstOut.add("Geographic distance" +
				",Time difference");
		for(int i=0;i<lstDist.size();i++){
			lstOut.add(lstDist.get(i) + "," + lstTime.get(i));
		}
		return lstOut;
	}
	
	
	/**
	 * Outputs table of data that was analyzed
	 * @param bio1 BIOM object
	 */
	private static ArrayList<String> outputDataTable(BiomIO bio1){
		
		//lstOut = output
		//sbl1 = current output line
		
		ArrayList<String> lstOut;
		StringBuilder sbl1;
		
		lstOut = new ArrayList<String>(bio1.axsSample.size()+1);
		sbl1 = new StringBuilder();
		sbl1.append("SampleID,Latitude,Longitude,DateTime");
		for(int j=0;j<bio1.axsObservation.size();j++){
			sbl1.append("," + bio1.axsObservation.getID(j));
		}
		lstOut.add(sbl1.toString());
		for(int i=0;i<bio1.axsSample.size();i++){
			sbl1 = new StringBuilder();
			sbl1.append(bio1.axsSample.getID(i) + "," + bio1.axsSample.getMetadata(i).get("latitude") + "," + bio1.axsSample.getMetadata(i).get("longitude") + "," + bio1.axsSample.getMetadata(i).get("datetime"));
			for(int j=0;j<bio1.axsObservation.size();j++){
				sbl1.append("," + bio1.getValueByIndices(j, i));
			}
			lstOut.add(sbl1.toString());
		}
		return lstOut;
	}
}
