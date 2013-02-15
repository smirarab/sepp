package phylolab.taxonamic;

import com.google.gson.Gson;
import com.google.gson.GsonBuilder;
import com.google.gson.JsonArray;
import com.google.gson.JsonElement;
import com.google.gson.JsonParser;
import com.google.gson.JsonSyntaxException;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.FilenameFilter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;
import java.util.Set;
import java.util.SortedSet;
import java.util.TreeSet;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import net.sf.json.JSONArray;
import net.sf.json.JSONObject;
import net.sf.json.JSONSerializer;
import phylolab.NewickTokenizer;

public class PPlacerJSONMerger {
    private Pattern seqNamePattern;
    private Pattern edgeNumPattern;

    private HashMap<String,JSONArray> nameToAllPlacements = new HashMap<String, JSONArray>();
    private HashMap<String,Double> nameToCummulativeLWR = new HashMap<String, Double>();

    static String join(Collection<String> s, String delimiter) {
	StringBuilder builder = new StringBuilder();
	Iterator<String> iter = s.iterator();
	while (iter.hasNext()) {
	    builder.append(iter.next());
	    if (!iter.hasNext()) {
		break;
	    }
	    builder.append(delimiter); 
	}
	return builder.toString();
    }

    public PPlacerJSONMerger() {
	this.edgeNumPattern = Pattern.compile(".*[\\[\\{]([0-9]*)[\\]\\}]");
	this.seqNamePattern = Pattern.compile("'*([^:']*)'*(:.*)*");
    }

    public HashMap<String, String> readTree(String tree, Set<String> sequences) {
	HashMap<String, String> res = new HashMap<String, String>();
	LinkedList<SortedSet<String>> stack = new LinkedList<SortedSet<String>>();

	NewickTokenizer tokenizer = new NewickTokenizer(tree, false);
	if (!"(".equals(tokenizer.nextToken())) {
	    throw new RuntimeException("The tree does not start with a (");
	}
	do {
	    String token = tokenizer.nextToken();
	    if ("(".equals(token)) {
		stack.addLast(new TreeSet<String>());
	    } else if (token.startsWith(")")) {
		if (stack.size() > 0) {
		    String edgeNum = getEdgeNum(token);
		    Collection<String> top = stack.getLast();
		    addBipartition(res, top, edgeNum);
		    stack.removeLast();
		    if (stack.size() > 0){ 
			stack.getLast().addAll(top);
		    }
		}
	    } else {
		if (";".equals(token)) {
		    continue;
		}
		String seqId = this.seqNamePattern.matcher(token).replaceAll("$1");
		//System.out.println(seqId);
		if (stack.size() > 0) {
		    stack.getLast().add(seqId);
		}
		String edgeNum = getEdgeNum(token);
		sequences.add(seqId);
		addBipartition(res, Arrays.asList(new String[] { seqId }),
			edgeNum);
	    }
	} while (tokenizer.hasNext());

	return res;
    }

    private HashMap<String, String> addInverseBipartitions(
	    HashMap<String, String> bipartitions, Set<String> sequences) {
	HashMap<String, String> newBipartitions = (HashMap<String,String>) bipartitions.clone();
	for (String bipartition : bipartitions.keySet()) {
	    String[] seqs = bipartition.split("\\(\\)");
	    if (seqs.length == 1) {
		continue;
	    }
	    TreeSet<String> seqSet = new TreeSet<String>();
	    Collections.addAll(seqSet, seqs);
	    TreeSet<String> newBipartition = new TreeSet<String>();
	    newBipartition.removeAll(seqSet);
	    addBipartition(newBipartitions, newBipartition,
		    (String) bipartitions.get(bipartition));
	}
	return newBipartitions;
    }

    private String getEdgeNum(String token) {
	Matcher match = this.edgeNumPattern.matcher(token);
	if (match.find()) {
	    return match.group(1);
	}
	return "?";
    }

    private void addBipartition(HashMap<String, String> res,
	    Collection<String> top, String edgeNum) {
	String name = nameBipartition(top);
	res.put(name, edgeNum);
    }

    private String nameBipartition(Collection<String> top) {
	return join(top, "()");
    }

    /**
     * Called from relabelJson to map between main tree and an individual tree
     * @param from
     * @param to
     * @return
     */
    private HashMap<String, String> mapTreeBranchNames(String from, String to) {
	HashMap<String, String> res = new HashMap<String, String>();
	Set<String> sequences =  new TreeSet<String>();
	HashMap<String, String> baseBipartitions = readTree(to,sequences);
	baseBipartitions = addInverseBipartitions(baseBipartitions, sequences);
	HashMap<String, String> fromBipartitions = readTree(from, sequences);
	for (String bipartition : fromBipartitions.keySet()) {
	    String fromLabel = (String) fromBipartitions.get(bipartition);
	    String toLabel = (String) baseBipartitions.get(bipartition);
	    res.put(fromLabel, toLabel);
	}
	return res;
    }

    private void relabelAndUpdatePlacements(String baseTree, JSONObject json,
	    HashMap<String, Double> mainEdgeLen, int rmUnerscore) {

	String jsonTree = json.getString("tree");
	HashMap<String, String> labelMap = mapTreeBranchNames(jsonTree, baseTree);

	JSONArray placements = json.getJSONArray("placements");

	for (Iterator<JSONObject> iterator = placements.iterator(); iterator.hasNext();) {	    
	    JSONObject placement = iterator.next();
	    /*
	     * replace nm with n
	     */
	    if (placement.containsKey("nm")) {
		JSONArray n = new JSONArray();
		for (Iterator nmIt = placement.getJSONArray("nm").iterator(); nmIt.hasNext();) {
		    n.add(((JSONArray) nmIt.next()).getString(0));
		}
		placement.put("n", n); 
		placement.discard("nm");
	    }
	    JSONArray n = placement.getJSONArray("n");
	    for (int i = 0; i < n.size(); i++) {
		String name = n.getString(i);
		String newName = name;
		Double prior = 1.;
		if (rmUnerscore != 0) {
		    String [] nameSplit = name.split("_");
		    if (nameSplit.length < rmUnerscore + 1) {
			throw new RuntimeException("Fragments names should have at least " + (rmUnerscore) + " underscores.");
		    }
		    /*
		     * Find and update new Name
		     */
		    StringBuilder newNameBuilder = new StringBuilder(nameSplit[0]);
		    for (int j = 1; j < nameSplit.length - 4; j++) {
			newNameBuilder.append("_");
			newNameBuilder.append(nameSplit[j]);
		    }		
		    newName = newNameBuilder.toString();

		    n.set(i, newName);
		    // Find the prior probability
		    prior = Integer.parseInt(nameSplit[nameSplit.length - 1])/1000000.0;
		} 

		JSONArray p = placement.getJSONArray("p");
		Double sum = 0.;		
		JSONArray current = nameToAllPlacements.containsKey(newName) ?  nameToAllPlacements.get(newName) : new JSONArray();
		for (Iterator<JSONArray> itp = p.iterator(); itp.hasNext();) {
		    JSONArray precord = itp.next();
		    JSONArray newRecord = new JSONArray();
		    newRecord.addAll(precord);
		    /*
		     * Adjust the placement edge label
		     */
		    String newLab = (String) labelMap.get(precord.getString(0));
		    newRecord.set(0, new Integer(newLab));
		    /*
		     * Adjust edge length values to correspond to somewhere on the main tree.
		     * TODO: This is pretty bad. We should fix this. 
		     */
		    if (precord.getDouble(3) > ((Double) mainEdgeLen.get(newLab)).doubleValue())
			newRecord.set(3, Double.valueOf(((Double) mainEdgeLen.get(newLab))
				.doubleValue() * 0.99D));
		    /*
		     * Adjust the weighted likelihood by alignment probability as a prior.
		     */
		    Double newLWR = new Double(prior * precord.getDouble(2));
		    newRecord.set(2, newLWR);
		    sum += newLWR;

		    current.add(newRecord);
		}	    
		nameToAllPlacements.put(newName, current);		   	    

		nameToCummulativeLWR.put(newName, nameToCummulativeLWR.containsKey(newName) ?
			nameToCummulativeLWR.get(newName) + sum : sum);	   
	    }
	}
    }

    private void mergePlacementsOFSameFragments(JSONArray all) {

	for (String fragment : this.nameToAllPlacements.keySet()) {
	    JSONArray placements = nameToAllPlacements.get(fragment);
	    Double sum = nameToCummulativeLWR.get(fragment);
	    System.err.println(fragment + " " + placements.size());
	    /*
	     * Normalize weighted likelihood ratios
	     */
	    for (Iterator<JSONArray> itp = placements.iterator(); itp.hasNext();) {
		JSONArray pr = itp.next();		    
		pr.set(2, new Double(pr.getDouble(2) / sum));
	    }		
	    JSONObject placement = new JSONObject();
	    placement.put("p", placements);
	    JSONArray n = new JSONArray();
	    n.add(fragment);
	    placement.put("n", n);
	    all.add(placement);

	}
    }

    private void writeGSONFile(String fileName, JSONObject json) throws IOException{
	FileWriter writer = new FileWriter(fileName);
	String string = json.toString();
	Gson gson = new GsonBuilder().setPrettyPrinting().create();
	JsonElement parsed = new JsonParser().parse(string);
	string = gson.toJson(parsed).replaceAll("([0-9],)\\n\\s*",
		"$1 ").replaceAll("(\\[)\\n\\s*", "$1").replaceAll(
			"\\n(\\s*)\\]", "\\]");
	writer.write(string);
	writer.close();
    }


    public static void errout(){
	System.out
	.println("Usage: merge.jar <json files directory> <base tree file> <output> [-r N] [-s]\n" +
		"\t\t<json files directory>: the directory with all pplacer results (.json files)\n" +
		"\t\t<base tree file>: The base tree file\n" +
		"\t\t<output>: output json file name\n" +
		"\t\t-s: (optional) sort the fragments by name.\n" +
		"\t\t-r N: (optional) rename fragments to remove everything after Nth _ from the end.\n" +
		"\t\t       Merge multiple placement of the same fragment into one entry, considering prior probabilities.\n" +
		"\t\t       Treat everything after the last _ as a prior probability out of 1000,000.\n\n" + 
		"\t\t ( NOTE: Instead of providing json directory and base tree files, you can use a '-'.\n"+
		"\t\t In this case base trees and json files are read from standard input.\n"+
		"\t\t First line of standard input should give the global base tree. Subsequent lines\n" +
		"\t\t should give a labeled tree followed by the location of .json file for each subset.\n" +
		"\t\t After these pairs of lines are given for all subsets, an empty line should indicate end of input. )");

	System.exit(1);
    }

    public JSONObject mergeJsonFiles(String mainTree, List<String> trees, 
	    List<String> jsonLocations, boolean sorted, int rmUnderscore) throws IOException {
	HashMap<String, Double> mainEdgeLen = new HashMap<String, Double>();
	/*
	 * Find the length of individual edges in the main tree
	 */
	mainTree = mainTree.replaceAll("'","");
	Matcher edgeLenMatcher = Pattern.compile(":([^\\[]*)\\[([^\\]]*)\\]")
		.matcher(mainTree);
	while (edgeLenMatcher.find()) {
	    mainEdgeLen.put(edgeLenMatcher.group(2), new Double(edgeLenMatcher.group(1)));
	}

	JSONObject resultsJson = new JSONObject();
	/*
	* Make the main tree, the "tree" of output json
	*/
	resultsJson.put("tree", mainTree);
	JSONArray resultsPlacements = new JSONArray();
	JSONArray fields = null;

	/*
	 * Read json files one by one, relabel them, and merge their placements into the output
	 */
	for (int i = 0; i < trees.size(); i++) {
	    try {					
		String jsonFile = jsonLocations.get(i);
		BufferedReader in = new BufferedReader(new FileReader(jsonFile));
		StringBuffer jsonString = new StringBuffer();
		String str;
		while ((str = in.readLine()) != null) {
		    jsonString.append(str);
		}
		JSONObject json = JSONObject.fromObject(jsonString.toString());

		String baseTree = trees.get(i);

		// Update the placement section of the json			
		this.relabelAndUpdatePlacements(baseTree, json, mainEdgeLen, rmUnderscore);							
		// UPDATE the global merge file
		// resultsPlacements.addAll(json.getJSONArray("placements"));

		fields = json.getJSONArray("fields");

		// Unnecessary IO
		//writeGSONFile(jsonFile.replace(".json", ".merged.json"),json);

	    } catch (JsonSyntaxException e) {
		System.err.println(e.getLocalizedMessage());
		System.err.println("The above warnning ignored. continue ...");
	    }
	}	
	
	/*
	 * Merge multiple placements for the same fragment
	 */
	this.mergePlacementsOFSameFragments(resultsPlacements);

	if (sorted) {
	    TreeSet<JSONObject> sortedPlacements = new TreeSet<JSONObject>(new Comparator<JSONObject>() {
		@Override
		public int compare(JSONObject o1, JSONObject o2) {					
		    String name1 = o1.optString("n") + o1.optString("nm");
		    String name2 = o2.optString("n") + o2.optString("nm");
		    return name1.compareTo(name2);
		}
	    });
	    sortedPlacements.addAll(resultsPlacements);
	    resultsPlacements = new JSONArray();
	    resultsPlacements.addAll(sortedPlacements);
	}
	resultsJson.put("placements", resultsPlacements);
	resultsJson.put("metadata", JSONObject.fromObject("{\"invocation\":" +
		"\"SEPP-generated json file (sepp 2).\"}"));
	resultsJson.put("version", 1);
	resultsJson.put("fields", fields);
	return resultsJson;
    }
    
    public static void main(String[] args) {			    

	if (args.length < 3) { 
	    errout();
	}

	String jsonDir = args[0];
	String baseFn = args[1];
	String outfilename = args[2];
	boolean sorted = false;
	int rmUnderscore = 0;
	String mainTree = "";
	List<String> trees = new ArrayList<String>();
	List<String> jsonLocations = new ArrayList<String>();			


	/*
	 * Parse optional arguments
	 */
	for (int i = 3; i < args.length; i++) {
	    if (args[i].equals("-s")) {
		sorted = true;
	    } else if (args[i].equals("-r")) {
		if (i+1 >= args.length) {
		    System.out.println("-r needs to be followd by a number.");
		    System.exit(1);
		}
		i++;
		rmUnderscore = Integer.parseInt(args[i]);
	    }
	}

	try {
	    /*
	     * Read the main tree from the file
	     */
	    BufferedReader inm = new BufferedReader(new InputStreamReader(
		    "-".equals(baseFn)? System.in : new FileInputStream(baseFn)));
	    mainTree = inm.readLine();
	    if (!"-".equals(baseFn)) {inm.close();}

	    /*
	     * Read in json files and their associated labeled tree
	     */
	    if (!"-".equals(jsonDir)) {
		File[] files = new File(jsonDir).listFiles(new FilenameFilter() {
		    public boolean accept(File dir, String name) {
			return (name.indexOf(".json") >= 0)
				&& (name.indexOf("merged") < 0);
		    }
		});

		for (int i = 0; i < files.length; i++) {
		    File jsonFile = files[i];
		    try {
			String baseTreeFn = jsonFile.getAbsolutePath().replace("json", "labeled.tree");
			BufferedReader in = new BufferedReader(new InputStreamReader(
				new FileInputStream(baseTreeFn)));
			trees.add(in.readLine()); 
			in.close();

			jsonLocations.add(jsonFile.getAbsolutePath());
		    } catch (FileNotFoundException e) {
			System.err.println(e.getLocalizedMessage());
			System.err.println("The above warnning ignored. continue ...");
		    } catch (IOException e) {
			System.err.println(e.getLocalizedMessage());
			System.err.println("The above warnning ignored. continue ...");
		    }
		}
	    } else {
		String line = "";
		while (true){					
		    line = inm.readLine();				
		    if (line.length() == 0) {
			break;
		    }
		    trees.add(line);
		    jsonLocations.add(inm.readLine());
		}
	    }
	    //System.err.println("json locations: " + jsonLocations);


	    PPlacerJSONMerger merger = new PPlacerJSONMerger();
	    JSONObject merged = merger.mergeJsonFiles(mainTree, trees, jsonLocations, sorted, rmUnderscore);
	    merger.writeGSONFile(outfilename, merged);

	} catch (IOException e) {
	    System.err.println("I/O Error: \n" + e.getMessage());
	    System.exit(1);
	}		
    }
}
