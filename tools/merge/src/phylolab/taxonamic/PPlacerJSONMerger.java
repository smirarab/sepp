package phylolab.taxonamic;

import com.google.gson.Gson;
import com.google.gson.GsonBuilder;
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
import phylolab.NewickTokenizer;

public class PPlacerJSONMerger {
	private Pattern seqNamePattern;
	LinkedList<SortedSet<String>> stack = new LinkedList<SortedSet<String>>();
	private Pattern edgeNumPattern;
	private Set<String> sequences = new TreeSet<String>();

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

	public HashMap<String, String> readTree(String tree) {
		HashMap<String, String> res = new HashMap<String, String>();
		NewickTokenizer tokenizer = new NewickTokenizer(tree, false);
		if (!"(".equals(tokenizer.nextToken())) {
			throw new RuntimeException("The tree does not start with a (");
		}
		do {
			String token = tokenizer.nextToken();
			if ("(".equals(token)) {
				this.stack.addLast(new TreeSet<String>());
			} else if (token.startsWith(")")) {
				if (this.stack.size() > 0) {
					String edgeNum = getEdgeNum(token);
					Collection<String> top = this.stack.getLast();
					addBipartition(res, top, edgeNum);
					this.stack.removeLast();
					if (this.stack.size() > 0){
						this.stack.getLast().addAll(top);
					}
				}
			} else {
				if (";".equals(token)) {
					continue;
				}
				String seqId = this.seqNamePattern.matcher(token).replaceAll(
						"$1");
                                //System.out.println(seqId);
				if (this.stack.size() > 0) {
					this.stack.getLast().add(seqId);
				}
				String edgeNum = getEdgeNum(token);
				this.sequences.add(seqId);
				addBipartition(res, Arrays.asList(new String[] { seqId }),
						edgeNum);
			}
		} while (tokenizer.hasNext());

		return res;
	}

	private HashMap<String, String> addInverseBipartitions(
			HashMap<String, String> bipartitions) {
		HashMap<String, String> newBipartitions = (HashMap<String,String>) bipartitions.clone();
		for (String bipartition : bipartitions.keySet()) {
			String[] seqs = bipartition.split("\\(\\)");
			if (seqs.length == 1) {
				continue;
			}
			TreeSet<String> seqSet = new TreeSet<String>();
			Collections.addAll(seqSet, seqs);
			TreeSet<String> newBipartition = new TreeSet<String>(this.sequences);
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

	private static void writeGSONFile(String fileName, JSONObject json) throws IOException{
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
	public HashMap<String, String> mapTrees(String relabel, String base) {
		HashMap<String, String> res = new HashMap<String, String>();
		HashMap<String, String> baseBipartitions = readTree(base);
		baseBipartitions = addInverseBipartitions(baseBipartitions);
		HashMap<String, String> relabelBipartitions = readTree(relabel);
		for (String bipartition : relabelBipartitions.keySet()) {
			String oldLabel = (String) relabelBipartitions.get(bipartition);
			String newLabel = (String) baseBipartitions.get(bipartition);
			res.put(oldLabel, newLabel);
		}
		return res;
	}

	private void relabelJson(String baseTree, JSONObject json,
			HashMap<String, Double> mainEdgeLen) {
		String jsonTree = json.getString("tree");
		HashMap<String, String> labelMap = mapTrees(jsonTree, baseTree);

		JSONArray placements = json.getJSONArray("placements");
		JSONArray fields = json.getJSONArray("fields");
                int locEdgeNum=0;
                while (! fields.getString(locEdgeNum).equals("edge_num")) {
                        locEdgeNum++;
                }
		
		for (Iterator<JSONObject> iterator = placements.iterator(); iterator.hasNext();) {
			JSONObject placement = iterator.next();
			JSONArray p = placement.getJSONArray("p");
			for (Iterator<JSONArray> itp = p.iterator(); itp.hasNext();) {
				JSONArray pr = itp.next();
				String newLab = (String) labelMap.get(pr.getString(locEdgeNum));
				pr.set(locEdgeNum, new Integer(newLab));

				if (pr.getDouble(3) > ((Double) mainEdgeLen.get(newLab))
						.doubleValue())
					pr.set(3, Double.valueOf(((Double) mainEdgeLen.get(newLab))
							.doubleValue() * 0.99D));
			}
		}
	}

	public static void main(String[] args) {
		boolean sorted = false;

		if (args.length < 3) { 
			System.out
					.println("Usage: merge.jar <json files directory> <base tree file> <output> [-s]\n" +
							"\t\t<json files directory>: the directory with all pplacer results (.json files)\n" +
							"\t\t<base tree file>: The base tree file\n" +
							"\t\t<output>: output json file name\n" +
							"\t\t-s: (optional) sort the fragments by name.\n\n" + 
							"\t\t ( NOTE: Instead of providing json directory and base tree files, you can use a '-'.\n"+
							"\t\t In this case base trees and json files are read from standard input.\n"+
							"\t\t First line of standard input should give the global base tree. Subsequent lines\n" +
							"\t\t should give a labeled tree followed by the location of .json file for each subset.\n" +
							"\t\t After these pairs of lines are given for all subsets, an empty line should indicate end of input. )");

			System.exit(1);
		}
		
		if ( (args.length ==4) && (args[3].equals("-s")) ) {
			sorted = true;
		}
		String mainTree = "";
		HashMap<String, Double> mainEdgeLen = new HashMap<String, Double>();
		String name;
		try {
			
			String baseFn = args[1];
			BufferedReader inm = new BufferedReader(new InputStreamReader(
					"-".equals(baseFn)? System.in : new FileInputStream(baseFn)));
			mainTree = inm.readLine();
			if (!"-".equals(baseFn)) {inm.close();}
			
			List<String> trees = new ArrayList<String>();
			List<String> jsonLocations = new ArrayList<String>();			
			String jsonDir = args[0];						
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
						String baseTreeFn = jsonFile.getAbsolutePath().replace("json",
								"labeled.tree");
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
			System.err.println("json locations: " + jsonLocations);
			
			mainTree = mainTree.replaceAll("'","");
			//System.out.println(mainTree);
			Matcher matcher = Pattern.compile(":([^\\[]*)\\[([^\\]]*)\\]")
					.matcher(mainTree);

			while (matcher.find()) {
				Double len = new Double(matcher.group(1));
				name = matcher.group(2);
				mainEdgeLen.put(name, len);
			}
			
			JSONObject resultsJson = new JSONObject();
			resultsJson.put("tree", mainTree);
			JSONArray resultsPlacements = new JSONArray();
			JSONArray fields = null;
		
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
					
					PPlacerJSONMerger merger = new PPlacerJSONMerger();
					merger.relabelJson(baseTree, json, mainEdgeLen);
	
					json.put("tree", mainTree);								
	
					// UPDATE the global merge file
					resultsPlacements.addAll(json.getJSONArray("placements"));
					fields = json.getJSONArray("fields");
					
					// Unnecessary IO
					if (false) {
						writeGSONFile(jsonFile.replace(".json", ".merged.json"),json);
					}
	
				} catch (JsonSyntaxException e) {
					System.err.println(e.getLocalizedMessage());
					System.err.println("The above warnning ignored. continue ...");
				}
			}	
			if (sorted) {
			    TreeSet<JSONObject> sortedPlacements = new TreeSet<JSONObject>(new Comparator<JSONObject>() {
					@Override
					public int compare(JSONObject o1, JSONObject o2) {					
						String name1 = o1.optString("n") + o1.optString("nm");
						String name2 = o2.optString("n") + o2.optString("nm");
						return name1.compareTo(name2);
					}
				});
                            for (int i=0 ; i<resultsPlacements.size() ; i++) {
                              sortedPlacements.add(resultsPlacements.getJSONObject(i));
                            }
			    resultsPlacements = new JSONArray();
			    resultsPlacements.addAll(sortedPlacements);
			}
			resultsJson.put("placements", resultsPlacements);
			resultsJson.put("metadata", JSONObject.fromObject("{\"invocation\":" +
					"\"SEPP-generated json file (sepp 2).\"}"));
			resultsJson.put("version", 1);
			resultsJson.put("fields", fields);
			writeGSONFile(args[2], resultsJson);
		} catch (FileNotFoundException e) {
			System.err.println("File not found: \n" + e.getMessage());
			System.exit(1);
		} catch (IOException e) {
			System.err.println("I/O Error: \n" + e.getMessage());
			System.exit(1);
		}		
	}
}
