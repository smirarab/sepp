package phylolab.taxonamic;

import com.google.gson.Gson;
import com.google.gson.GsonBuilder;
import com.google.gson.JsonElement;
import com.google.gson.JsonParser;
import com.google.gson.JsonSyntaxException;

import edu.rice.cs.bioinfo.programs.phylonet.algos.lca.SchieberVishkinLCA;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.io.ParseException;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.TNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.sti.STINode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.sti.STITree;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.FilenameFilter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Hashtable;
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

public class JSONMerger {
  private Pattern seqNamePattern = Pattern.compile("'*([^:']*)'*(:.*)*");
  private Pattern edgeNumPattern = Pattern.compile(".*[\\[\\{]([0-9]*)[\\]\\}]");

  private HashMap < String, JSONArray > nameToAllPlacements = new HashMap < String, JSONArray > ();
  private HashMap < String, Double > nameToCummulativeLWR = new HashMap < String, Double > ();
  private String mainTree;
  private List < String > trees;
  private List < String > jsonLocations;
  private boolean sorted;
  private int rmUnderscore;
  private HashMap < String, Double > mainEdgeLen;
  private STITree < TaxonomyData > taxonomy;
  private Hashtable < String, String > jsonNameToTaxonId;
  private Double threshold;
  private FileWriter cw;
  private boolean distribution;
  private boolean pushDown;
  private Double cutoff = 0D;

  static String join(Collection < String > s, String delimiter) {
    StringBuilder builder = new StringBuilder();
    Iterator < String > iter = s.iterator();
    while (iter.hasNext()) {
      builder.append(iter.next());
      if (!iter.hasNext()) {
        break;
      }
      builder.append(delimiter);
    }
    return builder.toString();
  }

  public JSONMerger(String mainTree, List < String > trees,
    List < String > jsonLocations, boolean sorted, int rmUnderscore,
    STITree < TaxonomyData > taxonomy, Hashtable < String, String > jsonNameToTaxonId, Double threshold, FileWriter classificationWriter, boolean pushUp, boolean distribution, Double cutoff) {
    mainEdgeLen = new HashMap < String, Double > ();
    this.mainTree = mainTree;
    this.trees = trees;
    this.jsonLocations = jsonLocations;
    this.sorted = sorted;
    this.rmUnderscore = rmUnderscore;
    this.taxonomy = taxonomy;
    this.jsonNameToTaxonId = jsonNameToTaxonId;
    this.threshold = threshold;
    this.cw = classificationWriter;
    this.pushDown = !pushUp;
    this.distribution = distribution;
    this.cutoff = cutoff;
  }

  /**
   * Reads the tree, and adds the sequence names to sequences
   * @param tree
   * @param sequences
   * @return
   */
  public HashMap < String, String > readTree(String tree, Set < String > sequences) {
    HashMap < String, String > res = new HashMap < String, String > ();
    LinkedList < SortedSet < String >> stack = new LinkedList < SortedSet < String >> ();
    //System.out.println("Reading\n" + tree);
    NewickTokenizer tokenizer = new NewickTokenizer(tree, false);
    if (!"(".equals(tokenizer.nextToken())) {
      throw new RuntimeException("The tree does not start with a (");
    }
    do {
      String token = tokenizer.nextToken();
      if ("(".equals(token)) {
        stack.addLast(new TreeSet < String > ());
      } else if (token.startsWith(")")) {
        if (stack.size() > 0) {
          String edgeNum = getEdgeNum(token);
          Collection < String > top = stack.getLast();
          addBipartition(res, top, edgeNum);
          stack.removeLast();
          if (stack.size() > 0) {
            if (top == null) {
              System.out.println();
            }
            stack.getLast().addAll(top);
          }
        }
      } else {
        if (";".equals(token)) {
          continue;
        }
        String seqId = this.seqNamePattern.matcher(token).replaceAll("$1");
        if (stack.size() > 0) {
          stack.getLast().add(seqId);
        }
        String edgeNum = getEdgeNum(token);
        sequences.add(seqId);
        addBipartition(res, Arrays.asList(new String[] {
            seqId
          }),
          edgeNum);
      }
    } while (tokenizer.hasNext());

    return res;
  }

  private HashMap < String, String > addInverseBipartitions(
    HashMap < String, String > bipartitions, Set < String > sequences) {
    HashMap < String, String > newBipartitions = (HashMap < String, String > ) bipartitions.clone();
    for (String bipartition: bipartitions.keySet()) {
      String[] seqs = bipartition.split("\\(\\)");
      if (seqs.length == 1) {
        continue;
      }
      TreeSet < String > seqSet = new TreeSet < String > ();
      Collections.addAll(seqSet, seqs);
      TreeSet < String > newBipartition = new TreeSet < String > ();
      newBipartition.removeAll(seqSet);
      addBipartition(newBipartitions, newBipartition, (String) bipartitions.get(bipartition));
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

  private void addBipartition(HashMap < String, String > res,
    Collection < String > top, String edgeNum) {
    String name = nameBipartition(top);
    res.put(name, edgeNum);
  }

  private String nameBipartition(Collection < String > top) {
    return join(top, "()");
  }

  /**
   * Called from relabelJson to map between main tree and an individual tree
   * @param from
   * @param to
   * @return
   */
  private HashMap < String, String > mapTreeBranchNames(String from, String to) {
    HashMap < String, String > res = new HashMap < String, String > ();
    Set < String > sequences = new TreeSet < String > ();
    HashMap < String, String > baseBipartitions = readTree(to, sequences);
    baseBipartitions = addInverseBipartitions(baseBipartitions, sequences);
    HashMap < String, String > fromBipartitions = readTree(from, sequences);
    for (String bipartition: fromBipartitions.keySet()) {
      String fromLabel = (String) fromBipartitions.get(bipartition);
      String toLabel = (String) baseBipartitions.get(bipartition);
      res.put(fromLabel, toLabel);
    }
    return res;
  }

  private void processJson(String originalTree, JSONObject json) {
    String jsonTree = json.getString("tree");
    HashMap < String, String > labelMap = mapTreeBranchNames(jsonTree, originalTree);

    JSONArray placements = json.getJSONArray("placements");
    JSONArray fields = json.getJSONArray("fields");
    
    int locEdgeNum=0;
    while (! fields.getString(locEdgeNum).equals("edge_num")) {
      locEdgeNum++;
    }

    for (Iterator < JSONObject > iterator = placements.iterator(); iterator.hasNext();) {
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
        if (rmUnderscore != 0) {
          String[] nameSplit = name.split("_");
          if (nameSplit.length < rmUnderscore + 1) {
            throw new RuntimeException("Fragments names should have at least " + (rmUnderscore) + " underscores.");
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
          prior = Integer.parseInt(nameSplit[nameSplit.length - 1]) / 1000000.0;
        }

        JSONArray p = placement.getJSONArray("p");
        Double sum = 0.;
        JSONArray current = nameToAllPlacements.containsKey(newName) ? nameToAllPlacements.get(newName) : new JSONArray();
        for (Iterator < JSONArray > itp = p.iterator(); itp.hasNext();) {
          JSONArray precord = itp.next();
          JSONArray newRecord = new JSONArray();
          newRecord.addAll(precord);
          /*
           * Adjust the placement edge label
           */
          String newLab = (String) labelMap.get(precord.getString(locEdgeNum));
          newRecord.set(locEdgeNum, new Integer(newLab));
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

  /**
   * Maps from json tree edges to the taxonomic tree
   * @return
   * @throws IOException
   * @throws ParseException
   */

  private HashMap < String, STINode < TaxonomyData >> mapJsonTreeToTaxonomy() throws IOException, ParseException {

    SchieberVishkinLCA lookup = new SchieberVishkinLCA(taxonomy);
    HashMap < String, STINode < TaxonomyData >> jsonToTaxonomy = new HashMap < String, STINode < TaxonomyData >> ();
    // stack 
    LinkedList < Set < TNode >> taxonomyNodeStack = new LinkedList < Set < TNode >> ();
    LinkedList < String > unmappedEdgeStack = new LinkedList < String > ();
    NewickTokenizer tokenizer = new NewickTokenizer(this.mainTree, false);

    if (!"(".equals(tokenizer.nextToken())) {
      throw new RuntimeException("The main tree does not start with a (");
    }
    // Read the json tree, and use LCA lookup to map between the two trees.
    do {
      String token = tokenizer.nextToken();
      if ("(".equals(token)) {
        taxonomyNodeStack.addLast(new HashSet < TNode > ());
      } else if (token.startsWith(")")) {
        if (taxonomyNodeStack.size() > 0) {
          String edgeNum = getEdgeNum(token);
          Set < TNode > top = taxonomyNodeStack.getLast();
          try {
            lookup.getLCA(top);
          } catch (Exception e) {
            System.out.println(e.toString());
          }
          TNode taxonomyNode = lookup.getLCA(top);
          if (pushDown) {
            jsonToTaxonomy.put(edgeNum, (STINode < TaxonomyData > ) taxonomyNode);
          } else {
            for (String edge: unmappedEdgeStack) {
              jsonToTaxonomy.put(edge, (STINode < TaxonomyData > ) taxonomyNode);
            }
            unmappedEdgeStack.remove();
            unmappedEdgeStack.addLast(edgeNum);
          }
          taxonomyNodeStack.removeLast();
          if (taxonomyNodeStack.size() > 0) {
            if (taxonomyNode == null) {
              System.out.println();
            }
            taxonomyNodeStack.getLast().add(taxonomyNode);
          }
        }
      } else {
        if (";".equals(token)) {
          continue;
        }
        String seqId = this.seqNamePattern.matcher(token).replaceAll("$1");
        String taxonID = this.jsonNameToTaxonId.get(seqId);
        TNode taxonomyNode = this.taxonomy.getNode(taxonID);

        if (taxonomyNodeStack.size() > 0) {
          if (taxonomyNode == null) {
            System.out.println();
          }
          taxonomyNodeStack.getLast().add(taxonomyNode);
        }
        String edgeNum = getEdgeNum(token);

        if (pushDown) {
          jsonToTaxonomy.put(edgeNum, (STINode < TaxonomyData > ) taxonomyNode);
        } else {
          unmappedEdgeStack.push(edgeNum);
        }
      }
    } while (tokenizer.hasNext());

    System.err.println(jsonToTaxonomy);

    return jsonToTaxonomy;
  }


  private void createMergedPlacements(JSONArray all) throws IOException, ParseException {

    HashMap < String, STINode < TaxonomyData >> jsonTreeIDToTaxonomyNode = null;
    if (taxonomy != null) {
      jsonTreeIDToTaxonomyNode = mapJsonTreeToTaxonomy();
    }

    for (String fragment: this.nameToAllPlacements.keySet()) {
      JSONArray placements = nameToAllPlacements.get(fragment);
      Double sum = nameToCummulativeLWR.get(fragment);
      //If returning distribution, only take as many placements necessary to reach some threshold
      // renormalize remaining LWR so that placements sum up to 1
      if (distribution) {
        System.out.println("Doing probability stuff\n");
        TreeSet < JSONArray > sortedPlacements = new TreeSet < JSONArray > (new Comparator < JSONArray > () {@
          Override
          public int compare(JSONArray o1, JSONArray o2) {
            double o = o1.getDouble(2) - o2.getDouble(2);
            if (o < 0) {
              o = 1;
            } else if (o > 0) {
              o = -1;
            } else {
              o = -1;
            }
            return (int) o;
          }
        });
        for (Iterator<JSONArray> it=placements.iterator() ; it.hasNext() ;) {
          sortedPlacements.add(it.next());
        }
        double total = 0;
        ArrayList < JSONArray > list = new ArrayList < JSONArray > ();
        for (Iterator < JSONArray > itp = sortedPlacements.iterator(); threshold > total && itp.hasNext();) {          
          JSONArray placementRecord = itp.next();
          if (placementRecord.getDouble(2) <= cutoff) {
            break;
          }
          total += placementRecord.getDouble(2);
          list.add(placementRecord);
        }
        placements = new JSONArray();
        placements.addAll(list);
        sum = sum * total;
      }

      Set < STINode < TaxonomyData >> lineages = new HashSet < STINode < TaxonomyData >> ();

      //System.err.println(fragment + " " + placements.size());
      /*
       * 1- Normalize weighted likelihood ratios
       * 2- update the probabilities for this fragment on the taxonomic tree.
       */
      for (Iterator < JSONArray > itp = placements.iterator(); itp.hasNext();) {
        JSONArray placementRecord = itp.next();
        Double probability = placementRecord.getDouble(2) / sum;
        if (jsonTreeIDToTaxonomyNode != null) {
          String edgeNumber = placementRecord.getString(0);
          //System.out.println(fragment + " " + edgeNumber + " " + probability);
          STINode < TaxonomyData > node = jsonTreeIDToTaxonomyNode.get(edgeNumber);
          while (node != null) {
            node.getData().probability += probability;
            lineages.add(node);
            node = node.getParent();
          }
        }
        placementRecord.set(2, probability);

      }
      JSONObject placement = new JSONObject();
      placement.put("p", placements);
      JSONArray n = new JSONArray();
      n.add(fragment);
      placement.put("n", n);
      all.add(placement);

      /*
       * Write out classification results for a current fragment
       */
      if (jsonTreeIDToTaxonomyNode != null) {
        for (STINode < TaxonomyData > lineage: lineages) {
          TaxonomyData data = lineage.getData();
          if (data.probability >= threshold || distribution) {

            cw.write(join(Arrays.asList(new String[] {
              fragment,
              data.id,
                data.name,
                data.rank,
                new DecimalFormat("##.0000").format(data.probability)
            }), ",") + "\n");
          }
          data.probability = 0D;
        }
      }
    }
  }

  private void writeGSONFile(String fileName, JSONObject json) throws IOException {
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


  public static void errout() {
    System.err
      .println("Usage: merge.jar <json files directory> <base tree file> <output> [-r N] [-s]\n" +
        "\t\t<json files directory>: the directory with all pplacer results (.json files)\n" +
        "\t\t<base tree file>: The base tree file\n" +
        "\t\t<output>: output json file name\n" +
        "\t\t-s: (optional) sort the fragments by name.\n" +
        "\t\t-u: (optional) push the fragments up the placement edge instead of pushing them down.\n" +
        "\t\t-t <taxonomy file>: (optional) The name of a taxonomy file. If provided, classification is also performed,\n" +
        "\t\t                     and results are written out to standard output.\n" +
        "\t\t-p N: (optional) A value between 0 and 1. When given with -t option, this specifies\n" +
        "\t\t      the minimum probability threshold for lineages appearing in ouput.\n" +
        "\t\t-m <name mapping file>: (optional) when given with -t option, this provides a comma-seperated\n" +
        "\t\t                        mapping between names in given taxonomy and json files.\n" +
        "\t\t-c <classification output file>: (optional) when given with -t option, results of classification are written\n" +
        "\t\t                        to this file.\n" +
        "\t\t-r N: (optional) rename fragments to remove everything after Nth _ from the end.\n" +
        "\t\t       Merge multiple placement of the same fragment into one entry, considering prior probabilities.\n" +
        "\t\t       Treat everything after the last _ as a prior probability out of 1000,000.\n\n" +
        "\t\t ( NOTE: Instead of providing json directory and base tree files, you can use a '-'.\n" +
        "\t\t In this case base trees and json files are read from standard input.\n" +
        "\t\t First line of standard input should give the global base tree. Subsequent lines\n" +
        "\t\t should give a labeled tree followed by the location of .json file for each subset.\n" +
        "\t\t After these pairs of lines are given for all subsets, an empty line should indicate end of input. )");

    System.exit(1);
  }

  public JSONObject mergeJsonFiles() throws IOException, ParseException {
    /*
     * Find the length of individual edges in the main tree
     */
    mainTree = mainTree.replaceAll("'", "");
    //System.out.println("main\n" + mainTree + "");
    Matcher edgeLenMatcher = Pattern.compile(":([^\\[]*)\\[([^\\]]*)\\]")
      .matcher(mainTree);
    //System.out.println("Matching");
    while (edgeLenMatcher.find()) {
      //System.out.println(edgeLenMatcher.group(1) + " main "+ edgeLenMatcher.group(2));
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
        while ((str = in .readLine()) != null) {
          jsonString.append(str);
        }
        JSONObject json = JSONObject.fromObject(jsonString.toString());

        String baseTree = trees.get(i);

        this.processJson(baseTree, json);

        fields = json.getJSONArray("fields");

      } catch (JsonSyntaxException e) {
        System.err.println(e.getLocalizedMessage());
        System.err.println("The above warnning ignored. continue ...");
      }
    }

    /*
     * Merge multiple placements for the same fragment
     */
    this.createMergedPlacements(resultsPlacements);

    if (sorted) {
      TreeSet < JSONObject > sortedPlacements = new TreeSet < JSONObject > (new Comparator < JSONObject > () {@
        Override
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
    boolean pushUp = false;
    boolean distribution = false;
    int rmUnderscore = 0;
    String mainTree = "";
    List < String > trees = new ArrayList < String > ();
    List < String > jsonLocations = new ArrayList < String > ();
    STITree < TaxonomyData > taxonomy = null;
    Hashtable < String, String > jsonNameToTaxonId = null;
    Double threshold = 0.95D;
    Double cutoff = 0D;
    FileWriter classificationWriter = null;

    /*
     * Parse optional arguments
     */
    for (int i = 3; i < args.length; i++) {
      if (args[i].equals("-s")) {
        sorted = true;
      } else if (args[i].equals("-u")) {
        pushUp = true;
      } else if (args[i].equals("-C")) {
        i++;
        cutoff = new Double(args[i]);        
      } else if (args[i].equals("-d")) {
        distribution = true;
      } else if (args[i].equals("-r")) {
        if (i + 1 >= args.length) {
          System.out.println("-r needs to be followed by a number.");
          System.exit(1);
        }
        i++;
        rmUnderscore = Integer.parseInt(args[i]);
      } else if (args[i].equals("-t")) {
        if (i + 1 >= args.length) {
          System.out.println("-t needs to be followed by a file name (taxonomy).");
          System.exit(1);
        }
        i++;
        try {
          taxonomy = TaxonomyData.readTaxonomy(args[i]);
        } catch (IOException ioe) {
          System.err.println("ERROR: Unable to read file from " + args[i]);
          return;
        }
      } else if (args[i].equals("-c")) {
        if (i + 1 >= args.length) {
          System.out.println("-c needs to be followed by a file name (classification output).");
          System.exit(1);
        }
        i++;
        try {
          classificationWriter = new FileWriter(args[i]);
        } catch (IOException ioe) {
          System.err.println("ERROR: Unable to read file from " + args[i]);
          return;
        }
      } else if (args[i].equals("-p")) {
        if (i + 1 >= args.length) {
          System.out.println("-p needs to be followd by a number.");
          System.exit(1);
        }
        i++;
        threshold = new Double(args[i]);

      } else if (args[i].equals("-m")) {
        if (i + 1 >= args.length) {
          System.out.println("-m needs to be followd by a file name (mapping between taxonomy ids and json names).");
          System.exit(1);
        }
        i++;
        try {
          jsonNameToTaxonId = TaxonomyData.readMapping(args[i]);
        } catch (IOException ioe) {
          System.err.println("ERROR: Unable to read file from " + args[i]);
          return;
        }
      }

    }

    try {
      /*
       * Read the main tree from the file
       */
      BufferedReader inm = new BufferedReader(new InputStreamReader(
        "-".equals(baseFn) ? System.in : new FileInputStream(baseFn)));
      mainTree = inm.readLine();
      if (!"-".equals(baseFn)) {
        inm.close();
      }

      /*
       * Read in json files and their associated labeled tree
       */
      if (!"-".equals(jsonDir)) {
        File[] files = new File(jsonDir).listFiles(new FilenameFilter() {
          public boolean accept(File dir, String name) {

            return ((name.indexOf(".json") >= 0) || (name.indexOf(".jplace") >= 0)) && (name.indexOf("merged") < 0);
          }
        });

        for (int i = 0; i < files.length; i++) {
          File jsonFile = files[i];
          System.out.println("Reading json " +  files[i].getAbsolutePath());
          try {
            String baseTreeFn = jsonFile.getAbsolutePath().replace("json", "labeled.tree");
            baseTreeFn = baseTreeFn.replace("jplace", "labeled.tree");
            System.out.println("Reading " + baseTreeFn);
            BufferedReader in = new BufferedReader(new InputStreamReader(
              new FileInputStream(baseTreeFn)));
            trees.add( in .readLine()); in .close();

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
        while (true) {
          line = inm.readLine();
          if (line.length() == 0) {
            break;
          }
          trees.add(line);
          jsonLocations.add(inm.readLine());
        }
      }
      //System.err.println("json locations: " + jsonLocations);


      JSONMerger merger = new JSONMerger(mainTree, trees,
        jsonLocations, sorted, rmUnderscore, taxonomy, jsonNameToTaxonId,
        threshold, classificationWriter, pushUp, distribution, cutoff);
      JSONObject merged = merger.mergeJsonFiles();
      merger.writeGSONFile(outfilename, merged);

      if (classificationWriter != null) {
        classificationWriter.close();
      }

    } catch (IOException e) {
      System.err.println("I/O Error: \n" + e.getMessage());
      System.exit(1);
    } catch (ParseException e) {
      System.err.println("Newick Parse Error: \n" + e.getMessage());
      System.exit(1);
    }
  }
}
