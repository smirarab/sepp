package phylolab.taxonamic;

import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.HashMap;
import java.util.Hashtable;

import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.sti.STINode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.sti.STITree;

class TaxonomyData {
    String rank;
    String id;
    String name;
    double probability;

    public TaxonomyData(String rank, String id, String name) {
	this.rank = rank;
	this.id = id;
	this.name = name;
	this.probability = 0.D;
    }

    public static STITree<TaxonomyData> readTaxonomy(String fileName) throws IOException {
	HashMap<String, STINode<TaxonomyData>> nodes = new HashMap<String, STINode<TaxonomyData>>();
	HashMap<String, String> parent = new HashMap<String, String>();
	STITree<TaxonomyData> tree = new STITree<TaxonomyData>();
	STINode<TaxonomyData> root = tree.getRoot();
	BufferedReader in = new BufferedReader(new InputStreamReader(new FileInputStream(fileName)));
	in.readLine(); // skip the header
	String line;
	boolean rootSeen = false;
	while ((line = in.readLine()) != null) {
	    String[] parts = line.split(",(?=([^\"]*\"[^\"]*\")*[^\"]*$)", -1);
	    String pid = parts[1].trim().replace("\"", ""), id = parts[0].trim().replace("\"", ""), name = parts[3].trim().replace("\"", ""), rank = parts[2].trim().replace("\"", "");
	    if (pid.equals("") || rank.equals("root")) {
		root.setName(id);
		root.setData(new TaxonomyData(rank, id, name));
		rootSeen = true;
	    } else {
		STINode<TaxonomyData> newNode, pNode;
		if ((pNode = tree.getNode(pid)) == null) {
		    pNode = root;
		    parent.put(id, pid);
		}
		newNode = pNode.createChild(id);
		newNode.setData(new TaxonomyData(rank, id, name));
		nodes.put(id, newNode);
	    }
	}
	in.close();
	if (!rootSeen) {
	    throw new RuntimeException("A root entry is expected in the taxonomy");
	}
	for (String child : parent.keySet()) {
	    nodes.get(parent.get(child)).adoptChild(nodes.get(child));
	}
	return tree;
    }

    public static Hashtable<String, String> readMapping(String fileName) throws IOException {
	Hashtable<String, String> ret = new Hashtable<String, String>();
	BufferedReader in = new BufferedReader(new InputStreamReader(new FileInputStream(fileName)));
	in.readLine(); // skip the header
	String line;
	while ((line = in.readLine()) != null) {
	    String [] parts = line.split(",");
	    ret.put(parts[0], parts[1]);
	}
	return ret ;
    }
}