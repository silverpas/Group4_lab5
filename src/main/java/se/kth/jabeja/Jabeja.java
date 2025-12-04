package se.kth.jabeja;

import org.apache.log4j.Logger;
import se.kth.jabeja.config.Config;
import se.kth.jabeja.config.NodeSelectionPolicy;
import se.kth.jabeja.io.FileIO;
import se.kth.jabeja.rand.RandNoGenerator;

import java.io.File;
import java.io.IOException;
import java.util.*;

public class Jabeja {
  final static Logger logger = Logger.getLogger(Jabeja.class);
  private final Config config;
  private final HashMap<Integer/*id*/, Node/*neighbors*/> entireGraph;
  private final List<Integer> nodeIds;
  private int numberOfSwaps;
  private int round;
  // private float T; for Task1
  private double T; //Task2 required 

  private boolean resultFileCreated = false;

  //-------------------------------------------------------------------
  public Jabeja(HashMap<Integer, Node> graph, Config config) {
    this.entireGraph = graph;
    this.nodeIds = new ArrayList(entireGraph.keySet());
    this.round = 0;
    this.numberOfSwaps = 0;
    this.config = config;
    this.T = config.getTemperature();
  }


  //-------------------------------------------------------------------
  public void startJabeja() throws IOException {
    for (round = 0; round < config.getRounds(); round++) {
      for (int id : entireGraph.keySet()) {
        sampleAndSwap(id);
      }

      //one cycle for all nodes have completed.
      //reduce the temperature
      saCoolDown();
      report();
    }
  }

  /**
   * Simulated analealing cooling function
   * TASK2 version
   */
    private void saCoolDown(){
    // Task 2: exponential cooling
    double alpha = 1.0 - config.getDelta();   // es. delta=0.003 -> alpha=0.997

    if (alpha < 0.0) {
      alpha = 0.0;  // sicurezza, ma in pratica non dovrebbe succedere
    }

    T = T * alpha;

    // stesso lower bound di prima per coerenza con il Task 1
    if (T < 1.0) {
      T = 1.0;
    }
  }

  /* TASK1 version
  private void saCoolDown(){
    // TODO for second task
    if (T > 1)
      T -= config.getDelta();
    if (T < 1)
      T = 1;
  }*/

  /**
   * Sample and swap algorith at node p
   * @param nodeId
   * Try to find a swap partner for this node and swap colors.
   * The sampling depends on the node selection policy:
   * LOCAL  -> sample neighbors only
   * RANDOM -> sample random nodes in the whole graph
   * HYBRID -> first neighbors, then random nodes if needed
   */
  private void sampleAndSwap(int nodeId) {
    Node nodep = entireGraph.get(nodeId);
    Node partner = null;

    // First, try to find a partner among neighbors (LOCAL or HYBRID)
    if (config.getNodeSelectionPolicy() == NodeSelectionPolicy.HYBRID
        || config.getNodeSelectionPolicy() == NodeSelectionPolicy.LOCAL) {

      Integer[] neighbors = getNeighbors(nodep);
      partner = findPartner(nodeId, neighbors);
    }

    // If no partner was found, try random nodes in the whole graph (RANDOM or HYBRID)
    if (partner == null &&
        (config.getNodeSelectionPolicy() == NodeSelectionPolicy.HYBRID
         || config.getNodeSelectionPolicy() == NodeSelectionPolicy.RANDOM)) {

      Integer[] randomSample = getSample(nodeId);
      partner = findPartner(nodeId, randomSample);
    }

    // If a good partner exists, swap the colors
    if (partner != null) {
      int colorP = nodep.getColor();
      int colorQ = partner.getColor();

      nodep.setColor(colorQ);
      partner.setColor(colorP);

      numberOfSwaps++;
    }
  }

  /* 
   * Find the best partner for a color swap among a set of candidate nodes.
   * A partner is accepted only if the swap improves the local objective,
   * according to the JA-BE-JA formula and the current temperature T.
   * TASK2 version
   */
    public Node findPartner(int nodeId, Integer[] nodes){

    Node nodep = entireGraph.get(nodeId);
    int colorP = nodep.getColor();

    Node bestPartner = null;
    double bestNewValue = 0.0;

    double alpha = config.getAlpha();

    if (nodes == null) {
      return null;
    }

    for (Integer candidateId : nodes) {
      if (candidateId == null) {
        continue;
      }

      Node nodeq = entireGraph.get(candidateId);
      if (nodeq == null) {
        continue;
      }

      int colorQ = nodeq.getColor();

      // se hanno lo stesso colore, swap inutile
      if (colorP == colorQ) {
        continue;
      }

      // gradi con colori attuali
      int d_pp = getDegree(nodep, colorP);
      int d_qq = getDegree(nodeq, colorQ);

      // gradi se swappassimo i colori
      int d_pq = getDegree(nodep, colorQ);
      int d_qp = getDegree(nodeq, colorP);

      // valore prima dello swap
      double oldValue = Math.pow(d_pp, alpha) + Math.pow(d_qq, alpha);
      // valore dopo lo swap
      double newValue = Math.pow(d_pq, alpha) + Math.pow(d_qp, alpha);

      double delta = newValue - oldValue;
      boolean accept = false;

      if (delta > 0) {
        // miglioramento: accetto sempre
        accept = true;
      } else {
        // peggioramento: accetto con probabilità exp(delta / T)
        double prob = Math.exp(delta / T);      // delta <= 0 => 0 < prob <= 1
        double rand = RandNoGenerator.nextInt(1000000) / 1000000.0;
        if (rand < prob) {
          accept = true;
        }
      }

      // tra tutti quelli accettati, prendo quello col valore nuovo più alto
      if (accept && newValue > bestNewValue) {
        bestNewValue = newValue;
        bestPartner = nodeq;
      }
    }

    return bestPartner;
  }


  /**
   * TASK1 version
   * public Node findPartner(int nodeId, Integer[] nodes){

    Node nodep = entireGraph.get(nodeId);
    int colorP = nodep.getColor();

    Node bestPartner = null;
    double bestNewValue = 0.0;

    double alpha = config.getAlpha();

    if (nodes == null) {
      return null;
    }

    // Check each candidate node as a possible swap partner
    for (Integer candidateId : nodes) {
      if (candidateId == null) {
        continue;
      }

      Node nodeq = entireGraph.get(candidateId);
      if (nodeq == null) {
        continue;
      }

      int colorQ = nodeq.getColor();

      // Swapping with a node of the same color has no effect
      if (colorP == colorQ) {
        continue;
      }

      // Compute degrees for current colors and for swapped colors
      int d_pp = getDegree(nodep, colorP);
      int d_qq = getDegree(nodeq, colorQ);
      int d_pq = getDegree(nodep, colorQ);
      int d_qp = getDegree(nodeq, colorP);

      // Old score (before swap)
      double oldValue = Math.pow(d_pp, alpha) + Math.pow(d_qq, alpha);

      // New score (after swap)
      double newValue = Math.pow(d_pq, alpha) + Math.pow(d_qp, alpha);

      // Accept this partner if new * T > old (simulated annealing rule)
      // and keep the best one (highest newValue)
      if (newValue * T > oldValue && newValue > bestNewValue) {
        bestNewValue = newValue;
        bestPartner = nodeq;
      }
    }

    return bestPartner;
  }
  */

  /**
   * The the degreee on the node based on color
   * @param node
   * @param colorId
   * @return how many neighbors of the node have color == colorId
   */
  private int getDegree(Node node, int colorId){
    int degree = 0;
    for(int neighborId : node.getNeighbours()){
      Node neighbor = entireGraph.get(neighborId);
      if(neighbor.getColor() == colorId){
        degree++;
      }
    }
    return degree;
  }

  /**
   * Returns a uniformly random sample of the graph
   * @param currentNodeId
   * @return Returns a uniformly random sample of the graph
   */
  private Integer[] getSample(int currentNodeId) {
    int count = config.getUniformRandomSampleSize();
    int rndId;
    int size = entireGraph.size();
    ArrayList<Integer> rndIds = new ArrayList<Integer>();

    while (true) {
      rndId = nodeIds.get(RandNoGenerator.nextInt(size));
      if (rndId != currentNodeId && !rndIds.contains(rndId)) {
        rndIds.add(rndId);
        count--;
      }

      if (count == 0)
        break;
    }

    Integer[] ids = new Integer[rndIds.size()];
    return rndIds.toArray(ids);
  }

  /**
   * Get random neighbors. The number of random neighbors is controlled using
   * -closeByNeighbors command line argument which can be obtained from the config
   * using {@link Config#getRandomNeighborSampleSize()}
   * @param node
   * @return
   */
  private Integer[] getNeighbors(Node node) {
    ArrayList<Integer> list = node.getNeighbours();
    int count = config.getRandomNeighborSampleSize();
    int rndId;
    int index;
    int size = list.size();
    ArrayList<Integer> rndIds = new ArrayList<Integer>();

    if (size <= count)
      rndIds.addAll(list);
    else {
      while (true) {
        index = RandNoGenerator.nextInt(size);
        rndId = list.get(index);
        if (!rndIds.contains(rndId)) {
          rndIds.add(rndId);
          count--;
        }

        if (count == 0)
          break;
      }
    }

    Integer[] arr = new Integer[rndIds.size()];
    return rndIds.toArray(arr);
  }


  /**
   * Generate a report which is stored in a file in the output dir.
   *
   * @throws IOException
   */
  private void report() throws IOException {
    int grayLinks = 0;
    int migrations = 0; // number of nodes that have changed the initial color
    int size = entireGraph.size();

    for (int i : entireGraph.keySet()) {
      Node node = entireGraph.get(i);
      int nodeColor = node.getColor();
      ArrayList<Integer> nodeNeighbours = node.getNeighbours();

      if (nodeColor != node.getInitColor()) {
        migrations++;
      }

      if (nodeNeighbours != null) {
        for (int n : nodeNeighbours) {
          Node p = entireGraph.get(n);
          int pColor = p.getColor();

          if (nodeColor != pColor)
            grayLinks++;
        }
      }
    }

    int edgeCut = grayLinks / 2;

    logger.info("round: " + round +
            ", edge cut:" + edgeCut +
            ", swaps: " + numberOfSwaps +
            ", migrations: " + migrations);

    saveToFile(edgeCut, migrations);
  }

  private void saveToFile(int edgeCuts, int migrations) throws IOException {
    String delimiter = "\t\t";
    String outputFilePath;

    //output file name
    File inputFile = new File(config.getGraphFilePath());
    outputFilePath = config.getOutputDir() +
            File.separator +
            inputFile.getName() + "_" +
            "NS" + "_" + config.getNodeSelectionPolicy() + "_" +
            "GICP" + "_" + config.getGraphInitialColorPolicy() + "_" +
            "T" + "_" + config.getTemperature() + "_" +
            "D" + "_" + config.getDelta() + "_" +
            "RNSS" + "_" + config.getRandomNeighborSampleSize() + "_" +
            "URSS" + "_" + config.getUniformRandomSampleSize() + "_" +
            "A" + "_" + config.getAlpha() + "_" +
            "R" + "_" + config.getRounds() + ".txt";

    if (!resultFileCreated) {
      File outputDir = new File(config.getOutputDir());
      if (!outputDir.exists()) {
        if (!outputDir.mkdir()) {
          throw new IOException("Unable to create the output directory");
        }
      }
      // create folder and result file with header
      String header = "# Migration is number of nodes that have changed color.";
      header += "\n\nRound" + delimiter + "Edge-Cut" + delimiter + "Swaps" + delimiter + "Migrations" + delimiter + "Skipped" + "\n";
      FileIO.write(header, outputFilePath);
      resultFileCreated = true;
    }

    FileIO.append(round + delimiter + (edgeCuts) + delimiter + numberOfSwaps + delimiter + migrations + "\n", outputFilePath);
  }
}
