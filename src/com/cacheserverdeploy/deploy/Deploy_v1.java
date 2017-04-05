package com.cacheserverdeploy.deploy;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.TreeSet;

import org.omg.Messaging.SYNC_WITH_TRANSPORT;

import com.filetool.util.FileUtil;
import com.sun.org.apache.bcel.internal.generic.NEW;
import com.sun.org.apache.xpath.internal.operations.And;
import com.sun.xml.internal.fastinfoset.algorithm.IEEE754FloatingPointEncodingAlgorithm;


public class Deploy_v1
{
    // maximal iteration of sub-gradient search at the root node
    private static final int ROOT_MAX_ITER = 100;
    
    // maximal iteration of sub-gradient search at the other nodes
    private static final int BRANCH_MAX_ITER = 5;
    
    // maximal iteration of sub-gradient search if backtracking happens
    private static final int BACKTRACK_MAX_ITER = 10;
    
    // time limit in milliseconds
    private static final long TIME_LIMIT = 88000L;
    
    // tolerance used to control the iteration of sub-gradient
    private static final double EPSILON = 0.01;
    
    private static final double INITIAL_LAMBDA = 1.5;
    
    private static final double THETA = 0.7;
    
    private static final double ETA = 1.05;
    
    private static final int STALL_LIMIT = 5;
    
    private static int upperBound;
    private static HashSet<Integer> bestKnown = new HashSet<Integer>();
    private static int[] bestPrice;
    private static HashSet<HashSet<Integer>> checked = new HashSet<HashSet<Integer>>();
    private static HashSet<TreeNode> checkedNodes = new HashSet<TreeNode>();
    
    private static Network network;
    private static int numOfNodes;
    private static int numOfClients;
    private static int deviceCost;
    private static int[] supply;
    private static int[] mappingOfClients;
    private static int[] reverseMapping;
    
    // Open nodes after the solution of dual sub-problem
    private static int[] integerInDS;
    
    
    public static String[] deployServer(String[] graphContent)
    {
        long startTime = System.currentTimeMillis();
        
        network = new Network(graphContent);
        numOfNodes = network.getNumOfNodes();
        numOfClients = network.getNumOfClients();
        deviceCost = network.getDeviceCost();
        supply = network.getSupply();
        mappingOfClients = network.getMappingOfClients();
        reverseMapping = network.getReverseMapping();
        bestPrice = new int[numOfNodes + 1];
        
        HashSet<Integer> directLinkNodes = new HashSet<Integer>();
        for (int e : mappingOfClients) {
            directLinkNodes.add(e);
        }
        //System.out.println(directLinkNodes);
        
        BranchAndBoundTree tree = new BranchAndBoundTree(numOfNodes);
        
        int[] price = new int[numOfNodes + 1];
        
        double[][] multipliers = new double[numOfNodes + 1][numOfClients];
        double[][] subgradient = new double[numOfNodes + 1][numOfClients];
        double[][] direction = new double[numOfNodes + 1][numOfClients];
        
        upperBound = deviceCost * numOfClients;
        upperBound = 110000;
        int iterLimit = ROOT_MAX_ITER;
        
        int totalIter = 1;
        
        while (System.currentTimeMillis() - startTime <= TIME_LIMIT) {
            
            TreeNode curNode = tree.getLast();
            
            //System.out.println(curNode);
            boolean needBranch = false;
            
            HashSet<Integer> union = curNode.getUnion();
            
            // set of nodes that have flow after solving minimum cost flow
            HashSet<Integer> usedNode = new HashSet<Integer>();
            
            /*
            while (union.containsAll(directLinkNodes)) {
                for (int i = 0; i < numOfNodes; i++) {
                    if (curNode.getUndecided().contains(i)) {
                        tree.branch(i, 0, 0);
                    }  
                }
                
            }
            */
            
            //System.out.println(curNode.getUndecided().containsAll(directLinkNodes));
            int bandCost = MinimumCostFlow.solve(network, curNode, price, usedNode);
            
            if (bandCost != MinimumCostFlow.INFEASIBLE) {
                if (bandCost + deviceCost * usedNode.size() < upperBound) {
                    upperBound = bandCost + deviceCost * usedNode.size();
                    bestKnown.clear();
                    bestKnown.addAll(usedNode);
                    bestPrice = Arrays.copyOf(price, numOfNodes + 1);
                }
                
                if (!curNode.getUndecided().isEmpty()) {
                    int openDeviceCost = deviceCost * curNode.getOpen().size();
                    
                    // From top to bottom, the band cost is increasing. At the current node,
                    // if the band cost plus the open device cost is greater than the upper bound,
                    // then there's no need to consider this node and its sub-nodes.
                    if  (bandCost + openDeviceCost < upperBound) {
                        double lambda = INITIAL_LAMBDA;
                        int iter = 1;
                        int stall = 0;
                        initializeMultipliers(multipliers, price);
                        
                        while (true) {
                            double dualResult = 0;
                            HashSet<Integer> openDSOdd = new HashSet<Integer>();
                            HashSet<Integer> openDSEven = new HashSet<Integer>();
                            if (iter % 2 == 1 && (dualResult = dualSubproblem(multipliers, curNode, openDSOdd)) > curNode.lowerBound) {
                                curNode.lowerBound = dualResult;
                            }
                            if (iter % 2 == 0 && (dualResult = dualSubproblem(multipliers, curNode, openDSEven)) > curNode.lowerBound) {
                                curNode.lowerBound = dualResult;
                            }
                            //System.out.println("Dual result: " + iter + " "+ dualResult);
                            //System.out.println("Lower bound: " + curNode.lowerBound);
                            //System.out.println("upper bound: " + upperBound);
                            if (curNode.lowerBound <= upperBound - 1) {
                                // If dualResult is not improved for successive STALL_LIMIT times,
                                // lambda is halved.
                                if (dualResult < curNode.lowerBound) {
                                    stall++;
                                } else {
                                    stall = 0;
                                }
                                if (stall >= STALL_LIMIT) {
                                    lambda /= 2;
                                }
                                
                                calcSubgradient(subgradient);
                                penaltyTest(curNode, dualResult);
                                
                                // After penalty test, the set can be changed, therefore,
                                // we need to check if all the nodes are fixed
                                if (!curNode.getUndecided().isEmpty()) {
                                    if (iter % 2 == 0) {
                                        // Construct a new tree node based on two successive dual sub-problems
                                        HashSet<Integer> open = new HashSet<Integer>();
                                        HashSet<Integer> closed = new HashSet<Integer>();
                                        HashSet<Integer> undecided = new HashSet<Integer>();
                                        open.addAll(openDSOdd);
                                        open.addAll(openDSEven);
                                        if (!checked.contains(open)) {
                                            // Do not need to check repeated nodes
                                            checked.add(open);
                                            //System.out.println("checked: " + checked);
                                            for (int i = 0; i < numOfNodes; i++) {
                                                if (!open.contains(i)) {
                                                    closed.add(i);
                                                }
                                            }
                                            
                                            TreeNode newNode = new TreeNode(closed, open, undecided);
                                            bandCost = MinimumCostFlow.solve(network, newNode, price, usedNode);
                                            //System.out.println(Arrays.toString(price));
                                            //network.printFlow();
                                            //System.out.println("bandCost: " + bandCost);
                                            //System.out.println("usedNode: " + usedNode + "\n");
                                            if (bandCost != MinimumCostFlow.INFEASIBLE) {
                                                if (bandCost + deviceCost * usedNode.size() < upperBound) {
                                                    upperBound = bandCost + deviceCost * usedNode.size();
                                                    bestKnown.clear();
                                                    bestKnown.addAll(usedNode);
                                                    bestPrice = Arrays.copyOf(price, numOfNodes + 1);
                                                    lambda = INITIAL_LAMBDA;
                                                }
                                            }
                                        }
                                    }
                                    if (curNode.lowerBound <= upperBound - 1 && norm2(subgradient) != 0) {
                                        if (iter == 1) {
                                            initializeDirection(direction, subgradient);
                                        }
                                        
                                        double stepsize = updateMultipliers(multipliers, subgradient, direction, lambda, dualResult);
                                        if (norm2(direction) <= EPSILON || stepsize <= EPSILON || iter >= iterLimit) {
                                            //System.out.println(norm2(direction));
                                            //System.out.println(stepsize);
                                            needBranch = true;
                                            break;
                                        }
                                        iter++;
                                    } else {
                                        needBranch = false;
                                        break;
                                    }
                                } else {
                                    bandCost = MinimumCostFlow.solve(network, curNode, price, usedNode);
                                    //System.out.println("bandCost: " + bandCost);
                                    //System.out.println("usedNode: " + usedNode + "\n");
                                    if (bandCost != MinimumCostFlow.INFEASIBLE) {
                                        if (bandCost + deviceCost * usedNode.size() < upperBound) {
                                            bestKnown.clear();
                                            bestKnown.addAll(usedNode);
                                            bestPrice = Arrays.copyOf(price, numOfNodes + 1);
                                            upperBound = bandCost + deviceCost * usedNode.size();
                                        }
                                    }
                                    needBranch = false;
                                    break;
                                }
                                
                            } else {
                                needBranch = false;
                                break;
                            }
                        }
                        checkedNodes.add(curNode);
                        if (needBranch) {
                            //System.out.println("Branch on: " + chooseBranchNode(network, curNode));
                            double parentLB = curNode.lowerBound;
                            int brachIndex = chooseBranchNode(network, curNode);
                            tree.branch(brachIndex, parentLB, 1);
                            iterLimit = BRANCH_MAX_ITER;
                            continue;
                        }
                    }
                }
                
            }
            
            if (totalIter % 1 == 0) {
                System.out.printf("%8s%8s%12s%12s%12s\n", "Iter", "Node", "Nodes left", "Best Known", "Lower bound");
                System.out.printf("%8d%8d%12d%12d%12.2f\n", totalIter, checkedNodes.size(), tree.size(), upperBound, curNode.lowerBound);
            }
            
            tree.removeLast();
            if (tree.isEmpty()) {
                break;
            }
            iterLimit = BACKTRACK_MAX_ITER;
            totalIter++;
        }
        System.out.println("Best: " + upperBound);
        System.out.println("Best node" + bestKnown);
        return outputResult();
    }
    
    private static void initializeMultipliers(double[][] multipliers, int[] price) {
        
        for (int i = 0; i < numOfNodes + 1; i++) {
            for (int k = 0; k < numOfClients; k++) {
                multipliers[i][k] = price[i];
            }
        }
        
        /*
        for (int k = 0; k < numOfClients; k++) {
            multipliers[mappingOfClients[k]][k] = price[mappingOfClients[k]];
            multipliers[numOfNodes][k] = price[mappingOfClients[k]];
        }
        */
    }
    
    private static void initializeDirection(double[][] direction, double[][] subgradient) {
        for (int i = 0; i < numOfNodes + 1; i++) {
            for (int k = 0; k < numOfClients; k++) {
                direction[i][k] = subgradient[i][k];
            }
        }
    }
    
    private static double updateMultipliers(double[][] multipliers, double[][] subgradient,
            double[][] direction, double lambda, double dualResult) {
        
        integerInDS = new int[numOfNodes];
        for (int i = 0; i < numOfNodes + 1; i++) {
            for (int k = 0; k < numOfClients; k++) {
                direction[i][k] = (subgradient[i][k] + THETA * direction[i][k]) / (1 + THETA);
            }
        }
        double estimatedUB = (ETA * upperBound + dualResult) / 2;
        //double estimatedUB = upperBound;
        double stepSize =  lambda * (estimatedUB - dualResult) / Math.pow(norm2(subgradient), 2);
        
        //System.out.print("In updateMultipliers, multipliers: ");
        for (int i = 0; i < numOfNodes + 1; i++) {
            for (int k = 0; k < numOfClients; k++) {
                multipliers[i][k] += stepSize * direction[i][k];
                //System.out.print(multipliers[i][k] + " ");
            }
            //System.out.println();
        }
        
        //System.out.println("norm(d): " + norm2(direction));
        //System.out.println("stepsize: " + stepSize);
        return stepSize;
    }
    
    private static double norm2(double[][] data) {
        double norm = 0;
        for (int i = 0; i < numOfNodes + 1; i++) {
            for (int k = 0; k < numOfClients; k++) {
                norm += data[i][k] * data[i][k];
            }
        }
        return Math.sqrt(norm);
    }
    
    /* dual subproblem calculation */
    private static double dualSubproblem(double[][] multipliers, TreeNode tNode, HashSet<Integer> openDS) {
        NetworkNode[] colHead = network.getColHead();
        double dualResult = 0;
        double[] modCost = new double[numOfClients];
        for (int i = 0; i < numOfNodes + 1; i++) {
            NetworkNode nNode = colHead[i];
            while (nNode != null) {
                for (int k = 0; k < numOfClients; k++) {
                    modCost[k] = nNode.cost + multipliers[nNode.to][k] - multipliers[i][k];
                }
                //System.out.println("In dualSubProblem, modCost, " + i + "->" + nNode.to + ": "
                //        + Arrays.toString(modCost));
                double curResult = dualSubproblemDecom(nNode, modCost, tNode);
                dualResult += curResult;
                
                //System.out.println("In dualSubProblem, modDeviceCost, " + i + "->" + nNode.to + ": "
                //        + nNode.modDeviceCost);
                //System.out.println("In dualSubProblem, dualResult, " + i + "->" + nNode.to + ": "
                //        + curResult);
                //System.out.println("In dualSubProblem, commodity flow, " + i + "->" + nNode.to + ": "
                //        + Arrays.toString(nNode.commodityFlow));
                
                nNode = nNode.right;
            }
        }
        
        for (int k = 0; k < numOfClients; k++) {
            //System.out.println("supply: " + (multipliers[mappingOfClients[k]][k] - multipliers[numOfNodes][k]) * supply[mappingOfClients[k]]);
            dualResult += (multipliers[mappingOfClients[k]][k] - multipliers[numOfNodes][k]) * supply[mappingOfClients[k]];
            //System.out.println("supply: " + supply[mappingOfClients[k]]);
        }
        
        openDS.addAll(tNode.getOpen());
        NetworkNode nNode = colHead[numOfNodes];
        while (nNode != null) {
            for (int k = 0; k < numOfClients; k++) {
                if (nNode.commodityFlow[k] != 0) {
                    openDS.add(nNode.to);
                    integerInDS[nNode.to] = 1;
                    break;
                }
            }
            nNode = nNode.right;
        }
        
        return dualResult;
    }
    
    
    private static double dualSubproblemDecom(NetworkNode nNode, double[] modCost, TreeNode tNode) {
        
        // reset modified device cost
        if (nNode.from == numOfNodes) {
            nNode.modDeviceCost = deviceCost;
        } else {
            nNode.modDeviceCost = 0;
        }
        nNode.commodityFlow = new int[numOfClients];
        
        HashSet<Integer> open = tNode.getOpen();
        HashSet<Integer> closed = tNode.getClosed();
        if (nNode.from == numOfNodes && closed.contains(nNode.to)) {
            return 0;
        }
        
        //System.out.println("In dualSubproblemDecom, modCostCopy, " + nNode.from + "->" + nNode.to + ": "
        //               + Arrays.toString(modCostCopy));
        
        int capacity = nNode.capacity;
        
        while (true) {
            double smallestModCost = modCost[0];
            int indexOfSmallest = 0;
            for (int k = 1; k < numOfClients; k++) {
                if (modCost[k] < smallestModCost) {
                    smallestModCost = modCost[k];
                    indexOfSmallest = k;
                }
            }
            
            //System.out.println("In dualSubproblemDecom, smallestModCost: " + smallestModCost);
            
            if (smallestModCost >= -0.000001) {
                break;
            }
            
            int demand = Math.abs(supply[mappingOfClients[indexOfSmallest]]);
            //System.out.println("demand: " + demand);
            if (demand < capacity) {
                nNode.commodityFlow[indexOfSmallest] = demand;
                capacity -= demand;
                nNode.modDeviceCost += modCost[indexOfSmallest] * demand;
                modCost[indexOfSmallest] = Integer.MAX_VALUE;
            } else {
                nNode.commodityFlow[indexOfSmallest] = capacity;
                nNode.modDeviceCost += modCost[indexOfSmallest] * capacity;
                break;
            }
        }
        
        if ((nNode.from == numOfNodes && open.contains(nNode.to)) || nNode.modDeviceCost <= 0) {
            //System.out.println("In dualSubproblemDecom, modDeviceCost: " + nNode.modDeviceCost);
            return nNode.modDeviceCost;
        } else {
            return 0;
        }
    }
    
    private static void dualAscent(double[][] multipliers, int[] fixedCost, TreeNode tNode, HashSet<Integer> openInDS) {
        
        double[][] multipliersCopy = new double[numOfNodes + 1][numOfClients];
        double[] tmpForSorting = new double[numOfNodes];
        HashSet<Integer> closed = tNode.getClosed();
        
        for (int k = 0; k < numOfClients; k++ ) {
            for (int i = 0; i < numOfNodes; i++) {
                if (closed.contains(i)) {
                    tmpForSorting[i] = Double.MAX_VALUE;
                } else {
                    tmpForSorting[i] = multipliers[i][k];
                }
            }
            Arrays.sort(tmpForSorting);
            for (int i = 0; i < numOfNodes; i++) {
                multipliersCopy[i][k] = tmpForSorting[i];
            }
            multipliersCopy[numOfNodes][k] = Double.POSITIVE_INFINITY;
        }
        
        double[] dualVar = new double[numOfClients];
        int[] scanned = new int[numOfClients];
        HashSet<Integer> needToScan = new HashSet<Integer>();
        for (int k = 0; k < numOfClients; k++) {
            dualVar[k] = multipliersCopy[0][k];
            scanned[k] = 1;
            needToScan.add(k);
        }
        
        double[] slack = new double[numOfNodes];
        for (int i = 0; i < numOfNodes; i++) {
            //slack[i] = deviceCost;
            slack[i] = fixedCost[i];
            for (int k = 0; k < numOfClients; k++) {
                slack[i] -= Math.max(0, dualVar[k] - multipliers[i][k]);
            }
        }
        
        while (true) {
            int iter = 0;
            int delta = 0;
            while (true) {
                if (needToScan.contains(iter)) {
                    double ascent = Double.POSITIVE_INFINITY;
                    for (int i = 0; i < numOfNodes; i++) {
                        if (dualVar[iter] >= multipliers[i][iter] && ascent > slack[i]) {
                            ascent = slack[i];
                        }
                    }
                    
                    if (ascent > multipliersCopy[scanned[iter]][iter] - dualVar[iter]) {
                        ascent = multipliersCopy[scanned[iter]][iter] - dualVar[iter];
                        delta = 1;
                        scanned[iter]++;
                    }
                    
                    for (int i = 0; i < numOfNodes; i++) {
                        if (dualVar[iter] >= multipliers[i][iter]) {
                            slack[i] -= ascent;
                        }
                    }
                    dualVar[iter] += ascent;
                }
                if (iter < numOfClients - 1) {
                    iter++;
                    continue;
                } else {
                    break;
                }
            }
            if (delta != 1) {
                break;
            }
        }
        System.out.println("Slack: " + Arrays.toString(slack));
        System.out.println("Dual vars: " + Arrays.toString(dualVar));
        double cost = 0;
        for (int i = 0; i < numOfClients; i++) {
            cost += dualVar[i];
        }
        System.out.println("Cost: " + cost);
        
        for (int i = 0; i < numOfNodes; i++) {
            boolean needOpen = false;
            if (slack[i] >= -0.000001 && slack[i] <= 0.000001) {
                for (int k = 0; k < numOfClients; k++) {
                    if (dualVar[k] >= multipliers[i][k] - 0.000001) {
                        needOpen = true;
                        break;
                    }
                }
            }
            if (needOpen) {
                openInDS.add(i);
            }
        }
        System.out.println(openInDS);
        
        // Calculate commodity flow
        int[][] x = new int[numOfNodes][numOfClients];
        for (int k = 0; k < numOfClients; k++) {
            int smallestIndex = 0;
            double smallestCost = Double.POSITIVE_INFINITY;
            for (int i = 0; i < numOfNodes; i++) {
                if (openInDS.contains(i) && multipliers[i][k] < smallestCost) {
                    smallestCost = multipliers[i][k];
                    smallestIndex = i;
                }
            }
            x[smallestIndex][k] = 1;
        }
        
        for (int[] ea : x) {
            for (int e : ea) {
                System.out.print(e + " ");
            }
            System.out.println();
        }
        
        //network.printFlow();
        
        // Check complementary conditions
        int[] openArray = new int[numOfNodes];
        for (int i = 0; i < numOfNodes; i++) {
            if (openInDS.contains(i)) {
                openArray[i] = 1;
            }
        }
        boolean satisfied = true;
        for (int k = 0; k < numOfClients; k++) {
            for (int i = 0; i < numOfNodes; i++) {
                if ((openArray[i] - x[i][k]) * Math.max(0, dualVar[k] - multipliers[i][k]) != 0) {
                    satisfied = false;
                    break;
                }
            }
            if (!satisfied) {
                break;
            }
        }
        
        System.out.println("Optimal: " + satisfied);
    }
    
    /* sub-gradient calculation*/
    private static void calcSubgradient(double[][] subgradient) {
        
        // Clear previous results
        for (int i = 0; i < numOfNodes + 1; i++) {
            for (int k = 0; k < numOfClients; k++) {
                subgradient[i][k] = 0;
            }
        }
        
        NetworkNode[] colHead = network.getColHead();
        NetworkNode[] rowHead = network.getRowHead();
        
        for (int i = 0; i < numOfNodes + 1; i++) {
            // outgoing
            NetworkNode nNode = colHead[i];
            while (nNode != null) {
                for (int k = 0; k < numOfClients; k++) {
                    subgradient[i][k] -= nNode.commodityFlow[k];
                }
                nNode = nNode.right;
            }
            
            // incoming
            nNode = rowHead[i];
            while (nNode != null) {
                for (int k = 0; k < numOfClients; k++) {
                    subgradient[i][k] += nNode.commodityFlow[k];
                }
                nNode = nNode.down;
            }
            
        }
        
        // deal with supply
        for (int k = 0; k < numOfClients; k++) {
            subgradient[mappingOfClients[k]][k] += supply[mappingOfClients[k]];
            subgradient[numOfNodes][k] -= supply[mappingOfClients[k]];
        }
    }
    
    
    private static void penaltyTest(TreeNode tNode, double dualResult) {
        HashSet<Integer> undecided = tNode.getUndecided();
        NetworkNode nNode = network.getColHead(numOfNodes);
        while (nNode != null) {
            if (undecided.contains(nNode.to)) {
                if (nNode.modDeviceCost > 0 && dualResult + nNode.modDeviceCost >= upperBound) {
                    tNode.close(nNode.to);
                } else if (nNode.modDeviceCost < 0 && dualResult + Math.abs(nNode.modDeviceCost) >= upperBound) {
                    tNode.open(nNode.to);
                }
            }
            nNode = nNode.right;
        }

    }
    
    private static int chooseBranchNode(Network network, TreeNode tNode) {
        //System.out.println("In chooseBranchNode: " + tNode);
        NetworkNode nNode = network.getColHead(numOfNodes);
        int branchIndex = -1;
        double maxModDeviceCost = -1;
        while (nNode != null) {
            if (tNode.getUndecided().contains(nNode.to) && Math.abs(nNode.modDeviceCost) > maxModDeviceCost) {
                maxModDeviceCost = Math.abs(nNode.modDeviceCost);
                branchIndex = nNode.to;
            }
            nNode = nNode.right;
        }
        return branchIndex;
    }
    
    private static String[] outputResult() {
        HashSet<Integer> closed = new HashSet<Integer>();
        HashSet<Integer> undecided = new HashSet<Integer>();
        HashSet<Integer> usedNode = new HashSet<Integer>();
        for (int i = 0; i < numOfNodes; i++) {
            if (!bestKnown.contains(i)) {
                closed.add(i);
            }
        }
        TreeNode tNode = new TreeNode(closed, bestKnown, undecided);
        int bandCost = MinimumCostFlow.solve(network, tNode, bestPrice, usedNode);
        //int totalCost = bandCost + deviceCost * bestKnown.size();
        //network.printFlow();
        
        LinkedList<String> output = new LinkedList<String>();

        NetworkNode[] colHead = network.getColHead();
        NetworkNode nNode = colHead[numOfNodes];
        ArrayList<Integer> path = new ArrayList<Integer>();
        
        int pathCount = 0;
        
        while (nNode != null) {
            int pathFlow = Integer.MAX_VALUE;
            if (nNode.flow != 0) {
                while (supply[nNode.to] == 0) {
                    path.add(nNode.to);
                    if (nNode.flow < pathFlow) {
                        pathFlow = nNode.flow;
                    }
                    nNode = colHead[nNode.to];
                    while (nNode.flow == 0) {
                        nNode = nNode.right;
                    }
                    if (nNode.flow < pathFlow) {
                        pathFlow = nNode.flow;
                    }
                }
                path.add(nNode.to);
                path.add(reverseMapping[nNode.to]);
                
                if (-supply[nNode.to] < pathFlow) {
                    pathFlow = -supply[nNode.to];
                }
                
                network.getNode(numOfNodes, path.get(0)).flow -= pathFlow;
                for (int i = 0; i < path.size() - 2; i++) {
                    network.getNode(path.get(i), path.get(i + 1)).flow -= pathFlow;
                }
                supply[path.get(path.size() - 2)] += pathFlow;
                path.add(pathFlow);
                pathCount++;
                
                StringBuilder sb = new StringBuilder();
                for (int e : path) {
                    sb.append(e).append(" ");
                }
                output.add(sb.toString());
                
                path.clear();
                nNode = colHead[numOfNodes];
            } else {
                nNode = nNode.right;
            }
        }
        output.add(0, Integer.toString(pathCount));
        output.add(1, "");
        
        return (String[]) output.toArray(new String[0]);
    }
    
    public static void main(String[] args) {
        //String[] graphContent = new String[]{"4 5 2", "", "100", "", "0 1 10 2", 
        //        "0 2 30 5", "1 2 5 1", "1 3 15 3", "2 3 10 4", "", "0 2 20", "1 3 15"};
        //deployServer(graphContent);
        
        /*
        HashSet<Integer> closed = new HashSet<Integer>();
        HashSet<Integer> open = new HashSet<Integer>();
        HashSet<Integer> undecided = new HashSet<Integer>();
        HashSet<Integer> usedNode = new HashSet<Integer>();
        closed.add(1);
        closed.add(3);
        undecided.add(0);
        undecided.add(2);
        TreeNode tNode = new TreeNode(closed, open, undecided);
        network = new Network(graphContent);
        
        System.out.println(MinimumCostFlow.solve(network, tNode, new int[network.getNumOfNodes() + 1], usedNode));
        network.printFlow();
        System.out.println(usedNode);
        */
        
        
        /*** Test dual ascent ***/
        /*
        numOfNodes = 5;
        numOfClients = 8;
        double[][] cost = {{120, 180, 100, 10000, 60, 10000, 180, 10000},
                 {210, 10000, 150, 240, 55, 210, 110, 165},
                 {180, 190, 110, 195, 50, 10000, 10000, 195},
                 {210, 190, 150, 180, 65, 120, 160, 120},
                 {170, 150, 110, 150, 70, 195, 200, 10000}};
        
        for (int i = 0; i < numOfNodes; i++) {
            for (int k = 0; k < numOfClients; k++) {
                cost[i][k] = 0;
            }
        }
        int[] fixedCost = new int[] {50, 50, 50, 50, 50};
        
        HashSet<Integer> open = new HashSet<Integer>();
        HashSet<Integer> closed = new HashSet<Integer>();
        HashSet<Integer> undecided = new HashSet<Integer>();
        HashSet<Integer> openInDS = new HashSet<Integer>();
        dualAscent(cost, fixedCost, new TreeNode(closed, open, undecided), openInDS);
        */
        
        String[] graphContent = FileUtil.read("case4.txt", null);
        Network network = new Network(graphContent);
        HashSet<Integer> closed = new HashSet<Integer>();
        HashSet<Integer> open = new HashSet<Integer>();
        HashSet<Integer> undecided = new HashSet<Integer>();
        HashSet<Integer> usedNode = new HashSet<Integer>();
        
        int[] openArray = new int[] {48, 20, 37, 22, 26, 12, 15};
        for (int i = 0; i < openArray.length; i++) {
            open.add(openArray[i]);
        }
        
        for (int i = 0; i < 50; i++) {
            if (!open.contains(i)) {
                closed.add(i);
            }
        }
        
        int[] price = new int[51];
        
        TreeNode tNode = new TreeNode(closed, open, undecided);
        System.out.println(tNode);
        int bandCost = MinimumCostFlow.solve(network, tNode, price, usedNode);
        int totalCost = bandCost + network.getDeviceCost() * open.size(); 
        System.out.println(totalCost);
        network.printFlow();
        outputResult();
    }
}
 