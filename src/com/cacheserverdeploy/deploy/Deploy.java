package com.cacheserverdeploy.deploy;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.LinkedList;


public class Deploy
{
    // maximal iteration of sub-gradient search at the root node
    private static final int ROOT_MAX_ITER = 100;
    
    // maximal iteration of sub-gradient search at the other nodes
    private static final int BRANCH_MAX_ITER = 5;
    
    // maximal iteration of sub-gradient search if backtracking happens
    private static final int BACKTRACK_MAX_ITER = 10;
    
    // time limit in milliseconds
    private static final long TIME_LIMIT = 90000L;
    
    // tolerance used to control the iteration of sub-gradient
    private static final double EPSILON = 0.01;
    
    private static final double INITIAL_LAMBDA = 2;
    
    private static final double THETA = 0.7;
    
    private static final double ETA = 1.05;
    
    private static final int STALL_LIMIT = 5;
    
    private static final double TOLERANCE = 0.000001;
    
    private static int upperBound;
    private static HashSet<Integer> bestKnown = new HashSet<Integer>();
    private static int[] bestPrice;
    private static HashSet<TreeNode> checkedNodes = new HashSet<TreeNode>();
    
    private static Network network;
    private static int numOfNodes;
    private static int numOfClients;
    private static int deviceCost;
    private static int[] supply;
    private static int[] demand;
    private static int[] mappingOfClients;
    private static int[] reverseMapping;
    
    
    public static String[] deployServer(String[] graphContent)
    {
        long startTime = System.currentTimeMillis();
        
        // get network data
        network = new Network(graphContent);
        numOfNodes = network.getNumOfNodes();
        numOfClients = network.getNumOfClients();
        deviceCost = network.getDeviceCost();
        supply = network.getSupply();
        demand = network.getDemand();
        mappingOfClients = network.getMappingOfClients();
        reverseMapping = network.getReverseMapping();
        bestPrice = new int[numOfNodes + 1];
        upperBound = deviceCost * numOfClients;
        
        // initialize the price, multipliers, sub-gradient and direction
        int[] price = new int[numOfNodes + 1];
        double[][] multipliers = new double[numOfNodes][numOfClients];
        double[][] subgradient = new double[numOfNodes][numOfClients];
        double[][] direction = new double[numOfNodes][numOfClients];
        
        int iterLimit = ROOT_MAX_ITER;
        int totalIter = 1;
        
        BranchAndBoundTree tree = new BranchAndBoundTree(numOfNodes);
        
        testLagrangean(tree.getLast());
        /*
        while (System.currentTimeMillis() - startTime <= TIME_LIMIT) {

            TreeNode curNode = tree.getLast();
            boolean needBranch = false;
            
            // set of nodes that have flow after solving minimum cost flow
            HashSet<Integer> usedNode = new HashSet<Integer>();
            
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
                                //penaltyTest(curNode, dualResult);
                                
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
                                    if (curNode.lowerBound <= upperBound - 1 && norm2(subgradient) != 0) {
                                        if (iter == 1) {
                                            initializeDirection(direction, subgradient);
                                        }
                                        
                                        double stepsize = updateMultipliers(multipliers, subgradient, direction, lambda, dualResult, iter);
                                        if (norm2(direction) <= EPSILON || stepsize <= EPSILON || iter >= iterLimit) {
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
        */
        return new String[]{"NA"};
    }
    
    private static void initializeMultipliers(double[][] multipliers, int[] price) {
        
        for (int i = 0; i < numOfNodes; i++) {
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
        for (int i = 0; i < numOfNodes; i++) {
            for (int k = 0; k < numOfClients; k++) {
                direction[i][k] = subgradient[i][k];
            }
        }
    }
    
    private static double updateMultipliers(double[][] multipliers, double[][] subgradient,
            double[][] direction, double lambda, double dualResult, int num) {
        
        if (num == 1) {
            initializeDirection(direction, subgradient);
        } else {
            for (int i = 0; i < numOfNodes; i++) {
                for (int k = 0; k < numOfClients; k++) {
                    direction[i][k] = (subgradient[i][k] + THETA * direction[i][k]) / (1 + THETA);
                }
            }
        }
        
        double estimatedUB = (ETA * upperBound + dualResult) / 2;
        //double estimatedUB = upperBound;
        double stepSize =  lambda * (estimatedUB - dualResult) / Math.pow(norm2(subgradient), 2);
        
        for (int i = 0; i < numOfNodes; i++) {
            for (int k = 0; k < numOfClients; k++) {
                multipliers[i][k] += stepSize * direction[i][k];
            }
        }
        
        return stepSize;
    }
    
    private static double norm2(double[][] data) {
        double norm = 0;
        for (int i = 0; i < numOfNodes; i++) {
            for (int k = 0; k < numOfClients; k++) {
                norm += data[i][k] * data[i][k];
            }
        }
        return Math.sqrt(norm);
    }
    
    /* dual subproblem calculation */
    private static double dualSubproblem(double[][] multipliers, TreeNode tNode, HashSet<Integer> IPlus) {
        NetworkNode[] colHead = network.getColHead();
        double dualResult = 0;
        double[] modCost = new double[numOfClients];
        for (int i = 0; i < numOfNodes; i++) {
            NetworkNode nNode = colHead[i];
            while (nNode != null) {
                for (int k = 0; k < numOfClients; k++) {
                    modCost[k] = nNode.cost + multipliers[nNode.to][k] - multipliers[i][k];
                }

                double curResult = dualSubproblemDecom(nNode, modCost);
                dualResult += curResult;
                nNode = nNode.right;
            }
        }
        
        for (int k = 0; k < numOfClients; k++) {
            dualResult += multipliers[mappingOfClients[k]][k] * supply[mappingOfClients[k]];
        }
        
        // Initialize dual ascent procedure
        double[] dualVar = new double[numOfClients];
        double[] slack = new double[numOfNodes];
        double[] tmpForSorting = new double[numOfNodes];
        double[][] multipliersMod = new double[numOfNodes + 1][numOfClients];
        
        // Sort every column of multipliers in ascending order. 
        for (int k = 0; k < numOfClients; k++ ) {
            for (int i = 0; i < numOfNodes; i++) {
                tmpForSorting[i] = multipliers[i][k] * demand[k];
            }
            Arrays.sort(tmpForSorting);
            for (int i = 0; i < numOfNodes; i++) {
                multipliersMod[i][k] = tmpForSorting[i];
            }
            // Add a high-cost dummy source
            multipliersMod[numOfNodes][k] = Double.POSITIVE_INFINITY;
        }
        
        // Initialize the dual variables and next scanned index
        int[] nextScan = new int[numOfClients];
        for (int k = 0; k < numOfClients; k++) {
            dualVar[k] = multipliersMod[0][k];
            nextScan[k] = 1;
        }
        
        // Initialize the slacks
        HashSet<Integer> closed = tNode.getClosed();
        HashSet<Integer> open = tNode.getOpen();
        for (int i = 0; i < numOfNodes; i++) {
            if (open.contains(i)) {
                slack[i] = 0;
            } else if (closed.contains(i)){
                slack[i] = Double.POSITIVE_INFINITY;
            } else {
                slack[i] = deviceCost;
            }
            for (int k = 0; k < numOfClients; k++) {
                slack[i] -= Math.max(0, dualVar[k] - multipliers[i][k] * demand[k]);
            }
        }
        
        
        IPlus.clear();
        
        double ascentResult = dualAscent(multipliers, multipliersMod, dualVar, slack, nextScan, null) + deviceCost * tNode.getOpen().size();
        while (dualAdjustment(multipliers, multipliersMod, dualVar, slack, nextScan, IPlus) > ascentResult) {
            ascentResult = dualAdjustment(multipliers, multipliersMod, dualVar, slack, nextScan, IPlus) + deviceCost * tNode.getOpen().size();
        }
        IPlus.addAll(tNode.getOpen());
        
        return dualResult + ascentResult;
    }
    
    
    private static double dualSubproblemDecom(NetworkNode nNode, double[] modCost) {
        nNode.commodityFlow = new int[numOfClients];
        nNode.modDeviceCost = 0;
        
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
            
            if (smallestModCost >= -TOLERANCE) {
                break;
            }
            
            int demand = Math.abs(supply[mappingOfClients[indexOfSmallest]]);
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
        
        if (nNode.modDeviceCost <= 0) {
            return nNode.modDeviceCost;
        } else {
            return 0;
        }
    }
    
    private static double dualAscent(double[][] multipliers, double[][] multipliersMod, 
            double[] dualVar, double[] slack, int[] nextScan, HashSet<Integer> KPlus) {
        
        // Dual ascent procedure
        while (true) {
            int delta = 0;
            for (int k = 0; k < numOfClients; k++) {
                if (KPlus == null || KPlus.contains(k)) {
                    double ascent = Double.POSITIVE_INFINITY;
                    for (int i = 0; i < numOfNodes; i++) {
                        if (dualVar[k] >= multipliers[i][k] * demand[k] - TOLERANCE && ascent > slack[i]) {
                            ascent = slack[i];
                        }
                    }
                    
                    if (ascent > multipliersMod[nextScan[k]][k] - dualVar[k] + TOLERANCE) {
                        ascent = multipliersMod[nextScan[k]][k] - dualVar[k];
                        delta = 1;
                        nextScan[k]++;
                    }
                    
                    for (int i = 0; i < numOfNodes; i++) {
                        if (dualVar[k] >= multipliers[i][k] * demand[k] - TOLERANCE) {
                            slack[i] -= ascent;
                        }
                    }
                    dualVar[k] += ascent;
                }
            }
            if (delta != 1) {
                break;
            }
        }
        
        double dualCost = 0;
        for (int k = 0; k < numOfClients; k++) {
            dualCost += dualVar[k];
        }
        
        return dualCost;
    }
    
    private static void primalSolution(double[][] multipliers, double[] dualVar, double[] slack, 
            HashSet<Integer> IStar, HashSet<Integer> IPlus) {
        
        for (int i = 0; i < numOfNodes; i++) {
            if (slack[i] >= -TOLERANCE && slack[i] <= TOLERANCE) {
                IStar.add(i);
            }
        }
        
        
        // Calculate commodity flow
        int[] iPlusk = new int[numOfClients];
        int[] iPrimek = new int[numOfClients];
        for (int k = 0; k < numOfClients; k++) {
            int smallestIndex = 0;
            double smallestCost = Double.POSITIVE_INFINITY;
            NetworkNode nNode = network.getColHead(numOfNodes);
            while (nNode != null) {
                nNode.commodityFlow[k] = 0;
                if (IPlus.contains(nNode.to) && multipliers[nNode.to][k] * demand[k] < smallestCost) {
                    smallestCost = multipliers[nNode.to][k] * demand[k];
                    smallestIndex = nNode.to;
                }
                nNode = nNode.right;
            }
            iPlusk[k] = smallestIndex;
            network.getNode(numOfNodes, smallestIndex).commodityFlow[k] = (int) multipliers[smallestIndex][k] * demand[k];
        
            nNode = network.getColHead(numOfNodes);
            double secondBest = Double.POSITIVE_INFINITY;
            while (nNode != null) {
                if (nNode.to != iPlusk[k] && IPlus.contains(nNode.to) && multipliers[nNode.to][k] * demand[k] < smallestCost) {
                    secondBest = multipliers[nNode.to][k] * demand[k];
                    smallestIndex = nNode.to;
                }
                nNode = nNode.right;
            }
            iPrimek[k] = smallestIndex;
        }
        
        //network.printFlow();
        
        // Check complementary conditions
        int[] openArray = new int[numOfNodes];
        for (int i = 0; i < numOfNodes; i++) {
            if (IPlus.contains(i)) {
                openArray[i] = 1;
            }
        }
        boolean satisfied = true;
        for (int k = 0; k < numOfClients; k++) {
            for (int i = 0; i < numOfNodes; i++) {
                if ((openArray[i] - network.getNode(numOfNodes, i).commodityFlow[k]) * Math.max(0, dualVar[k] - multipliers[i][k] * demand[k]) != 0) {
                    satisfied = false;
                    break;
                }
            }
            if (!satisfied) {
                break;
            }
        }
    }
    
    private static double dualAdjustment(double[][] multipliers, double[][] multipliersMod,
            double[] dualVar, double[] slack, int[] nextScan, HashSet<Integer> IPlus) {
        
        // Set of zero-slacks, eligible facility
        HashSet<Integer> IStar = new HashSet<Integer>();
        for (int i = 0; i < numOfNodes; i++) {
            if (slack[i] == 0) {
                IStar.add(i);
            }
        }
        
        // Set of essential facility
        for (int k = 0; k < numOfClients; k++) {
            int count = 0;
            int essentialIndex = 0;
            for (int e : IStar) {
                if (multipliers [e][k] * demand[k] <= dualVar[k]) {
                    count++;
                    essentialIndex = e;
                }
            }
            if (count == 1) {
                IPlus.add(essentialIndex);
            }
        }
        
        // Augment I+
        for (int k = 0; k < numOfClients; k++) {
            int count = 0;
            double smallestMultiplier = Double.POSITIVE_INFINITY;
            int smallestIndex = 0;
            for (int e : IPlus) {
                if (multipliers[e][k] * demand[k] <= dualVar[k]) {
                    count++;
                }
                
            }
            
            if (count == 0) {
                for (int e : IStar) {
                    if (!IPlus.contains(e) && multipliers[e][k] * demand[k] < smallestMultiplier) {
                        smallestMultiplier = multipliers[e][k] * demand[k];
                        smallestIndex = e;
                    }
                }
                IPlus.add(smallestIndex);
            }
        }
        
        
        HashSet<Integer>[] IkStar = new HashSet[numOfClients];
        HashSet<Integer>[] IkPlus = new HashSet[numOfClients];
        for (int k = 0; k < numOfClients; k++) {
            // Initialize IkStar
            IkStar[k] = new HashSet<Integer>();
            for (int e : IStar) {
                if (dualVar[k] >= multipliers[e][k] * demand[k]) {
                    IkStar[k].add(e);
                }
            }
             // Initialize IkPlus
            IkPlus[k] = new HashSet<Integer>();
            for (int e : IPlus) {
                if (dualVar[k] > multipliers[e][k] * demand[k]) {
                    IkPlus[k].add(e);
                }
            }
        }
        
        // Initialize KiPlus
        HashSet<Integer>[] KiPlus = new HashSet[numOfNodes];
        for (int i = 0; i < numOfNodes; i++) {
            KiPlus[i] = new HashSet<Integer>();
            for (int k = 0; k < numOfClients; k++) {
                if (IkStar[k].size() == 1 && IkStar[k].contains(i)) {
                    KiPlus[i].add(k);
                }
            }
        }
        
        
        for (int k = 0; k < numOfClients; k++) {
            if (IkPlus[k].size() <= 1) {
                continue;
            }
            
            // Initialize iPlusk
            int iPlusk = 0;
            double smallestCost = Double.POSITIVE_INFINITY;
            for (int e : IPlus) {
                if (multipliers[e][k] * demand[k] < smallestCost) {
                    iPlusk = e;
                    smallestCost = multipliers[e][k] * demand[k];
                }
            }
            // Initialize iPrimek
            int iPrimek = 0;
            double secondSmallestCost = Double.POSITIVE_INFINITY;
            for (int e : IkPlus[k]) {
                if (multipliers[e][k] * demand[k] != smallestCost &&
                        multipliers[e][k] * demand[k] < secondSmallestCost) {
                    iPrimek = e;
                    secondSmallestCost = multipliers[e][k] * demand[k];
                }
            }
            
            if (KiPlus[iPlusk].isEmpty() && KiPlus[iPrimek].isEmpty()) {
                continue;
            }
            
            // Initialize cMinus
            double cMinus = Double.NEGATIVE_INFINITY;
            for (int i = 0; i < numOfNodes; i++) {
                if (dualVar[k] > multipliers[i][k] * demand[k] && multipliers[i][k] * demand[k] > cMinus) {
                    cMinus = multipliers[i][k] * demand[k];
                }
            }
            
            for (int i = 0; i < numOfNodes; i++) {
                if (dualVar[k] > multipliers[i][k] * demand[k]) {
                    slack[i] += dualVar[k] - cMinus;
                }
            }
            dualVar[k] = cMinus;
            
            HashSet<Integer> KPlus = new HashSet<Integer>();
            KPlus.addAll(KiPlus[iPlusk]);
            KPlus.addAll(KiPlus[iPrimek]);
            dualAscent(multipliers, multipliersMod, dualVar, slack, nextScan, KPlus);
            KPlus.add(k);
            dualAscent(multipliers, multipliersMod, dualVar, slack, nextScan, KPlus);
            dualAscent(multipliers, multipliersMod, dualVar, slack, nextScan, null);
            
            
        }
        
        double dualResult = 0;
        for (int k = 0; k < numOfClients; k++) {
            dualResult += dualVar[k];
        }
        return dualResult;
    }
    
    
    /* sub-gradient calculation*/
    private static void calcSubgradient(double[][] subgradient) {
        
        // Clear previous results
        for (int i = 0; i < numOfNodes; i++) {
            for (int k = 0; k < numOfClients; k++) {
                subgradient[i][k] = 0;
            }
        }
        
        NetworkNode[] colHead = network.getColHead();
        NetworkNode[] rowHead = network.getRowHead();
        
        for (int i = 0; i < numOfNodes; i++) {
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
    
    private static void testLagrangean(TreeNode tNode) {
        double[][] multipliers = new double[numOfNodes][numOfClients];
        double[][] subgradient = new double[numOfNodes][numOfClients];
        double[][] direction = new double[numOfNodes][numOfClients];
        
        HashSet<Integer> IPlus = new HashSet<Integer>();
        
        double lowerBound = 0;
        
        int stall = 0;
        double lambda = INITIAL_LAMBDA;
        
        for (int i = 1; i < 100; i++) {
            double dualResult = dualSubproblem(multipliers, tNode, IPlus);
            if (dualResult < lowerBound) {
                stall++;
            } else {
                stall = 0;
                lowerBound = dualResult;
            }
            if (stall >= STALL_LIMIT) {
                lambda /= 2;
            }
            
            calcSubgradient(subgradient);
            updateMultipliers(multipliers, subgradient, direction, lambda, dualResult, i);
            System.out.println("lower bound: " + lowerBound);
        }
    }
    
    public static void main(String[] args) {
        String[] graphContent = new String[]{"4 5 2", "", "100", "", "0 1 10 2", 
                "0 2 30 5", "1 2 5 1", "1 3 15 3", "2 3 10 4", "", "0 2 20", "1 3 15"};
        deployServer(graphContent);
        
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
        //numOfNodes = 5;
        //numOfClients = 8;
        double[][] cost = {{120, 180, 100, 10000, 60, 10000, 180, 10000},
                 {210, 10000, 150, 240, 55, 210, 110, 165},
                 {180, 190, 110, 195, 50, 10000, 10000, 195},
                 {210, 190, 150, 180, 65, 120, 160, 120},
                 {170, 150, 110, 150, 70, 195, 200, 10000}};
        
        numOfNodes = 5;
        numOfClients = 8;
        demand = new int[numOfClients];
        for (int k = 0; k < numOfClients; k++) {
            demand[k] = 1;
        }
        
        int[] fixedCost = new int[] {200, 200, 200, 400, 300};
        
        HashSet<Integer> undecided = new HashSet<Integer>();
        for (int i = 0; i < numOfNodes; i++) {
            undecided.add(i);
        }
        TreeNode tNode = new TreeNode(new HashSet<Integer>(), new HashSet<Integer>(), undecided);
        // Initialize dual ascent procedure
        double[] dualVar = new double[numOfClients];
        double[] slack = new double[numOfNodes];
        double[] tmpForSorting = new double[numOfNodes];
        double[][] multipliers = cost;
        double[][] multipliersMod = new double[numOfNodes + 1][numOfClients];
        
        // Sort every column of multipliers in ascending order. 
        for (int k = 0; k < numOfClients; k++ ) {
            for (int i = 0; i < numOfNodes; i++) {
                tmpForSorting[i] = multipliers[i][k] * demand[k];
            }
            Arrays.sort(tmpForSorting);
            for (int i = 0; i < numOfNodes; i++) {
                multipliersMod[i][k] = tmpForSorting[i];
            }
            // Add a high-cost dummy source
            multipliersMod[numOfNodes][k] = Double.POSITIVE_INFINITY;
        }
        
        // Initialize the dual variables and next scanned index
        int[] nextScan = new int[numOfClients];
        for (int k = 0; k < numOfClients; k++) {
            dualVar[k] = multipliersMod[0][k];
            nextScan[k] = 1;
        }
        
        // Initialize the slacks
        HashSet<Integer> closed = tNode.getClosed();
        HashSet<Integer> open = tNode.getOpen();
        for (int i = 0; i < numOfNodes; i++) {
            if (open.contains(i)) {
                slack[i] = 0;
            } else if (closed.contains(i)){
                slack[i] = Double.POSITIVE_INFINITY;
            } else {
                slack[i] = fixedCost[i];
            }
            for (int k = 0; k < numOfClients; k++) {
                slack[i] -= Math.max(0, dualVar[k] - multipliers[i][k] * demand[k]);
            }
        }
        dualAscent(multipliers, multipliersMod, dualVar, slack, nextScan, null);
        System.out.println("Slack: " + Arrays.toString(slack));
        System.out.println("Dual Variable: " + Arrays.toString(dualVar));
        double dualResult = dualAdjustment(multipliers, multipliersMod, dualVar, slack, nextScan);
        System.out.println("Dual result: " + dualResult);
        dualResult = dualAdjustment(multipliers, multipliersMod, dualVar, slack, nextScan);
        System.out.println("Dual result: " + dualResult);
        */    
    }
}
