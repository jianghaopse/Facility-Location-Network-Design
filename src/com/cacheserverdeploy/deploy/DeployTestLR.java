package com.cacheserverdeploy.deploy;

import java.util.Arrays;
import java.util.HashSet;


public class DeployTestLR
{
    
    // upper bound of the total cost
    private static int upperBound;
    
    // lower bound of the total cost
    private static int lowerBound = 0;
    
    private static int totalBandCost;
    
    private static int totalDeviceCost;
    
    private static int branchNode;
    
    // maximal iteration of sub-gradient search at the root node
    private static final int ROOT_MAX_ITER = 100;
    
    // maximal iteration of sub-gradient search at the other nodes
    private static final int BRANCH_MAX_ITER = 5;
    
    // maximal iteration of sub-gradient search if backtracking happens
    private static final int BACKTRACK_MAX_ITER = 10;
    
    // time limit in milliseconds
    private static final long TIME_LIMIT = 1000L;
    
    // tolerance used to control the iteration of sub-gradient
    private static final double EPSILON = 0.01;
    
    private static double INITIAL_LAMBDA = 1.5;
    
    private static double THETA = 0;
    
    private static double ETA = 1.05;
    
    private static int STALL_LIMIT = 5;
    
    private static Network network;
    private static int numOfNodes;
    private static int numOfClients;
    private static int deviceCost;
    private static int[] supply;
    private static int[] mappingOfClients;
    
    
    
    
    public static String[] deployServer(String[] graphContent)
    {
        long startTime = System.currentTimeMillis();
        
        network = new Network(graphContent);
        numOfNodes = network.getNumOfNodes();
        numOfClients = network.getNumOfClients();
        deviceCost = network.getDeviceCost();
        supply = network.getSupply();
        mappingOfClients = network.getMappingOfClients();
        
        BranchAndBoundTree tree = new BranchAndBoundTree(numOfNodes);
        
        int[] price = new int[numOfNodes + 1];
        
        int iterLimit = ROOT_MAX_ITER;
        
        
        
        // open node of case0
        //int[] open0 = new int[] {7, 13, 15, 22, 37, 38, 43};
        
        TreeNode curNode = tree.getLast();
        /*
        HashSet<Integer> openHash = new HashSet<Integer>();
        for (int e : open0) {
            openHash.add(e);
        }
        
        for (int i = 0; i < numOfNodes + 1; i++) {
            if (!openHash.contains(i)) {
                curNode.close(i);
                continue;
            }
            curNode.open(i);
        }
        */
        //curNode.close(2);
        HashSet<Integer> usedNode = new HashSet<Integer>();
        upperBound = MinimumCostFlow.solve(network, curNode, price, usedNode) + deviceCost * numOfClients;
        upperBound = 200;
        System.out.println("UB: " + upperBound);
        
        
        
        System.out.println("******************iter1********************");
        
        double[][] multipliers = new double[numOfNodes + 1][numOfClients];
        double[][] subgradient = new double[numOfNodes + 1][numOfClients];
        
        for (int i = 0; i < numOfNodes + 1; i++) {
            for (int k = 0; k < numOfClients; k++) {
                multipliers[i][k] = price[i];
            }
        }
        
        System.out.println("Multipliers");
        for (int j = 0; j < numOfNodes + 1; j++) {
            for (int k = 0; k < numOfClients; k++) {
                System.out.printf("%8.4f", multipliers[j][k]);
            }
            System.out.println();
        }
        
        double dualResult = dualSubproblem(multipliers, curNode);
        calcSubgradient(subgradient);
        
        System.out.println("Subgradient");
        for (int j = 0; j < numOfNodes + 1; j++) {
            for (int k = 0; k < numOfClients; k++) {
                System.out.print(subgradient[j][k] + " ");
            }
            System.out.println();
        }
        
        
        double[][] direction = new double[numOfNodes + 1][numOfClients];
        for (int i = 0; i < numOfNodes + 1; i++) {
            for (int k = 0; k < numOfClients; k++) {
                direction[i][k] = subgradient[i][k];
            }
        }
        
        double lambda = INITIAL_LAMBDA;
        double bestDualResult = dualResult;
        int stall = 0;
        
        for (int i = 0; i < 20; i++) {
            
            System.out.println("******************iter" + (i + 2) + "********************");
            updateMultipliers(multipliers, subgradient, direction, lambda, dualResult);
            System.out.println("Multipliers");
            for (int j = 0; j < numOfNodes + 1; j++) {
                for (int k = 0; k < numOfClients; k++) {
                    System.out.printf("%8.4f", multipliers[j][k]);
                }
                System.out.println();
            }
            dualResult = dualSubproblem(multipliers, curNode);
            calcSubgradient(subgradient);
            if (dualResult <= bestDualResult) {
                stall++;
            } else {
                bestDualResult = dualResult;
                stall = 0;
            }
            if (stall > STALL_LIMIT) {
                lambda /= 2;
            }
            
            System.out.println("dualResult: " + dualResult);
            System.out.println("norm(gradient): " + Math.pow(norm2(subgradient), 2));
            System.out.println("lambda: " + lambda);
            
            System.out.println("Subgradient");
            for (int j = 0; j < numOfNodes + 1; j++) {
                for (int k = 0; k < numOfClients; k++) {
                    System.out.print(subgradient[j][k] + " ");
                }
                System.out.println();
            }
        }
        
        
        System.out.println("Best dual result: " + bestDualResult);
        
        return null;
    }
    
    private static void updateMultipliers(double[][] multipliers, double[][] subgradient,
            double[][] direction, double lambda, double dualResult) {
        
        for (int i = 0; i < numOfNodes + 1; i++) {
            for (int k = 0; k < numOfClients; k++) {
                direction[i][k] = (subgradient[i][k] + THETA * direction[i][k]) / (1 + THETA);
            }
        }
        //double estimatedUB = (ETA * upperBound + dualResult) / 2;
        double estimatedUB = upperBound;
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
    private static double dualSubproblem(double[][] multipliers, TreeNode tNode) {
        NetworkNode[] colHead = network.getColHead();
        double dualResult = 0;
        double[] modCost = new double[numOfClients];
        for (int i = 0; i < numOfNodes + 1; i++) {
            NetworkNode nNode = colHead[i];
            while (nNode != null) {
                for (int k = 0; k < numOfClients; k++) {
                    modCost[k] = nNode.cost + multipliers[nNode.to][k] - multipliers[i][k];
                }
                System.out.println("In dualSubProblem, modCost, " + i + "->" + nNode.to + ": "
                        + Arrays.toString(modCost));
                double curResult = dualSubproblemDecom(nNode, modCost, tNode);
                dualResult += curResult;
                
                System.out.println("In dualSubProblem, modDeviceCost, " + i + "->" + nNode.to + ": "
                        + nNode.modDeviceCost);
                System.out.println("In dualSubProblem, dualResult, " + i + "->" + nNode.to + ": "
                        + curResult);
                System.out.println("In dualSubProblem, commodity flow, " + i + "->" + nNode.to + ": "
                        + Arrays.toString(nNode.commodityFlow));
                
                nNode = nNode.right;
            }
        }
        
        for (int k = 0; k < numOfClients; k++) {
            //System.out.println("supply: " + (multipliers[mappingOfClients[k]][k] - multipliers[numOfNodes][k]) * supply[mappingOfClients[k]]);
            dualResult += (multipliers[mappingOfClients[k]][k] - multipliers[numOfNodes][k]) * supply[mappingOfClients[k]];
            //System.out.println("supply: " + supply[mappingOfClients[k]]);
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
    
    
    public static void penaltyTest() {
        
    }
    
    public static void main(String[] args) {
        String[] graphContent = new String[]{"4 5 2", "", "800", "", "0 1 10 2", 
                "0 2 30 5", "1 2 5 1", "1 3 15 3", "2 3 10 4", "", "0 2 20", "1 3 15"};
        deployServer(graphContent);
    }
}
