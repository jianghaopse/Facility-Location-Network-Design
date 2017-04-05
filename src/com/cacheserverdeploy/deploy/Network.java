package com.cacheserverdeploy.deploy;

import java.util.Arrays;
import java.util.HashSet;

/**
 * Use orthogonal list to store the network
 * @author Jiang Hao
 * @version 1.0
 */
public class Network {
    private int numOfNodes;
    private int numOfArcs;
    private int numOfClients;
    private int deviceCost;
    private int totalSupply = 0;
    
    private NetworkNode[] colHead;
    private NetworkNode[] rowHead;
    private int[] supply;
    private int[] demand;
    private int[] mappingOfClients;
    private int[] reverseMapping;
    
    public Network(String[] graphContent) {
        String[] networkBasicData = graphContent[0].split(" ");
        numOfNodes = Integer.parseInt(networkBasicData[0]);
        numOfArcs = Integer.parseInt(networkBasicData[1]);
        numOfClients = Integer.parseInt(networkBasicData[2]); 
        deviceCost = Integer.parseInt(graphContent[2]);
        
        colHead = new NetworkNode[numOfNodes + 1];
        rowHead = new NetworkNode[numOfNodes + 1];
        
        for (int i = 0; i < numOfArcs; i++) {
            String[] line = graphContent[i + 4].split(" ");
            int rowIndex = Integer.parseInt(line[0]);
            int colIndex = Integer.parseInt(line[1]);
            
            // forward path (same as the original data)
            NetworkNode oldRight = colHead[rowIndex];
            NetworkNode oldDown = rowHead[colIndex];
            colHead[rowIndex] = new NetworkNode();
            colHead[rowIndex].from = rowIndex;
            colHead[rowIndex].to = colIndex;
            colHead[rowIndex].capacity = Integer.parseInt(line[2]);
            colHead[rowIndex].cost = Integer.parseInt(line[3]);
            //colHead[rowIndex].modCost = new double[numOfClients];
            colHead[rowIndex].commodityFlow = new int[numOfClients];
            colHead[rowIndex].right = oldRight;
            colHead[rowIndex].down = oldDown;
            rowHead[colIndex] = colHead[rowIndex];
            
            // undirected graph, the backward path is also stored
            oldRight = colHead[colIndex];
            oldDown = rowHead[rowIndex];
            colHead[colIndex] = new NetworkNode();
            colHead[colIndex].from = colIndex;
            colHead[colIndex].to = rowIndex;
            colHead[colIndex].capacity = Integer.parseInt(line[2]);
            colHead[colIndex].cost = Integer.parseInt(line[3]);
            //colHead[colIndex].modCost = new double[numOfClients];
            colHead[colIndex].commodityFlow = new int[numOfClients];
            colHead[colIndex].right = oldRight;
            colHead[colIndex].down = oldDown;
            rowHead[rowIndex] = colHead[colIndex];
        }
        
        // initialize the supply of nodes
        supply = new int[numOfNodes + 1];
        demand = new int[numOfClients];
        mappingOfClients = new int[numOfClients];
        reverseMapping = new int[numOfNodes];
        for (int i = 0; i < numOfClients; i++) {
            String[] line = graphContent[numOfArcs + 5 + i].split(" ");
            int curDemand = Integer.parseInt(line[2]);
            supply[Integer.parseInt(line[1])] = -curDemand;
            demand[i] = curDemand;
            totalSupply += curDemand;
            mappingOfClients[i] = Integer.parseInt(line[1]);
            reverseMapping[Integer.parseInt(line[1])] = i;
        }
        supply[numOfNodes] = totalSupply;
        
        // super node from [numOfNum] to all the real nodes
        for (int i = 0; i < numOfNodes; i++) {
            NetworkNode oldRight = colHead[numOfNodes];
            NetworkNode oldDown = rowHead[i];
            colHead[numOfNodes] = new NetworkNode();
            colHead[numOfNodes].from = numOfNodes;
            colHead[numOfNodes].to = i;
            colHead[numOfNodes].capacity = totalSupply;
            colHead[numOfNodes].cost = 0;
            colHead[numOfNodes].modDeviceCost = deviceCost;
            //colHead[numOfNodes].modCost = new double[numOfClients];
            colHead[numOfNodes].commodityFlow = new int[numOfClients];
            colHead[numOfNodes].right = oldRight;
            colHead[numOfNodes].down = oldDown;
            rowHead[i] = colHead[numOfNodes];
        }
    }
    
    public void update(TreeNode treeNode, int[] price) {
        HashSet<Integer> closed = treeNode.getClosed();
        NetworkNode node = getColHead(numOfNodes);
        while (node != null) {
            if (closed.contains(node.to)) {
                node.flow = 0;
                node.capacity = 0;
            } else {
                node.capacity = totalSupply;
            }
            node = node.right;
        }
        for (int i = 0; i < numOfNodes + 1; i++) {
            node = colHead[i];
            while (node != null) {
                // inactive
                if (price[node.from] - price[node.to] < node.cost) {
                    node.flow = 0;
                } else if (price[node.from] - price[node.to] > node.cost) {
                    //active
                    node.flow = node.capacity;
                }
                node = node.right;
            }
        }
        //printFlow();
    }
    
    public void reset() {
        for (int i = 0; i < numOfNodes + 1; i++) {
            NetworkNode nNode = colHead[i];
            while (nNode != null) {
                nNode.flow = 0;
                nNode = nNode.right;
            }
        }
    }
    
    public int getNumOfNodes() {
        return numOfNodes;
    }
    
    public int getNumOfArcs() {
        return numOfArcs;
    }
    
    public int getNumOfClients() {
        return numOfClients;
    }
    
    public int getDeviceCost() {
        return deviceCost;        
    }
    
    public int[] getSupply() {
        return supply;
    }
    
    public int[] getDemand() {
        return demand;
    }
    
    public int[] getMappingOfClients() {
        return mappingOfClients;
    }
    
    public int[] getReverseMapping() {
        return reverseMapping;
    }
    
    public NetworkNode[] getColHead() {
        return colHead;
    }
    
    public NetworkNode getColHead(int from) {
        return colHead[from];
    }
    
    public NetworkNode[] getRowHead() {
        return rowHead;
    }
    
    public NetworkNode getRowHead(int to) {
        return rowHead[to];
    }
    
    public NetworkNode getNode(int from, int to) {
        NetworkNode node = getColHead(from);
        while (node != null) {
            if (node.to == to) {
                return node;
            }
            node = node.right;
        }
        return null;
    }
    
    public void printFlow() {
        NetworkNode[] colHead = getColHead();
        for (int i = 0; i < numOfNodes + 1; i++) {
            NetworkNode node = colHead[i];
            while (node != null) {
                if (node.flow != 0) {
                    System.out.println(node.from + "->" + node.to + ": " + node.flow);
                }
                node = node.right;
            }
        }
    }
    
    public static void main(String[] args) {
        String[] graphContent = new String[]{"4 5 2", "", "100", "", "0 1 10 2", 
                "0 2 30 5", "1 2 5 1", "1 3 15 3", "2 3 10 4", "", "0 2 20", "1 3 15"};
        Network network = new Network(graphContent);
        network.printFlow();
        System.out.println(Arrays.toString(network.getSupply()));

    }
    
}


