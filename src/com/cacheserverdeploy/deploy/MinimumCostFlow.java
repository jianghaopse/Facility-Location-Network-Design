package com.cacheserverdeploy.deploy;

import java.util.HashSet;


public class MinimumCostFlow {
    public static final int INFEASIBLE = -1;
    public static final int FEASIBLE = -2;
    public static final int OPTIMAL = -3;
    public static final int NONOPTIMAL = -4;
    public static final int NO_AUGMENT_PATH = -5;
    
    // surplus of a node, surplus = incoming - outgoing + supply
    private static int[] surplus;
    private static HashSet<Integer> labeled = new HashSet<Integer>();
    private static HashSet<Integer> scanned = new HashSet<Integer>();
    private static HashSet<Integer> newAdded = new HashSet<Integer>();
    private static int[][] augmentPath;
    
    
    public static int solve(Network network, TreeNode treeNode, int[] price, HashSet<Integer> usedNode) {
        // Set the capacity of the closed path to be 0, 
        // and modify the flow to satisfy complementary condition
        network.update(treeNode, price);
        
        int numOfNodes = network.getNumOfNodes();
        surplus = new int[numOfNodes + 1];
        newAdded.clear();
        
        augmentPath = new int[numOfNodes + 1][2];
        
        while (true) {
            
            while (true) {
                scanned.clear();
                labeled.clear();
                int indexOfPosSurplus = 0;
                if ((indexOfPosSurplus = getPosSurplus(network)) == INFEASIBLE) {
                    return INFEASIBLE;
                } else if (indexOfPosSurplus == OPTIMAL) {
                    usedNode.clear();
                    //usedNode.addAll(treeNode.getOpen());
                    NetworkNode nNode = network.getColHead(numOfNodes);
                    while (nNode != null) {
                        if (nNode.flow > 0) {
                            usedNode.add(nNode.to);
                        }
                        nNode = nNode.right;
                    }
                    return calcCost(network);
                } else {
                    labeled.add(indexOfPosSurplus);
                    
                    int indexOfNegSurplus = 0;
                    int ascent = 0;
                    boolean skipAug = false;
                    
                    do {
                        if (!scanned.equals(labeled)) {
                            int next = nextScanned();
                            
                            //System.out.println("scanned: " + scanned.toString());
                            ascent += labelNeighbor(network, price, next);
                            //System.out.println("ascent: " + ascent);
                            scanned.add(next);
                            if (ascent > 0) {
                                skipAug = true;
                            }
                        } else {
                            skipAug = true;
                        }
                    } while (!skipAug && (indexOfNegSurplus = getNegSurplus(network)) == NO_AUGMENT_PATH);
                    if (skipAug) {
                        break;
                    }
                    //System.out.println("neg" + indexOfNegSurplus);
                    flowAugmentation(network, indexOfPosSurplus, indexOfNegSurplus);
                } 
            }
            
            if (priceChange(network, price) == INFEASIBLE) {
                return INFEASIBLE;
            }
        }
        
    }
    
    private static int calcSurplus(Network network, int i) {
        int incoming = 0;
        int outgoing = 0;
        // incoming
        NetworkNode node = network.getRowHead(i);
        int[] supply = network.getSupply();
        while (node != null) {
            incoming += node.flow;
            node = node.down;
        }
        // outgoing
        node = network.getColHead(i);
        while (node != null) {
            outgoing += node.flow;
            node = node.right;
        }
        surplus[i] = incoming - outgoing + supply[i];
        return surplus[i];
    }
    
    private static int getPosSurplus(Network network) {
        int numOfNodes = network.getNumOfNodes();
        for (int i = numOfNodes; i >= 0; i--) {
            //System.out.println(i + "surplus:" + calcSurplus(network, surplus, i));
            if (calcSurplus(network, i) > 0) {
                return i;
            }
        }
        for (int e : surplus) {
            if (surplus[e] < 0) {
                return INFEASIBLE;
            }
        }
        return OPTIMAL;
    }
    
    
    private static int nextScanned() {
        for (int i : labeled) {
            if (!scanned.contains(i)) {
                return i;
            }
        }
        // this is impossible
        return 0;
    }
    
    private static int labelNeighbor(Network network, int[] price, int next) {
        int ascentIncrement = 0;
        //backward path
        NetworkNode node = network.getRowHead(next);
        while (node != null) {
            if (!labeled.contains(node.from) && price[node.from] - price[next] == node.cost
                    && node.flow > 0) {
                newAdded.add(node.from);
                augmentPath[node.from][0] = next;
                augmentPath[node.from][1] = 0;
            }
            if (scanned.contains(node.from) && price[node.from] - price[next] == node.cost) {
                ascentIncrement += node.capacity - node.flow;
            }
            if (!scanned.contains(node.from) && price[node.from] - price[next] == node.cost) {
                ascentIncrement -= node.flow;
            }
            
            node = node.down;
        }
        
        // forward path
        node = network.getColHead(next);
        while (node != null) {
            if (!labeled.contains(node.to) && price[next] - price[node.to] == node.cost
                    && node.flow < node.capacity) {
                newAdded.add(node.to);
                augmentPath[node.to][0] = next;
                augmentPath[node.to][1] = 1;
            }
            if (scanned.contains(node.to) && price[next] - price[node.to] == node.cost) {
                ascentIncrement += node.flow;
            }
            if (!scanned.contains(node.to) && price[next] - price[node.to] == node.cost) {
                ascentIncrement -= node.capacity - node.flow;
            }
            node = node.right;
        }
        labeled.addAll(newAdded);
        
        return ascentIncrement;
    }
    
    private static int getNegSurplus(Network network) {
        for (int e : newAdded) {
            if (calcSurplus(network, e) < 0) {
                newAdded.clear();
                return e;
            }
        }
        newAdded.clear();
        return NO_AUGMENT_PATH;
    }
    
    private static void flowAugmentation(Network network, int indexOfPosSurplus, int indexOfNegSurplus) {
        
        int delta = Math.min(surplus[indexOfPosSurplus], -surplus[indexOfNegSurplus]);
        int i = indexOfNegSurplus;
        
        // find the largest flow that can push
        while (true) {
            // backward path
            if (augmentPath[i][1] == 0) {
                int slack = network.getNode(i, augmentPath[i][0]).flow;
                delta = delta < slack? delta : slack;
            } else {
                // forward path
                int slack = network.getNode(augmentPath[i][0], i).capacity - network.getNode(augmentPath[i][0], i).flow;
                delta = delta < slack? delta : slack;
            }
            if (augmentPath[i][0] == indexOfPosSurplus) {
                break;
            }
            i = augmentPath[i][0];
        }
        //System.out.println("In flowAugmentation, delta: " + delta);
        
        // push the flow
        i = indexOfNegSurplus;
        while (true) {
            // backward path
            if (augmentPath[i][1] == 0) {
                network.getNode(i, augmentPath[i][0]).flow -= delta;
            } else {
                // forward path
                network.getNode(augmentPath[i][0], i).flow += delta;
            }
            if (augmentPath[i][0] == indexOfPosSurplus) {
                break;
            }
            i = augmentPath[i][0];
        }
    }
    
    private static int priceChange(Network network, int[] price) {
        boolean updated = false;
        for (int e : scanned) {
            // outgoing
            NetworkNode node = network.getColHead(e);
            //System.out.println(Arrays.toString(price));
            while (node != null) {
                if (!scanned.contains(node.to) && price[e] - price[node.to] == node.cost) {
                    node.flow = node.capacity;
                }
                node = node.right;
            }
            // incoming
            node = network.getRowHead(e);
            while (node != null) {
                if (!scanned.contains(node.from) && price[node.from] - price[e] == node.cost) {
                    node.flow = 0;
                }
                node = node.down;
            }
        }

        int gamma = Integer.MAX_VALUE;
        for (int e : scanned) {
            // outgoing
            NetworkNode node = network.getColHead(e);
            while (node != null) {
                if (!scanned.contains(node.to) && node.flow < node.capacity) {
                    gamma = Math.min(gamma, price[node.to] + node.cost - price[e]);
                    updated = true;
                }
                node = node.right;
            }
            // incoming
            node = network.getRowHead(e);
            while (node != null) {
                if ( !scanned.contains(node.from) && node.flow > 0) {
                    gamma = Math.min(gamma, price[node.from] - node.cost - price[e]);
                    updated = true;
                }
                node = node.down;
            }
        }
        
        if (!updated) {
            return INFEASIBLE;
        }
        for (int e : scanned) {
            price[e] += gamma;
        }
        
        //System.out.println("In priceChange: network");
        //network.printFlow();
        //System.out.println("In priceChange: " + gamma);
        return FEASIBLE;
    }
    
    private static int calcCost(Network network) {
        NetworkNode[] colHead = network.getColHead();
        int totalBandCost = 0;
        // omit the last row, because the band cost is 0
        for (int i = 0; i < colHead.length - 1; i++) {
            NetworkNode node = colHead[i];
            while (node != null) {
                totalBandCost += node.cost * node.flow;
                node = node.right;
            }
        }
        return totalBandCost;
    }
    
    public static void main(String[] args) {
        String[] graphContent = new String[]{"4 5 2", "", "100", "", "0 1 10 2", 
                "0 2 30 5", "1 2 5 1", "1 3 15 3", "2 3 10 4", "", "0 2 20", "1 3 15"};
        Network network = new Network(graphContent);
        int numOfNodes = network.getNumOfNodes();
        
        BranchAndBoundTree tree = new BranchAndBoundTree(numOfNodes);
        TreeNode root = tree.getLast();
        root.close(0);
        root.close(2);
        
        int[] price = new int[numOfNodes + 1];
        price = new int[] {0, 0, 0, 0, 0};
        HashSet<Integer> usedNode = new HashSet<Integer>();
        System.out.println(MinimumCostFlow.solve(network, root, price, usedNode) + network.getDeviceCost() * usedNode.size());
        network.printFlow();
    }

}
