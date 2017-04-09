package com.cacheserverdeploy.deploy;

import java.util.HashSet;


/**
 * Tree node in the branch and bound tree.
 * Every tree node have three sets, closed, open and undecided.
 * Closed set contains the closed paths from the super node to all the others.
 * Open set contains the open paths from the super node to all the others.
 * Undecided set contains undecided paths from the super node to all the others.
 * @author Jiang Hao
 * @version 1.0
 */
public class TreeNode {
    private HashSet<Integer> closed = null;
    private HashSet<Integer> open = null;
    private HashSet<Integer> undecided = null;
    public double lowerBound = 0;
    
    
    public TreeNode(HashSet<Integer> closed, HashSet<Integer> open, HashSet<Integer> undecided) {
        this.closed = closed;
        this.open = open;
        this.undecided = undecided;
    }
    
    // close a node
    public void close(int index) {
        if (undecided.contains(index)) {
            undecided.remove(index);
            closed.add(index);
        }
    }
    
    // open a node
    public void open(int index) {
        if (undecided.contains(index)) {
            undecided.remove(index);
            open.add(index);
        }
    }
    
    public HashSet<Integer> getClosed() {
        return closed;
    }
    
    public HashSet<Integer> getOpen() {
        return open;
    }
    
    public HashSet<Integer> getUndecided() {
        return undecided;
    }
    
    public HashSet<Integer> getUnion() {
        HashSet<Integer> union = new HashSet<Integer>();
        union.addAll(open);
        union.addAll(undecided);
        return union;
    }
    
    public double getLowerBound() {
        return lowerBound;
    }
    
    public void setLowerBound(double lowerBound) {
        this.lowerBound = lowerBound;
    }
    
    /**
     * Deep copy of tree node
     */
    @Override
    public TreeNode clone() {
        HashSet<Integer> newClosed = new HashSet<Integer>();
        HashSet<Integer> newOpen = new HashSet<Integer>();
        HashSet<Integer> newUndecided = new HashSet<Integer>();
        newClosed.addAll(closed);
        newOpen.addAll(open);
        newUndecided.addAll(undecided);
        TreeNode newNode = new TreeNode(newClosed, newOpen, newUndecided);
        newNode.lowerBound = lowerBound;
        return newNode;
    }
    
    @Override
    public String toString() {
        String output = "closed: " + closed.toString() + "\nopen: " + open.toString() +  "\nundecided" + undecided.toString();
        return output;
        
    }
    
    public boolean isLeaf() {
        return undecided.isEmpty();
    }
    
    public static void main(String[] args) {
        HashSet<Integer> open = new HashSet<Integer>();
        HashSet<Integer> closed = new HashSet<Integer>();
        HashSet<Integer> undecided = new HashSet<Integer>();
        for (int i = 0; i < 5; i++) {
            undecided.add(i);
        }
        
        TreeNode root = new TreeNode(closed, open, undecided);
        //TreeNode newRoot = root.clone();
        
        undecided.add(100);
        root.close(1);
        root.open(2);
        //HashSet<Integer> returned = root.getUndecided();
        //HashSet<Integer> newReturned = newRoot.getUndecided();
        /*
        for (int i : returned) {
            System.out.println(i);
        }
        
        for (int i : newReturned) {
            System.out.println(i);
        }
        */
        
        HashSet<Integer> union = root.getUnion();
        for (int e : union) {
            System.out.println(e);
        }

    }

}
