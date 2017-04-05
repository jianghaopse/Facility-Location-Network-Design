package com.cacheserverdeploy.deploy;

import java.util.ArrayList;
import java.util.HashSet;


/**
 * Branch and bound tree in the optimization of network deployment
 * @author Jiang Hao
 * @version 1.0
 */

public class BranchAndBoundTree {

    private ArrayList<TreeNode> tree = new ArrayList<TreeNode>();
    
    /**
     * Initialize a branch and bound tree with the root node, 
     * in which all paths are undecided.
     * @param n the total number of network nodes
     */
    public BranchAndBoundTree(int n) {
        HashSet<Integer> open = new HashSet<Integer>();
        HashSet<Integer> closed = new HashSet<Integer>();
        HashSet<Integer> undecided = new HashSet<Integer>();
        for (int i = 0; i < n; i++) {
            undecided.add(i);
        }
        TreeNode root = new TreeNode(closed, open, undecided);
        tree.add(root);
    }
    
    public boolean isEmpty() {
        return tree.isEmpty();
    }
    
    public void removeLast() {
        if (!isEmpty()) {
            tree.remove(tree.size() - 1);
        }
    }
    
    public TreeNode getLast() {
        return tree.get(tree.size() - 1);
    }
    
    public int size() {
        return tree.size();
    }
    
    /**
     * Branch on the path from the super node to node n.
     * A new branch is created by moving n from undecided set to closed or open sets.
     * After branching, the parent node is removed, and two new nodes are added.
     * @param n the path number that we want branch
     */
    public void branch(int n, double parentLB, int first) {
        TreeNode lastNode = getLast();
        TreeNode leftBranch = lastNode.clone();
        leftBranch.getClosed().add(n);
        leftBranch.getUndecided().remove(n);
        leftBranch.lowerBound = parentLB;
        TreeNode rightBranch = lastNode.clone();
        rightBranch.getOpen().add(n);
        rightBranch.getUndecided().remove(n);
        rightBranch.lowerBound = parentLB;
        removeLast();
        if (first == 1) {
            tree.add(leftBranch);
            tree.add(rightBranch);
        } else {
            tree.add(rightBranch);
            tree.add(leftBranch);
        }
        
    }
    
    public static void main(String[] args) {
        BranchAndBoundTree tree = new BranchAndBoundTree(5);
        TreeNode root = tree.getLast();
        HashSet<Integer> closedRoot = new HashSet<Integer>();
        HashSet<Integer> openRoot = new HashSet<Integer>();
        HashSet<Integer> undecidedRoot = new HashSet<Integer>();
        closedRoot = root.getClosed();
        openRoot = root.getOpen();
        undecidedRoot = root.getUndecided();
        
        System.out.println("root");
        System.out.println("Closed");
        for (int i : closedRoot) {
            System.out.println(i);
        }
        System.out.println("Open");
        for (int i : openRoot) {
            System.out.println(i);
        }
        System.out.println("Undecided");
        for (int i : undecidedRoot) {
            System.out.println(i);
        }
        
        //tree.branch(3);
        TreeNode right = tree.getLast();
        HashSet<Integer> closedRight = right.getClosed();
        HashSet<Integer> openRight = right.getOpen();
        HashSet<Integer> undecidedRight = right.getUndecided();
        
        System.out.println("Right branch");
        System.out.println("Closed");
        for (int i : closedRight) {
            System.out.println(i);
        }
        System.out.println("Open");
        for (int i : openRight) {
            System.out.println(i);
        }
        System.out.println("Undecided");
        for (int i : undecidedRight) {
            System.out.println(i);
        }
        
        tree.removeLast();
        TreeNode left = tree.getLast();
        HashSet<Integer> closedLeft = left.getClosed();
        HashSet<Integer> openLeft = left.getOpen();
        HashSet<Integer> undecidedLeft = left.getUndecided();
        
        System.out.println("Left branch");
        System.out.println("Closed");
        for (int i : closedLeft) {
            System.out.println(i);
        }
        System.out.println("Open");
        for (int i : openLeft) {
            System.out.println(i);
        }
        System.out.println("Undecided");
        for (int i : undecidedLeft) {
            System.out.println(i);
        }
        
        //tree.removeLast();
        System.out.println(tree.isEmpty());
    }
    
}
