package com.cacheserverdeploy.deploy;

public class NetworkNode {
    public int from;
    public int to;
    public int capacity;
    public int cost;
    //public double[] modCost;
    public double modDeviceCost = 0;
    public int flow = 0;
    public int[] commodityFlow;
    public boolean balanced = false;
    public NetworkNode right = null;
    public NetworkNode down = null;
}
