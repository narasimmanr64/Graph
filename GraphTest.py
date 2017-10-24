import GraphLib as gl
import GraphReaderLib as grl
import timeit
import json
import numpy as np
import pandas as pd
import os

def CreateGraph():
    r = grl.GraphReader()
    nodes = r.getNodeIds('./Data/SimpleGraphNodeIDs.txt')
    n = len(nodes)
    g = gl.Graph(len(nodes))
    g.AddNodes(nodes)

    nodeDetails = r.getCoord('./Data/SimpleGraphNodeDetails.txt')

    for i in range(g.V):
        nodeName = nodeDetails['NodeId'][i]
        x = nodeDetails['XCoord'][i]
        y = nodeDetails['YCoord'][i]
        g.XCoord[nodeName] = x
        g.YCoord[nodeName] = y

    edges = r.getEdge('./Data/SimpleGraphEdges.txt')
    n = edges.shape[0]
    print(edges.iloc[0])
    for i in range(n):
        u = edges.iloc[i]['StartNodeId']
        v = edges.iloc[i]['EndNodeId']
        w = edges.iloc[i]['Distance']
        bidirect = edges.iloc[i]['BiDirectional']
        g.AddEdge(u,v,w,bidirect)
    return g

g = CreateGraph()
allPaths = g.GetAllPaths('A','G',5)
print(allPaths)
print(len(allPaths))
