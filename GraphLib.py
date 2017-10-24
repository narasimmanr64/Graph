# Util function to create a graph from a text file.
# the text file should have the number of nodes in the first line.
# the second line should contains the nodeIDs (Strings) separated by comma.
# the third line onwards, the edges from node 'a', to node 'b' with a weight of 5 will be given as 'a','b',5.0
# There is an integer value, 1 or 2 followed by the distance to indicate whether it is unidirectional or bidirectional.
# The third line onwards it looks like 'a','b',5.0,1 or 'a','b',5.0,2.  The default is 2
from __future__ import generators
import json
import numpy as np

class priorityDictionary(dict):
    def __init__(self):
        '''Initialize priorityDictionary by creating binary heap
of pairs (value,key).  Note that changing or removing a dict entry will
not remove the old pair from the heap until it is found by smallest() or
until the heap is rebuilt.'''
        self.__heap = []
        dict.__init__(self)

    def smallest(self):
        '''Find smallest item after removing deleted items from heap.'''
        # if len(self) == 0:
        #     raise IndexError, "smallest of empty priorityDictionary"
        heap = self.__heap
        while heap[0][1] not in self or self[heap[0][1]] != heap[0][0]:
            lastItem = heap.pop()
            insertionPoint = 0
            while 1:
                smallChild = 2 * insertionPoint + 1
                if smallChild + 1 < len(heap) and \
                                heap[smallChild] > heap[smallChild + 1]:
                    smallChild += 1
                if smallChild >= len(heap) or lastItem <= heap[smallChild]:
                    heap[insertionPoint] = lastItem
                    break
                heap[insertionPoint] = heap[smallChild]
                insertionPoint = smallChild
        return heap[0][1]

    def __iter__(self):
        '''Create destructive sorted iterator of priorityDictionary.'''

        def iterfn():
            while len(self) > 0:
                x = self.smallest()
                yield x
                del self[x]

        return iterfn()

    def __setitem__(self, key, val):
        '''Change value stored in dictionary and add corresponding
pair to heap.  Rebuilds the heap if the number of deleted items grows
too large, to avoid memory leakage.'''
        dict.__setitem__(self, key, val)
        heap = self.__heap
        if len(heap) > 2 * len(self):
            self.__heap = [(v, k) for k, v in self.items()]
            self.__heap.sort()  # builtin sort likely faster than O(n) heapify
        else:
            newPair = (val, key)
            insertionPoint = len(heap)
            heap.append(None)
            while insertionPoint > 0 and \
                            newPair < heap[(insertionPoint - 1) // 2]:
                heap[insertionPoint] = heap[(insertionPoint - 1) // 2]
                insertionPoint = (insertionPoint - 1) // 2
            heap[insertionPoint] = newPair

    def setdefault(self, key, val):
        '''Reimplement setdefault to call our customized __setitem__.'''
        if key not in self:
            self[key] = val
        return self[key]

class Graph:
    def __init__(self, nVertices):
        self.V = nVertices
        self.Nodes = {}
        self.RemovedNodes = {}
        self.XCoord = {}
        self.YCoord = {}
        self.Paths = []
        self.LARGENUMBER = 10000

    def __getitem__(self, item):
        return self.Nodes[item]
    def AddNodes(self, nodeList):
        for u in nodeList:
            if u not in self.Nodes:
                self.Nodes[u] = {}
    def AddEdge(self, u, v, weight, bidirectional=2):
        if u not in self.Nodes:
            self.Nodes[u] = {}
        self.Nodes[u][v] = weight
        if bidirectional == 2:
            if v not in self.Nodes:
                self.Nodes[v] = {}
            self.Nodes[v][u] = weight
    def OutEdges(self, nodeId=None):
        if nodeId is None:
            return self.Nodes
        else:
            return self.Nodes[nodeId]
    def InEdges(self, nodeId=None):
        result = {}
        if nodeId is None:
            for nodeId in self.Nodes.keys():
                r = self.InEdges(nodeId)
                result[nodeId] = r
        else:
            for AllNodes in self.Nodes.keys():
                if nodeId in self.Nodes[AllNodes]:
                    result[AllNodes] = self.Nodes[AllNodes][nodeId]
        return result

    def ModifyEdge(self, fromNode, toNode, newWeight):

        if toNode in self.Nodes[fromNode]:
            self.Nodes[fromNode][toNode] = newWeight
        else:
            raise ValueError("No Edge from %s to %s.  Call addEdge method before modifying it" %(fromNode, toNode))

    def RemoveEdge(self, fromNode, toNode, removal='hard'):
        if fromNode in self.Nodes:
            if toNode in self.Nodes[fromNode]:
                if removal == 'hard':
                    del self.Nodes[fromNode][toNode]
                else:
                    if fromNode not in self.RemovedNodes:
                        self.RemovedNodes[fromNode] = {}
                    self.RemovedNodes[fromNode][toNode] = self.Nodes[fromNode][toNode]
                    self.Nodes[fromNode][toNode] = self.LARGENUMBER
            else:
                raise ValueError("No Edge from %s to %s." %(fromNode, toNode))
        else:
            raise ValueError("No Edge from %s to %s." % (fromNode, toNode))

    def RestoreEdge(self, fromNode, toNode):
        if fromNode in self.RemovedNodes:
            if toNode in self.RemovedNodes[fromNode]:
                self.Nodes[fromNode][toNode] = self.RemovedNodes[fromNode][toNode]
                del self.RemovedNodes[fromNode][toNode]
            else:
                raise ValueError("No Edge from %s to %s." %(fromNode, toNode))
        else:
            raise ValueError("No Edge from %s to %s." % (fromNode, toNode))


    def RemoveNodes(self,nodeIds, removal='hard'):
        for nodeId in nodeIds:
            if nodeId in self.Nodes:
                inedges = self.InEdges(nodeId)
                for startNode in inedges:
                    self.RemoveEdge(startNode,nodeId, removal)
                outedges = list(self.OutEdges(nodeId))
                for endNode in outedges:
                    self.RemoveEdge(nodeId, endNode, removal)
            else:
                raise ValueError("No Node by name %s" % (nodeId))
            if removal == 'hard':
                del self.Nodes[nodeId]

    def RestoreNodes(self, nodeIds):
        for nodeId in nodeIds:
            if nodeId in self.Nodes:
                inedges = self.InEdges(nodeId)
                for startNode in inedges:
                    self.RestoreEdge(startNode,nodeId)
                outedges = list(self.OutEdges(nodeId))
                for endNode in outedges:
                    self.RestoreEdge(nodeId, endNode)
            else:
                raise ValueError("No Node by name %s" % (nodeId))

    def GetConnectedNodes(self, nodeId):
        ie = self.InEdges(nodeId)

        oe = self.OutEdges(nodeId)
        return ie,oe

    def PrintGraph(self):
        print(json.dumps(self.Nodes, sort_keys=True))

    def GetDistance(self, fromNode, toNode):
        if fromNode in self.Nodes:
            x1 = self.XCoord[fromNode]
            y1 = self.YCoord[fromNode]
        else:
            raise ValueError("No Node by name %s" %(fromNode))
        if toNode is self.Nodes:
            x2 = self.XCoord[toNode]
            y2 = self.YCoord[toNode]
        else:
            raise ValueError("No Node by name %s" %(fromNode))
        return(np.sqrt( (x1-x2)*(x1-x2) + (y1-y2)*(y1-y2)))

    def AStar(self, start, end):
        D = {}
        P = {}

        Q = priorityDictionary()
        Q[start] = 0
        for v in Q:
            D[v] = Q[v]
            if v == end:
                break
            vwLength = []
            adjNodes = []
            for w in list(self.Nodes[v]):
                if w not in D:
                    vwLength.append(D[v] + self.GetDistance(v, w) + self.GetDistance(w,end))
                    adjNodes.append(w)
            n = np.argmin(vwLength)
            w1 = adjNodes[n]

            if w1 in D:
                continue
            elif w1 not in Q or vwLength < Q[w1]:
                Q[w1] = vwLength[n]
                P[w1] = v
        return (D, P)

    def Dijkstra(self, start, end=None):

        D = {}  # dictionary of final distances
        P = {}  # dictionary of predecessors
        Q = priorityDictionary()   # estimated distances of non-final vertices
        Q[start] = 0

        for v in Q:
            D[v] = Q[v]
            if v == end:
                break
            for w in list(self.Nodes[v]):
                vwLength = D[v] + self.Nodes[v][w]
                if w in D:
                    if vwLength < D[w]:
                        raise ValueError("Dijkstra: found better path to already-final vertex")
                elif w not in Q or vwLength < Q[w]:
                    Q[w] = vwLength
                    P[w] = v

        return (D, P)


    def shortestPath(self, start, end, method='Dijkstra'):
        """
        Find a single shortest path from the given start vertex to the given
        end vertex. The input has the same conventions as Dijkstra(). The
        output is a list of the vertices in order along the shortest path.
        """
        if method == 'Dijkstra':
            D, P = self.Dijkstra(start, end)
        else:
            D, P = self.AStar(start, end)
        Path = []
        while 1:
            Path.append(end)
            if end == start:
                break
            end = P[end]
        Path.reverse()
        Length = D[Path[-1]]

        return Path, Length


    def kshortestPath(self, start, end, method='Dijkstra'):
        kpath = []
        lengths = []
        if method == 'Dijkstra':
            DPath, length = self.shortestPath(start, end)
        else:
            DPath, length = self.shortestPath(start, end, method='AStar')

        kpath.append(DPath)
        lengths.append(length)
        lena = len(DPath)
        for n in range(2, lena):
            for i in range(lena - n + 1):
                tempList = []
                for j in range(n):
                    tempList.append(DPath[i+j])
                weights = []
                for j in range(len(tempList)-1):
                    u = tempList[j]
                    v = tempList[j+1]
                    w = self.Nodes[u][v]
                    weights.append(w)
                    self.RemoveEdge(u,v)
                NewPath, NewLength = self.shortestPath(start, end)
                if not self.IsSamePath(kpath, NewPath):
                    kpath.append(NewPath)
                    lengths.append(NewLength)
                for j in range(len(tempList) - 1):
                    u = tempList[j]
                    v = tempList[j + 1]
                    self.AddEdge(u,v,weights[j],1)
        return kpath, lengths

    def IsSamePath(self, ListOfPaths, Path):
        for IndPath in ListOfPaths:
            if IndPath == Path:
                return True
        return False

    def AllPathsUtil(self, u, d, visited, path, paths, minpathlength):

        # Mark the current node as visited and store in path
        visited[u] += 1
        path.append(u)

        if u == d:
            if len(path) >= minpathlength:
                paths.append(list(path))
        else:
            # If current vertex is not destination
            # Recur for all the vertices adjacent to this vertex
            ie,oe = self.GetConnectedNodes(u)
            if oe:
                for i in oe.keys():
                    if visited[i] < 1:
                        self.AllPathsUtil(i, d, visited, path, paths, minpathlength)
            if ie:
                for i in ie.keys():
                    if visited[i] < 1:
                        self.AllPathsUtil(i, d, visited, path, paths, minpathlength)

        # Remove current vertex from path[] and mark it as unvisited
        path.pop()
        visited[u] -= 1

    # Get all paths from 's' to 'd'
    def GetAllPaths(self, s, d, minpathlength):

        # Mark all the vertices as not visited
        visited = dict.fromkeys(self.Nodes,0)

        # Create an array to store paths
        paths = []
        path = []

        # Call the recursive helper function to print all paths
        self.AllPathsUtil(s, d, visited, path, paths, minpathlength)
        return(paths)
