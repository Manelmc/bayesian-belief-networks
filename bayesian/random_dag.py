import numpy as np
import networkx as nx
from random import *
import sys

def newRandomDAG(graphSize, pLink, maxParents):
    '''
    Expects
    graphSize: The number of nodes in the graph.
    pLink: The link probability, i.e. the probability for creating an edge between two nodes. Use this parameter to control the sparseness of the graph.
    maxParents: The maximum number of parents for a child node. 

    Returns
    graph: A random directed acyclic graph as a dictionary
    '''
    graph=generateRandomGraph(graphSize,pLink)
    graph=makeDG(graph,maxParents)
    graph=makeAcyclic(graph,graphSize)

    return(graph)

def generateRandomGraph(nNodes,pLink):
# generates a random directed graph, with given size and link probability
    seed()
    randomGraph=np.zeros((nNodes,nNodes))

    for n in range(nNodes): 
        for m in range(nNodes):
            linker=uniform(0,1)
            if pLink>linker and n!=m and randomGraph[n,m]==0:
                randomGraph[m,n]=1

            if randomGraph[n,m]==1:
                switcher=randint(0,1)
                if switcher==1:
                    randomGraph[n,m]=0
                    randomGraph[m,n]=1

    return randomGraph

def makeDG(graph,k):
# turns a random graph into a directed graph
    nNodes=graph.shape[0]

    for n in range(nNodes):
        #delete parents>k
        while np.sum(graph[n,:])>k:
            nParents=np.sum(graph[n,:]) 
            deleteP=randint(0,nParents-1)
            indexParents=np.where(graph[n,: ]==1)
            graph[n,indexParents[0][deleteP]]=0

        for m in range(nNodes):
            #delete self-connections
            if n==m and graph[n,m]==1:
                graph[n,m]=0

            #delete bi-directions (a random connection will remain)
            if  graph[m,n]==1 and graph[n,m]==1:
                deleteBi=randint(0,1)
                if deleteBi==1:
                    randomGraph[n,m]=0
                    randomGraph[m,n]=1

    return graph

def makeAcyclic(graph,graphSize):
# make a given directed graph acyclic

    graphDict=makeDict(graph)
    cycles=strongly_connected_components(graphDict)

    graphDAG=dict(graphDict)
    nrPops=0
    while len(cycles)!=0:
        nrPops+=1
        seed()
        realCycles=getRealCycles(cycles)
        for c in realCycles:
            cycleLength=len(c)
            deleteCycle=randint(0,cycleLength-1)
            edges=graphDAG[c[deleteCycle]]
            edgesUpdated=(np.array(edges)).tolist()

            if deleteCycle==0 and c[cycleLength-1] in edgesUpdated:
                indexDel=edgesUpdated.index(c[cycleLength-1])
                edgesUpdated.pop(indexDel)
            elif deleteCycle!=0 and c[deleteCycle-1] in edgesUpdated:
                indexDel=edgesUpdated.index(c[deleteCycle-1])
                edgesUpdated.pop(indexDel)
            else:
                indexDel=randint(0,len(edgesUpdated)-1)
                edgesUpdated.pop(indexDel)
            graphDAG.update({c[deleteCycle]:edgesUpdated})
        cycles=strongly_connected_components(graphDAG)
        realCycles=getRealCycles(cycles)
        cycles=realCycles
    graphDAG=makeAdj(graphDAG,graphSize)
    return graphDAG

def makeAcyclic2(graph,graphSize):
# make a given directed graph acyclic and give back the indexes of effected nodes
    graphDict=makeDict(graph)
    cycles=strongly_connected_components(graphDict)
    graphDAG=dict(graphDict)
    effectedNodes=[]
    while len(cycles)!=0:
        realCycles=getRealCycles(cycles)
        for c in realCycles:
            cycleLength=len(c)
            deleteCycle=randint(0,cycleLength-1)
            edges=graphDAG[c[deleteCycle]]
            edgesUpdated=(np.array(edges)).tolist()

            if deleteCycle not in effectedNodes:
                effectedNodes.append(deleteCycle)

            if deleteCycle==0 and c[cycleLength-1] in edgesUpdated:
                indexDel=edgesUpdated.index(c[cycleLength-1])
                edgesUpdated.pop(indexDel)
            elif deleteCycle!=0 and c[deleteCycle-1] in edgesUpdated:
                indexDel=edgesUpdated.index(c[deleteCycle-1])
                edgesUpdated.pop(indexDel)
            else:
                indexDel=randint(0,len(edgesUpdated)-1)
                edgesUpdated.pop(indexDel)
            graphDAG.update({c[deleteCycle]:edgesUpdated})
        cycles=strongly_connected_components(graphDAG)
        realCycles=getRealCycles(cycles)
        cycles=realCycles
    graphDAG=makeAdj(graphDAG,graphSize)

    return graphDAG,effectedNodes

#----------------- helper functions for graph generation ---------------------#
def strongly_connected_components(graph):
    """
    Tarjan's Algorithm (named for its discoverer, Robert Tarjan) is a graph theory algorithm
    for finding the strongly connected components of a graph.
    
    Based on: http://en.wikipedia.org/wiki/Tarjan%27s_strongly_connected_components_algorithm
    """

    index_counter = [0]
    stack = []
    lowlinks = {}
    index = {}
    result = []
    
    def strongconnect(node):
        # set the depth index for this node to the smallest unused index
        index[node] = index_counter[0]
        lowlinks[node] = index_counter[0]
        index_counter[0] += 1
        stack.append(node)
    
        # Consider successors of `node`
        try:
            successors = graph[node]
        except:
            successors = []
        for successor in successors:
            if successor not in lowlinks:
                # Successor has not yet been visited; recurse on it
                strongconnect(successor)
                lowlinks[node] = min(lowlinks[node],lowlinks[successor])
            elif successor in stack:
                # the successor is in the stack and hence in the current strongly connected component (SCC)
                lowlinks[node] = min(lowlinks[node],index[successor])
        
        # If `node` is a root node, pop the stack and generate an SCC
        if lowlinks[node] == index[node]:
            connected_component = []
            
            while True:
                successor = stack.pop()
                connected_component.append(successor)
                if successor == node: break
            component = tuple(connected_component)
            # storing the result
            result.append(component)
    
    for node in graph:
        if node not in lowlinks:
            strongconnect(node)
    
    return result


def getRealCycles(cycles):
# cycles=strongly_connected_components(graphDict) gives back all strongly connected components. Need to get rid of 1-node components.
    
    indexRealCycles=[]
    
    for index,c in enumerate(cycles):
        cycleLength=len(c)
        if cycleLength>1:
            indexRealCycles.append(index)
    realCycles=[]

    for i in indexRealCycles:
        realCycles.append(cycles[i])
    
    return realCycles

def makeDict(graph):
# transform a graph into a dictionary-form, given its adjacency matrix
    nNodes=graph.shape[0]

    graphDict={}
    for n in range(nNodes):
        children=[]
        for m in range(nNodes):
            if graph[m,n]==1:
                children.append(m)
        if len(children)!=0:
            graphDict.update({n:children})

    return graphDict


def makeAdj(graph,graphSize):
# transforms a graph into adjacency-matrix-form, given dict form
    graphAdj=np.zeros((graphSize,graphSize))

    for n in graph:
        edges=graph[n]
        if len(edges)!=0:
            for e in edges:
                graphAdj[e,n]=1

    return graphAdj

def makeAdjFromNX(graph,graphSize):
    graphAdj=np.zeros((graphSize,graphSize))

    for edge in nx.edges_iter(graph):
        graphAdj[edge[1],edge[0]]=1

    return graphAdj

def generateCPT(graph, categories):
# generate a conditional probability table for a graph
    nNodes=graph.shape[0]
    CPT=[]

    for n in range(nNodes):
        nParents=int(np.sum(graph[n,:]))
        Parents=np.where(graph[n,:]==1)[0]
        categoriesN=categories[n]
        
        if nParents>0:
            for index,p in enumerate(Parents):
                if index==0:
                    categoriesP=categories[p]
                elif index!=0 and categories[p]!=0:
                    categoriesP=categoriesP*categories[p]
        else:
            categoriesP=1

        CPTn=[]
        for c in range(categoriesP):
            CPTn.append(categoricPs(categoriesN))
        
        CPT.append(CPTn)

    return CPT

def categoricPs(categories):
# generate random probabilities for a given number of categories 
    probabilities=np.zeros(categories)

    for c in range(categories):
        p=randint(1,100)
        probabilities[c]=p
    probabilities=probabilities/sum(probabilities)

    for c in range(categories):
        if c<categories-1:
            probabilities[c]=round(probabilities[c],3)
        else:
            probabilities[c]=round(1-sum(probabilities[:c]),3)

    return probabilities


#--------------------- functions for modification ----------------------------#

def deleteEdge(graph):
    newGraph=nx.DiGraph(graph)

    numberOfEdges=nx.number_of_edges(newGraph)
    edgeNumber=randint(0,numberOfEdges-1)
    edges=nx.edges_iter(newGraph)
    edgesList=[edge for edge in edges]
    newGraph.remove_edge(edgesList[edgeNumber][0],edgesList[edgeNumber][1])
    effectedNodes=np.array([edgesList[edgeNumber][1]])
    return newGraph,effectedNodes

def addEdge(graph):
    newGraph=nx.DiGraph(graph)
    numberOfNodes=nx.number_of_nodes(newGraph)
    edge=[]
    edge.append(randint(0,numberOfNodes-1))
    edge.append(randint(0,numberOfNodes-1))
    newGraph.add_edge(edge[0],edge[1])

    maxIterations=10
    i=0
    while nx.is_directed_acyclic_graph(newGraph)==False and i<maxIterations:
        i+=1
        newGraph.remove_edge(edge[0],edge[1])
        edge=[]
        edge.append(randint(0,numberOfNodes-1))
        edge.append(randint(0,numberOfNodes-1))
        newGraph.add_edge(edge[0],edge[1])

    if i==10 and nx.is_directed_acyclic_graph(newGraph)==False:
        print 'couldnt add!'
        newGraph.remove_edge(edge[0],edge[1])

    effectedNodes=np.array([edge[0]])
    return newGraph,effectedNodes

def reverseEdge(graph):
    newGraph=nx.DiGraph(graph)
    numberOfEdges=nx.number_of_edges(newGraph)
    edgeNumber=randint(0,numberOfEdges-1)
    edges=nx.edges_iter(newGraph)
    edgesList=[edge for edge in edges]
    newGraph.remove_edge(edgesList[edgeNumber][0],edgesList[edgeNumber][1])
    newGraph.add_edge(edgesList[edgeNumber][1],edgesList[edgeNumber][0])

    maxIterations=10
    i=0
    while nx.is_directed_acyclic_graph(newGraph)==False and i<maxIterations:
        i+=1
        newGraph.remove_edge(edgesList[edgeNumber][1],edgesList[edgeNumber][0])
        newGraph.add_edge(edgesList[edgeNumber][0],edgesList[edgeNumber][1])

        edgeNumber=randint(0,numberOfEdges-1)
        edges=nx.edges_iter(newGraph)
        edgesList=[edge for edge in edges]
        newGraph.remove_edge(edgesList[edgeNumber][0],edgesList[edgeNumber][1])
        newGraph.add_edge(edgesList[edgeNumber][1],edgesList[edgeNumber][0])

    if i==10 and nx.is_directed_acyclic_graph(newGraph)==False:
        print 'couldnt reverse!'
        newGraph.remove_edge(edgesList[edgeNumber][1],edgesList[edgeNumber][0])
        newGraph.add_edge(edgesList[edgeNumber][0],edgesList[edgeNumber][1])

    effectedNodes=np.array(edgesList[edgeNumber])

    return newGraph,effectedNodes

def crossoverBest(survivors,BICs,graphSize):
    crossedoverGraph=[]
    crossedoverBIC=[]
    for node in range(graphSize):
        highest=-np.inf
        index=0
        for i in range(len(survivors)):
            if BICs[i][node]>highest:
                index=i
        crossedoverGraph.append(survivors[i][node,:])
        crossedoverBIC.append(BICs[i][node])
    crossedoverGraph=np.array(crossedoverGraph)
    crossedoverGraph,effectedNodes=makeAcyclic2(crossedoverGraph,graphSize)
    return crossedoverGraph,crossedoverBIC,effectedNodes


def crossoverRandom(survivors,BICs,graphSize):
    crossedoverGraph=[]
    crossedoverBIC=[]
    for node in range(graphSize):
        survivorNumber=randint(0,len(survivors)-1)
        crossedoverGraph.append(survivors[survivorNumber][node,:])
        crossedoverBIC.append(BICs[survivorNumber][node])

    crossedoverGraph=np.array(crossedoverGraph)
    crossedoverGraph,effectedNodes=makeAcyclic2(crossedoverGraph,graphSize)
    return crossedoverGraph,crossedoverBIC,effectedNodes



def make_samples(X, bin_size):
    samples = np.empty(X.shape)
    bins = []
    for i, node in enumerate(X.T):
        samples[:,i] = np.digitize(node, np.histogram(node, bins=bin_size)[1],right=True)-1
        bins.append(np.histogram(node, bins=bin_size)[1])
    return samples, bins
