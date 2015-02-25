import numpy as np
import networkx as nx
from random import *
import sys
import math

def newRandomDAG(graphSize, pLink, maxParents):
    '''
    Expects
    graphSize: The number of nodes in the graph.
    pLink: The link probability, i.e. the probability for creating an edge between two nodes. Use this parameter to control the sparseness of the graph.
    maxParents: The maximum number of parents for a child node. 

    Returns
    graph: A random directed acyclic graph in adjacency matrix form
    '''
    graph=generateRandomGraph(graphSize,pLink)
    graph=makeDG(graph,maxParents)
    graph, _=make_acyclic(graph,graphSize)

    return(graph)

def generateRandomGraph(nNodes,pLink):
# generates a random directed graph, with given size and link probability
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

def make_acyclic(graph, graphSize):
# make a given directed graph acyclic and give back the indexes of affected nodes
    graphDict=makeDict(graph)
    cycles=strongly_connected_components(graphDict)
    graphDAG=dict(graphDict)
    affected_nodes=[]
    print graphDict
    while len(cycles)!=0:
        realCycles=getRealCycles(cycles)
        print realCycles
        for c in realCycles:
            print c
            cycleLength=len(c)
            deleteCycle=randint(0,cycleLength-1)
            print deleteCycle
            edges=graphDAG[c[deleteCycle]]
            print edges
            exit()
            edgesUpdated=(np.array(edges)).tolist()

            if deleteCycle not in affected_nodes:
                affected_nodes.append(deleteCycle)

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

    return graphDAG, affected_nodes

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
#--------------------- functions for modification ----------------------------#

def deleteEdge(graph):
    newGraph=nx.DiGraph(graph)
    numberOfEdges=nx.number_of_edges(newGraph)
    assert numberOfEdges is not 0, 'There is no causal relation between the variables.'
    edgeNumber=randint(0,numberOfEdges-1)
    edges=nx.edges_iter(newGraph)
    edgesList=[edge for edge in edges]
    newGraph.remove_edge(edgesList[edgeNumber][0],edgesList[edgeNumber][1])
    # only BIC for child node changes. 
    affected_nodes=np.array([edgesList[edgeNumber][1]])
    print 'affected nodes del: ',affected_nodes
    return newGraph,affected_nodes

def addEdge(graph):
    '''Add an edge to the graph. If the new edge makes the graph cyclic remove the edge again. Try this procedure up to 10 times.'''
    newGraph=nx.DiGraph(graph)
    numberOfNodes=nx.number_of_nodes(newGraph)
    assert numberOfNodes is not 0, 'There is no causal relation between the variables.'

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

    # only BIC for parent node changes. 
    affected_nodes=np.array([edge[0]])
    print 'affected nodes add: ',affected_nodes
    return newGraph,affected_nodes

def reverseEdge(graph):
    '''Reverse an edge in the graph. If the new edge makes the graph cyclic remove the edge again. Try this procedure up to 10 times.'''
    newGraph=nx.DiGraph(graph)
    numberOfEdges=nx.number_of_edges(newGraph)
    assert numberOfEdges is not 0, 'There is no causal relation between the variables.'
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
    # BIC for parent and child nodes change. 
    affected_nodes=np.array(edgesList[edgeNumber])
    print 'affected nodes reverse: ',affected_nodes

    return newGraph,affected_nodes

def crossoverBest(survivors,BICs,graphSize):
    '''Crossover the nodes that yield the highest BIC amongst the fittest 
    individual graphs'''
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
    crossedoverGraph,affected_nodes=make_acyclic(crossedoverGraph,graphSize)
    return crossedoverGraph,crossedoverBIC,affected_nodes


def crossoverRandom(survivors,BICs,graphSize):
    '''Crossover survivors at random'''
    crossedoverGraph=[]
    crossedoverBIC=[]
    for node in range(graphSize):
        survivorNumber=randint(0,len(survivors)-1)
        crossedoverGraph.append(survivors[survivorNumber][node,:])
        crossedoverBIC.append(BICs[survivorNumber][node])

    crossedoverGraph=np.array(crossedoverGraph)
    crossedoverGraph,affected_nodes=make_acyclic(crossedoverGraph,graphSize)
    return crossedoverGraph,crossedoverBIC,affected_nodes

def make_samples(X, bin_size):
    '''Create samples from data by discretizing over specified bins. 
    Return the samples in form of a dict and the bins for each node'''
    samples = np.empty(X.shape)
    bins = []
    for i, node in enumerate(X.T):
        samples[:,i] = np.digitize(node, np.histogram(node, bins=bin_size)[1],right=True)-1
        bins.append(np.histogram(node, bins=bin_size)[1])
    return samples, bins

def getCombinations(graph, samples, categories, affected_nodes):

    combinationsList={}
    numberOfCombinationsList={}
        
    for node in affected_nodes:
        Parents=np.where(graph[node,:]==1)[0]
        #get the number of combinations by computing the product of categories in the parents
        numberOfCombinations = np.prod(categories[Parents])
        numberOfCombinationsList.update({node:numberOfCombinations})
        numberOfCombinationsIncludingNode=numberOfCombinations*categories[node]
        combinationCountsNode=np.zeros(numberOfCombinationsIncludingNode)

        for samplePieceN in samples:
            #get the combinationNr for one node in sample n
            combinationNr=0
            for index,parent in enumerate(Parents):
                if index==0:
                    combinationNr=samplePieceN[parent]*categories[node]
                    catParentOld=categories[parent]*categories[node]
                else:
                    combinationNr+=samplePieceN[parent]*catParentOld
                    catParentOld*=categories[parent]

            combinationNr=combinationNr+samplePieceN[node]
            combinationCountsNode[combinationNr]+=1

        combinationsList.update({node:combinationCountsNode})

    return combinationsList,numberOfCombinationsList


def computeBIC(estGraph, samples, categories, affected_nodes=[], BICold=[]):
# compute the Bayesian Information Criterion for a graph, given data
    if affected_nodes==[]:
        affected_nodes=np.arange(0,len(estGraph))
        BICnew=np.zeros(len(affected_nodes))
    else:
        BICnew=np.copy(BICold)


    Nijk,qi=getCombinations(estGraph,samples,categories,affected_nodes)
    Nij={}
    for node in Nijk:
        nodeNij=[]
        index=0
        for i in range(qi[node]):
            nodeNij.append(sum(Nijk[node][index:index+categories[node]]))
            index+=categories[node]
        Nij.update({node:nodeNij})

    LL={}
    for node in affected_nodes:
        index=0
        LLNode=0
        for combinationCount in range(len(Nij[node])):
            for categoryCount in range(categories[node]):
                if Nij[node][combinationCount]==0:
                    Nij[node][combinationCount]=1
                if Nijk[node][categoryCount+index]==0:
                    Nijk[node][categoryCount+index]=1
                LLNode+=Nijk[node][categoryCount+index]*math.log(Nijk[node][categoryCount+index]/Nij[node][combinationCount])
            index+=categories[node]
        LL.update({node:LLNode})

    categories=np.array(categories)
    qi=np.array(qi.items())[:,1]    

    complexityPenalty=0.5*math.log(len(samples))*abs((categories[affected_nodes]-1)*qi)
    for index,node in enumerate(affected_nodes):
        BICnew[node]=LL[node]-complexityPenalty[index]

    return BICnew

def hill_climber(self, max_iter):
    bestGraphHC=makeDict(self.graph)
    bestGraphHC = nx.DiGraph(bestGraphHC)
    bestBIC=self.BIC

    for i in range(max_iter):
        mode = randint(1,3)
        if mode==1:
            graphHC,affected_nodes=deleteEdge(bestGraphHC)
        elif mode==2:
            graphHC,affected_nodes=addEdge(bestGraphHC)
        elif mode==3:
            graphHC,affected_nodes=reverseEdge(bestGraphHC)

        graphHC=makeAdjFromNX(graphHC,self.graph_size)
        newBIC=computeBIC(graphHC, self.samples, self.categories, affected_nodes, bestBIC)

        if sum(newBIC)>sum(bestBIC):
            bestGraphHC=makeDict(graphHC)
            bestGraphHC=nx.DiGraph(bestGraphHC)
            bestBIC=newBIC[:]

    self.graph = makeAdjFromNX(bestGraphHC, self.graph_size)
    self.BIC = sum(bestBIC)
    print self.BIC
    print self.graph

def simulatedAnnealing(self, startingTemp, alpha, max_iter):
    bestGraphSA=makeDict(self.graph)
    bestGraphSA=nx.DiGraph(bestGraphSA)
    bestBIC=self.BIC

    temperature=startingTemp

    for i in range(max_iter):
        mode = randint(1,3)
        if mode==1:
            graphSA,affected_nodes=deleteEdge(bestGraphSA)
        elif mode==2:
            graphSA,affected_nodes=addEdge(bestGraphSA)
        elif mode==3:
            graphSA,affected_nodes=reverseEdge(bestGraphSA)

        graphSA=makeAdjFromNX(graphSA,self.graph_size)
        newBIC=computeBIC(graphSA, self.samples,self.categories,affected_nodes,bestBIC)

        if sum(newBIC) < sum(bestBIC):
            delta=abs(sum(newBIC))-abs(sum(bestBIC))
            acceptanceProb = np.exp(-(delta)/temperature)
            #print delta, acceptanceProb
        else:
            acceptanceProb=0

        if sum(newBIC)>sum(bestBIC) or uniform(0,1)<acceptanceProb:
            bestGraphSA=makeDict(graphSA)
            bestGraphSA=nx.DiGraph(bestGraphSA)
            bestBIC=newBIC[:]

        temperature*=alpha

    self.graph = makeAdjFromNX(bestGraphSA, self.graph_size)
    self.BIC = sum(bestBIC)
    print self.BIC
    print self.graph

def evolutionaryAlgorithm(self, max_iter, size_initial_population, nr_survivors, nr_crossovers, nr_mutations):

#---------------------------- initialization ---------------------------------#
# create an initial population of random graphs

    population={}
    BICs={} 
    BICsums={}

    for i in range(size_initial_population):
        graph=newRandomDAG(self.graph_size, self.p_link, self.max_parents)
        BIC=computeBIC(graph, self.samples, self.categories)
        population.update({i:graph})
        BICs.update({i:BIC})
        BICsums.update({i:sum(BIC)})
    print 'size of initial population:',i+1   
    nr_survivors=i
    indexBest=max(BICsums.iterkeys(), key=(lambda key: BICsums[key]))
    print 'Best BIC from initial population:', sum(BICs[indexBest])

#--------------------------- Evolutionary Algorithm --------------------------#
    generation=0
    for i in range(max_iter):
        generation+=1
        survivors={}
        # select promising individuals
        survivorKeys=[]
        tempBICsums=dict(BICsums)
        for i in range(nr_survivors):
            survivorKeys.append(max(tempBICsums.iterkeys(), key=(lambda key: tempBICsums[key])))
            survivors.update({i:population[survivorKeys[i]]})
            population.update({i:population[survivorKeys[i]]})
            # copy the scores of the survivors into the first slots of the dict
            BICs.update({i:BICs[survivorKeys[i]]})
            BICsums.update({i:sum(BICs[survivorKeys[i]])})
            tempBICsums.pop(survivorKeys[i],None)

        # crossover survivors to create off springs
        for i in range(nr_crossovers):
            graphNumber=i+nr_survivors
            if i==0:
            # crossover by using locally best scores of survivors
                crossedoverGraph,crossedoverBIC,affected_nodes=crossoverBest(survivors,BICs,self.graph_size)
            
                population.update({graphNumber:crossedoverGraph})
                BIC=computeBIC(crossedoverGraph,self.samples, self.categories,affected_nodes,crossedoverBIC)
            # evaluate modified individuals
                BICs.update({graphNumber:BIC})
                BICsums.update({graphNumber:sum(BIC)})
            else:
            # cross over by choosing local sequence randomly from survivors     
                crossedoverGraph,crossedoverBIC,affected_nodes=crossoverRandom(survivors,BICs,self.graph_size)
                population.update({graphNumber:crossedoverGraph})
                BIC=computeBIC(crossedoverGraph,self.samples,self.categories,affected_nodes,crossedoverBIC)
            # evaluate modified individuals
                BICs.update({graphNumber:BIC})
                BICsums.update({graphNumber:sum(BIC)})

        # mutate survivors until population size is reached
        for i in range(nr_mutations):
            survivorToMutate=randint(0,nr_survivors-1)
            graphNumber=i+nr_crossovers+nr_survivors
            graphToMutate=makeDict(survivors[survivorToMutate])
            graphToMutate=nx.DiGraph(graphToMutate)
            mode = randint(1,3)
            if mode==1: #Mutate by deleting edge
                graphMutated,affected_nodes=deleteEdge(graphToMutate)
                graphMutated=makeAdjFromNX(graphMutated,self.graph_size)
                population.update({graphNumber:graphMutated})
                # evaluate modified individuals
                BICs.update({graphNumber:computeBIC(graphMutated,self.samples,self.categories,affected_nodes,BICs[survivorToMutate])})
                BICsums.update({graphNumber:sum(BICs[graphNumber])})
            elif mode==2: #Mutate by adding edge
                graphMutated,affected_nodes=addEdge(graphToMutate)
                graphMutated=makeAdjFromNX(graphMutated,self.graph_size)
                population.update({graphNumber:graphMutated})
                # evaluate modified individuals
                BICs.update({graphNumber:computeBIC(graphMutated,self.samples,self.categories,affected_nodes,BICs[survivorToMutate])})
                BICsums.update({graphNumber:sum(BICs[graphNumber])})
            elif mode==3: #Mutate by reversing edge
                graphMutated,affected_nodes=reverseEdge(graphToMutate)
                graphMutated=makeAdjFromNX(graphMutated,self.graph_size)
                population.update({graphNumber:graphMutated})
                # evaluate modified individuals
                BICs.update({graphNumber:computeBIC(graphMutated,self.samples,self.categories,affected_nodes,BICs[survivorToMutate])})
                BICsums.update({graphNumber:sum(BICs[graphNumber])})

    indexBest=max(BICsums.iterkeys(), key=(lambda key: BICsums[key]))

    self.BIC = sum(BICs[indexBest])
    self.graph = population[indexBest]

    print self.BIC
    print self.graph
