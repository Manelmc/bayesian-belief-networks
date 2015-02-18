import numpy as np
import math
from random_dag import *
import networkx as nx

def getCategories(samples,maxCat):
    categories=[]
    for node in samples[0]:
        counts=np.zeros(maxCat)
        for n in samples:
            counts[n[node]]+=1

        #delete the categories without a count
        c=1
        k=0
        while c!=0 and k<maxCat:
            if counts[k]==0:
                c=0
            else:
                k+=1
        if k<4:
            counts=counts[0:k]
        #save number of categories for each node
        categories.append(len(counts))

    return categories

def getCombinations(graph, samples, categories, effectedNodes):

    combinationsList={}
    numberOfCombinationsList={}
        
    for node in effectedNodes:
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


def computeBIC(estGraph, samples, categories, effectedNodes=[], BICold=[]):
# compute the Bayesian Information Criterion for a graph, given data
    if effectedNodes==[]:
        effectedNodes=np.arange(0,len(estGraph))
        BICnew=np.zeros(len(effectedNodes))
    else:
        BICnew=np.copy(BICold)


    Nijk,qi=getCombinations(estGraph,samples,categories,effectedNodes)
    Nij={}
    for node in Nijk:
        nodeNij=[]
        index=0
        for i in range(qi[node]):
            nodeNij.append(sum(Nijk[node][index:index+categories[node]]))
            index+=categories[node]
        Nij.update({node:nodeNij})

    LL={}
    for node in effectedNodes:
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

    complexityPenalty=0.5*math.log(len(samples))*abs((categories[effectedNodes]-1)*qi)
    for index,node in enumerate(effectedNodes):
        BICnew[node]=LL[node]-complexityPenalty[index]

    return BICnew

def hill_climber(self, max_iter):
    bestGraphHC=makeDict(self.graph)
    bestGraphHC = nx.DiGraph(bestGraphHC)
    bestBIC=self.BIC

    for i in range(max_iter):
        if nx.number_of_edges(bestGraphHC)>2:
            mode = randint(1,3)
        else:
            mode = randint(2,3)
        
        if mode==1:
            graphHC,effectedNodes=deleteEdge(bestGraphHC)
        elif mode==2:
            graphHC,effectedNodes=addEdge(bestGraphHC)
        elif mode==3:
            graphHC,effectedNodes=reverseEdge(bestGraphHC)

        graphHC=makeAdjFromNX(graphHC,self.graph_size)
        newBIC=computeBIC(graphHC, self.samples, self.categories, effectedNodes, bestBIC)

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
        if nx.number_of_edges(bestGraphSA)>2:
            mode = randint(1,3)
        else:
            mode = randint(2,3)

        if mode==1:
            graphSA,effectedNodes=deleteEdge(bestGraphSA)
        elif mode==2:
            graphSA,effectedNodes=addEdge(bestGraphSA)
        elif mode==3:
            graphSA,effectedNodes=reverseEdge(bestGraphSA)

        graphSA=makeAdjFromNX(graphSA,self.graph_size)
        newBIC=computeBIC(graphSA, self.samples,self.categories,effectedNodes,bestBIC)

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

def evolutionaryAlgorithm(self, max_iter, size_initial_population, nrSurvivors, nrCrossovers, nrMutations):

#---------------------------- initialization ---------------------------------#
# create an initial population of random graphs

    pLink=0.3
    maxParents=3

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
    nrSurvivors=i
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
        for i in range(nrSurvivors):
            survivorKeys.append(max(tempBICsums.iterkeys(), key=(lambda key: tempBICsums[key])))
            survivors.update({i:population[survivorKeys[i]]})
            population.update({i:population[survivorKeys[i]]})
            # copy the scores of the survivors into the first slots of the dict
            BICs.update({i:BICs[survivorKeys[i]]})
            BICsums.update({i:sum(BICs[survivorKeys[i]])})
            tempBICsums.pop(survivorKeys[i],None)

        # crossover survivors to create off springs
        for i in range(nrCrossovers):
            graphNumber=i+nrSurvivors
            if i==0:
            # crossover by using locally best scores of survivors
                crossedoverGraph,crossedoverBIC,effectedNodes=crossoverBest(survivors,BICs,self.graph_size)
            
                population.update({graphNumber:crossedoverGraph})
                BIC=computeBIC(crossedoverGraph,self.samples, self.categories,effectedNodes,crossedoverBIC)
            # evaluate modified individuals
                BICs.update({graphNumber:BIC})
                BICsums.update({graphNumber:sum(BIC)})
            else:
            # cross over by choosing local sequence randomly from survivors     
                crossedoverGraph,crossedoverBIC,effectedNodes=crossoverRandom(survivors,BICs,self.graph_size)
                population.update({graphNumber:crossedoverGraph})
                BIC=computeBIC(crossedoverGraph,self.samples,self.categories,effectedNodes,crossedoverBIC)
            # evaluate modified individuals
                BICs.update({graphNumber:BIC})
                BICsums.update({graphNumber:sum(BIC)})

        # mutate survivors until population size is reached
        for i in range(nrMutations):
            survivorToMutate=randint(0,nrSurvivors-1)
            graphNumber=i+nrCrossovers+nrSurvivors
            mode=randint(1,3)
            if mode==1: #Mutate by deleting edge
                graphToMutate=makeDict(survivors[survivorToMutate])
                graphToMutate=nx.DiGraph(graphToMutate)
                graphMutated,effectedNodes=deleteEdge(graphToMutate)
                graphMutated=makeAdjFromNX(graphMutated,self.graph_size)
                population.update({graphNumber:graphMutated})
                # evaluate modified individuals
                BICs.update({graphNumber:computeBIC(graphMutated,self.samples,self.categories,effectedNodes,BICs[survivorToMutate])})
                BICsums.update({graphNumber:sum(BICs[graphNumber])})
            elif mode==2: #Mutate by adding edge
                graphToMutate=makeDict(survivors[survivorToMutate])
                graphToMutate=nx.DiGraph(graphToMutate)
                graphMutated,effectedNodes=addEdge(graphToMutate)
                graphMutated=makeAdjFromNX(graphMutated,self.graph_size)
                population.update({graphNumber:graphMutated})
                # evaluate modified individuals
                BICs.update({graphNumber:computeBIC(graphMutated,self.samples,self.categories,effectedNodes,BICs[survivorToMutate])})
                BICsums.update({graphNumber:sum(BICs[graphNumber])})
            elif mode==3: #Mutate by reversing edge
                graphToMutate=makeDict(survivors[survivorToMutate])
                graphToMutate=nx.DiGraph(graphToMutate)
                graphMutated,effectedNodes=reverseEdge(graphToMutate)
                graphMutated=makeAdjFromNX(graphMutated,self.graph_size)
                population.update({graphNumber:graphMutated})
                # evaluate modified individuals
                BICs.update({graphNumber:computeBIC(graphMutated,self.samples,self.categories,effectedNodes,BICs[survivorToMutate])})
                BICsums.update({graphNumber:sum(BICs[graphNumber])})

    indexBest=max(BICsums.iterkeys(), key=(lambda key: BICsums[key]))

    self.BIC = sum(BICs[indexBest])
    self.graph = population[indexBest]

    print self.BIC
    print self.graph