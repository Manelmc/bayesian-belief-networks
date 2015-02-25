import numpy as np
from random import *
import sys
import math
import pandas as pd
import copy as cp

class Graph:
    """docstring for ClassName"""
    def __init__(self, names, p_link, max_parents, bin_size):
        self.names = names
        self.graph_size = len(names)
        self.p_link = p_link
        self.max_parents = max_parents

        self.adj_matrix = np.zeros((self.graph_size, self.graph_size))
        self.new_rnd()
        self.categories = np.repeat(bin_size, self.graph_size)

        self.nodes_w_changed_BIC = []

    def is_DAG(self):
        return strongly_connected_components(self.adj_matrix) == []

    def BIC(self, X):
        return compute_BIC(self, X)

    def new_rnd(self):
        for n in range(self.graph_size): 
            for m in range(self.graph_size):
                linker = uniform(0, 1)
                if (self.p_link > linker and 
                    n != m and 
                    self.adj_matrix[n, m] == 0 and
                    np.sum(self.adj_matrix[n, :]) < self.max_parents):
                    self.adj_matrix[m, n] = 1

                # To increase the number of possible graphs
                # switch edge direction at random
                if (self.adj_matrix[n, m] == 1 and
                    np.sum(self.adj_matrix[m, :]) < self.max_parents):
                    switcher = randint(0, 1)
                    if switcher == 1:
                        self.adj_matrix[n, m] = 0
                        self.adj_matrix[m, n] = 1
        self.make_acyclic() 

    def make_acyclic(self):
        cycles = strongly_connected_components(self.adj_matrix)
        self.last_changed_nodes = []

        while len(cycles) != 0:
            # delete one parent of one node in each cycle
            # per iteration
            for c in cycles:
                node_in_cycle = randint(0, len(c) - 1)
                parents = np.where(self.adj_matrix[node_in_cycle,:] == 1)[0]
                # Only consider deleting parents if they are in the cycle
                parents = list(set(c)&set(parents))
                if len(parents) == 0: continue
                # A node can be in multiple cycles. 
                # Add to list of last changed nodes
                # only if not in it yet.
                if node_in_cycle not in self.last_changed_nodes:
                    self.last_changed_nodes.append(node_in_cycle)

                parent_to_delete = randint(0, len(parents) -1 )

                self.adj_matrix[node_in_cycle, parents[parent_to_delete]] = 0

            cycles = strongly_connected_components(self.adj_matrix)

    def add_edge(self):
        pass
    def delete_edge(self):
        pass
    def reverse_edge(self):
        pass
    def count_combinations(self):
        pass

def strongly_connected_components(graph):
    """
    Tarjan's Algorithm (named for its discoverer, Robert Tarjan)
    is a graph theory algorithm for finding the 
    strongly connected components of a graph.
    """

    index_counter = [0]
    stack = []
    lowlinks = {}
    index = {}
    cycles = []
    
    def strongconnect(node):
        index[node] = index_counter[0]
        lowlinks[node] = index_counter[0]
        index_counter[0] += 1
        stack.append(node)
        
        successors = np.where(graph[:,node] == 1)[0]

        for successor in successors:
            if successor not in lowlinks:
                strongconnect(successor)
                lowlinks[node] = min(lowlinks[node],lowlinks[successor])
            elif successor in stack:
                lowlinks[node] = min(lowlinks[node],index[successor])
        
        if lowlinks[node] == index[node]:
            connected_component = []
            
            while True:
                successor = stack.pop()
                connected_component.append(successor)
                if successor == node: break
            component = tuple(connected_component)

            cycles.append(component)
    
    for node in range(graph.shape[0]):
        if node not in lowlinks:
            strongconnect(node)

    cycles[:] = [x for x in cycles if len(x) > 1]
    
    return cycles

def count_combinations(self, X):

    combinationsList={}
    numberOfCombinationsList={}
        
    for node in self.nodes_w_changed_BIC:
        Parents=np.where(self.adj_matrix[node,:]==1)[0]
        #get the number of combinations by computing the product of categories in the parents
        numberOfCombinations = np.prod(self.categories[Parents])
        numberOfCombinationsList.update({node:numberOfCombinations})
        numberOfCombinationsIncludingNode=numberOfCombinations*self.categories[node]
        combinationCountsNode=np.zeros(numberOfCombinationsIncludingNode)

        for samplePieceN in X:
            #get the combinationNr for one node in sample n
            combinationNr=0
            for index,parent in enumerate(Parents):
                if index==0:
                    combinationNr=samplePieceN[parent]*self.categories[node]
                    catParentOld=self.categories[parent]*self.categories[node]
                else:
                    combinationNr+=samplePieceN[parent]*catParentOld
                    catParentOld*=self.categories[parent]

            combinationNr=combinationNr+samplePieceN[node]
            combinationCountsNode[combinationNr]+=1

        combinationsList.update({node:combinationCountsNode})

    return combinationsList,numberOfCombinationsList


def compute_BIC(self, X, BIC_old=[]):
# compute the Bayesian Information Criterion for a graph, given data
    if self.nodes_w_changed_BIC == []:
        self.nodes_w_changed_BIC = np.arange(0, self.graph_size)
        BIC_new = np.zeros(len(self.nodes_w_changed_BIC))
    else:
        BIC_new = np.copy(BIC_old)

    Nijk, qi = count_combinations(self, X)
    Nij = {}
    for node in Nijk:
        nodeNij=[]
        index=0
        for i in range(qi[node]):
            nodeNij.append(sum(Nijk[node][index:index+self.categories[node]]))
            index += self.categories[node]
        Nij.update({node:nodeNij})

    LL={}
    for node in self.nodes_w_changed_BIC:
        index=0
        LLNode=0
        for combinationCount in range(len(Nij[node])):
            for categoryCount in range(self.categories[node]):
                if Nij[node][combinationCount]==0:
                    Nij[node][combinationCount]=1
                if Nijk[node][categoryCount+index]==0:
                    Nijk[node][categoryCount+index]=1
                LLNode+=Nijk[node][categoryCount+index]*math.log(Nijk[node][categoryCount+index]/Nij[node][combinationCount])
            index+=self.categories[node]
        LL.update({node:LLNode})

    self.categories=np.array(self.categories)
    qi=np.array(qi.items())[:,1]    

    complexityPenalty=0.5*math.log(len(X))*abs((self.categories[self.nodes_w_changed_BIC]-1)*qi)
    for index, node in enumerate(self.nodes_w_changed_BIC):
        BIC_new[node]=LL[node]-complexityPenalty[index]

    self.nodes_w_changed_BIC = []

    return BIC_new

def make_samples(X, bin_size):
    '''Create samples from data by discretizing over specified bins. 
    Return the samples in form of a dict and the bins for each node'''
    samples = np.empty(X.shape)
    bins = []
    for i, node in enumerate(X.T):
        samples[:,i] = np.digitize(node, np.histogram(node, bins=bin_size)[1],right=True)-1
        bins.append(np.histogram(node, bins=bin_size)[1])
    return samples, bins


class StructureLearner:
    """docstring for StructureLearner"""
    def __init__(self, graph, X, learn_method, max_iter):

        self.X = X
        self.graph = graph
        self.learn_method = learn_method
        self.max_iter = max_iter

    def learn(self):
        iteration = 0
        converged = False
        self.learn_method._setup(self.graph, self.X)

        while iteration < self.max_iter and not converged:
            iteration += 1
            self.learn_method.step()
        print self.learn_method.best_graph.adj_matrix
        print self.learn_method.best_score


class LearningMethod:
    def _setup(self):
        raise NotImplementedError()

    def step(self):
        raise NotImplementedError()

class HillClimbing(LearningMethod):
    """ Add, delete or reverse an edge.
        If this improves the BIC save the graph. 
        Stop after maximum number of iterations."""
    def __init__(self):
        pass

    def _setup(self, graph, X):
        self.X = X
        self.best_graph = graph
        self.new_graph = cp.deepcopy(self.best_graph)
        self.best_score = graph.BIC(X)

    def step(self):
        self.new_graph.adj_matrix = cp.copy(self.best_graph.adj_matrix)
        mode = randint(1, 3)

        if mode == 1:
            self.new_graph.delete_edge()
        elif mode == 2:
            self.new_graph.add_edge()
        elif mode == 3:
            self.new_graph.reverse_edge()

        new_score = self.new_graph.BIC(self.X)

        if sum(new_score) > sum(self.best_score):
            self.best_score = new_score
            self.best_graph.adj_matrix = cp.copy(self.new_graph.adj_matrix)

class SimulatedAnnealing(LearningMethod):
    """ Add, delete or reverse an edge.
        Accept the change if the BIC is better,
        or if the BIC is worse with a certain
        changing acceptance probality"""

    def __init__(self, starting_temperature=1, alpha=0.95):
        self.temperature = starting_temperature
        self.alpha = alpha

    def _setup(self, graph, X):
        self.X = X
        self.best_graph = graph
        self.new_graph = cp.deepcopy(self.best_graph)
        self.best_score = graph.BIC(X)

    def step(self):
        self.new_graph.adj_matrix = cp.copy(self.best_graph.adj_matrix)
        mode = randint(1, 3)
        if mode == 1:
            self.new_graph.delete_edge()
        elif mode == 2:
            self.new_graph.add_edge()
        elif mode == 3:
            self.new_graph.reverse_edge()

        new_score = self.new_graph.BIC(self.X)

        if sum(new_score) < sum(self.best_score):
            delta = abs(sum(new_score)) - abs(sum(self.best_score))
            p_accept = np.exp(-(delta) / self.temperature)

            if uniform(0, 1) < p_accept:
                self.best_score = new_score
                self.best_graph.adj_matrix = cp.copy(self.new_graph.adj_matrix)  
        else:
            self.best_score = new_score
            self.best_graph.adj_matrix = cp.copy(self.new_graph.adj_matrix)

        self.temperature *= self.alpha        

class GeneticAlgorithm(LearningMethod):
    """docstring for GeneticAlgorithm"""
    def __init__(self, size_initial_population=20, nr_survivors=5, nr_crossovers=7, nr_mutations=8):
        self.size_population = size_initial_population
        self.nr_survivors = nr_survivors
        self.nr_crossovers = nr_crossovers
        self.nr_mutations = nr_mutations
        self.population = {}
        self.pop_BICs = {}
        self.pop_BICs_per_node = {} 

    def _setup(self):
        for i in range(self.size_population):
            self.population.update({i:Graph.new_rnd(self.graph_size)})
            self.pop_BICs_per_node.update({i:self.population[i].BIC()})
            self.pop_BICs.update({i:sum(self.pop_BICs_per_node[i])})
        index_best = max(self.pop_BICs.iterkeys(), key=(lambda key: self.pop_BICs[key]))
        self.best_score = self.pop_BICs[index_best]
        self.best_graph = self.population[index_best]

    def step(self):
        # the x best graphs survive unchanged
        # cross over best nodes amongst survivors and random nodes from survivors until number of cross overs is reached
        # mutate survivors until number of mutations is reached
        pass

if __name__ == '__main__':
    data = pd.read_excel('../../../data_bnn.xlsx', 'data_1', index_col=None, na_values=['NA'])
    data.dropna(inplace=True)

    names = [str(name) for name in data.columns.values]

    bin_size = 3
    max_parents = 2
    p_link = 0.4
    max_iter = 10

    graph = Graph(names, p_link, max_parents, bin_size)
    X,_ = make_samples(data.as_matrix(), bin_size)
    # learn_method = HillClimbing()
    learn_method = SimulatedAnnealing()
    SL = StructureLearner(graph, X, learn_method, max_iter)
    SL.learn()

        