import numpy as np
import random as rnd
import math
import copy as cp
from scipy import sparse as sp
import sys


class Graph:

    def __init__(self, names, p_link, max_parents, bin_size):
        self.names = names
        self.size = len(names)
        self.bin_size = bin_size
        self.p_link = p_link
        self.max_parents = max_parents

        self.adj_matrix = np.zeros((self.size, self.size))
        self.new_rnd()
        self.categories = np.repeat(bin_size, self.size)

        self.nodes_w_changed_BIC = []

    def remove_nodes(self, node_indices):
        all_nodes = np.arange(self.size)
        ind_keep = [i for j, i in enumerate(all_nodes) if j not in node_indices]

        self.names = list(self.names[i] for i in ind_keep)
        self.adj_matrix = self.adj_matrix[np.ix_(ind_keep, ind_keep)]
        self.categories = self.categories[np.ix_(ind_keep)]
        self.size = len(self.names)

    def is_DAG(self):
        return strongly_connected_components(self.adj_matrix) == []

    def BIC(self, X, BIC_old=[]):
        return compute_BIC(self, X, BIC_old)

    def new_rnd(self):
        for n in range(self.size): 
            for m in range(self.size):
                linker = rnd.uniform(0, 1)
                if (self.p_link > linker and 
                    n != m and 
                    self.adj_matrix[n, m] == 0 and
                    np.sum(self.adj_matrix[n, :]) < self.max_parents):
                    self.adj_matrix[m, n] = 1

                # To increase the number of possible graphs
                # switch edge direction at random
                if (self.adj_matrix[n, m] == 1 and
                    np.sum(self.adj_matrix[m, :]) < self.max_parents):
                    switcher = rnd.randint(0, 1)
                    if switcher == 1:
                        self.adj_matrix[n, m] = 0
                        self.adj_matrix[m, n] = 1
        try: 
            self.make_acyclic()
        except:
            self.new_rnd() 

    def make_acyclic(self):
        cycles = strongly_connected_components(self.adj_matrix)
        self.last_changed_nodes = []
        counter = 0
        while len(cycles) > 0:
            counter += 1
            assert counter < 20
            # delete one parent of one node in each cycle
            # per iteration
            for c in cycles:
                node_in_cycle = rnd.randint(0, len(c) - 1)
                parents = np.where(self.adj_matrix[node_in_cycle,:] == 1)[0]
                # Only consider deleting parents if they are in the cycle
                parents = list(set(c)&set(parents))
                if len(parents) == 0: continue
                # A node can be in multiple cycles. 
                # Add to list of last changed nodes
                # only if not in it yet.
                if node_in_cycle not in self.last_changed_nodes:
                    self.last_changed_nodes.append(node_in_cycle)

                parent_to_delete = rnd.randint(0, len(parents) -1 )

                self.adj_matrix[node_in_cycle, parents[parent_to_delete]] = 0

            cycles = strongly_connected_components(self.adj_matrix)

    def add_edge(self):
        non_edges = np.where(self.adj_matrix == 0)
        non_edges = zip(non_edges[0], non_edges[1])
        non_edges = [ne for ne in non_edges if ne[0]!=ne[1]]

        nr_non_edges = len(non_edges)

        edge_to_add = rnd.randint(0, nr_non_edges - 1)

        self.adj_matrix[non_edges[edge_to_add]] = 1

        counter = 0

        while not self.is_DAG() and counter < 20:
            counter += 1
            self.adj_matrix[non_edges[edge_to_add]] = 0

            edge_to_add = rnd.randint(0, nr_non_edges - 1)

            self.adj_matrix[non_edges[edge_to_add]] = 1  

        if self.is_DAG() and counter < 20:
            self.nodes_w_changed_BIC.append(non_edges[edge_to_add][0])
        else:
            self.adj_matrix[non_edges[edge_to_add]] = 0


    def delete_edge(self):
        edges = np.where(self.adj_matrix == 1)
        edges = zip(edges[0], edges[1])

        nr_edges = len(edges)
        assert nr_edges > 0, 'No edges in graph.'
        edge_to_delete = rnd.randint(0, nr_edges - 1)

        self.adj_matrix[edges[edge_to_delete]] = 0
        self.nodes_w_changed_BIC.append(edges[edge_to_delete][1])

    def reverse_edge(self):
        edges = np.where(self.adj_matrix == 1)
        edges = zip(edges[0], edges[1])

        nr_edges = len(edges)
        assert nr_edges is not 0, 'No edges in graph.'

        edge_to_reverse = rnd.randint(0, nr_edges - 1)

        self.adj_matrix[edges[edge_to_reverse]] = 0
        self.adj_matrix[edges[edge_to_reverse][::-1]] = 1

        counter = 0
        while not self.is_DAG() and counter < 20:
            counter += 1
            self.adj_matrix[edges[edge_to_reverse]] = 1
            self.adj_matrix[edges[edge_to_reverse][::-1]] = 0  

            edge_to_reverse = rnd.randint(0, nr_edges - 1)

            self.adj_matrix[edges[edge_to_reverse]] = 0
            self.adj_matrix[edges[edge_to_reverse][::-1]] = 1  

        if self.is_DAG() and counter < 20:
            self.nodes_w_changed_BIC.append(list(edges[edge_to_reverse])[0])
        else:
            self.adj_matrix[edges[edge_to_reverse]] = 1
            self.adj_matrix[edges[edge_to_reverse][::-1]] = 0 
        

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


def count_combinations(graph, X):
    '''
    - N_ijk: The number of instances in the data in which X_i takes the k-th value and the parents of X_i take configuration j.
    - q_i: Number of domain-combinations of parents of node i

    '''

    N_ijk = {}
    q_i = {}
        
    for node in graph.nodes_w_changed_BIC:
        parents = np.where(graph.adj_matrix[node, :] == 1)[0]

        q_i.update({node:np.prod(graph.categories[parents])})

        combination_counts = np.zeros(q_i[node] * graph.categories[node])

        for sample in X:
            combination_nr = 0
            for p, parent in enumerate(parents):
                if p == 0:
                    combination_nr = sample[parent] * graph.categories[node]
                    cat_prev_p = graph.categories[parent] * graph.categories[node]
                else:
                    combination_nr += sample[parent] * cat_prev_p
                    cat_prev_p *= graph.categories[parent]

            combination_nr = combination_nr + sample[node]
            combination_counts[combination_nr] += 1

        N_ijk.update({node:combination_counts})

    return N_ijk, q_i


def compute_BIC(graph, X, BIC_old = []):
    '''compute the Bayesian Information Criterion for a graph, 
    given data. Old BIC can be passed to speed up computation.'''

    if graph.nodes_w_changed_BIC == []:
        graph.nodes_w_changed_BIC = np.arange(0, graph.size)
        BIC_new = np.zeros(len(graph.nodes_w_changed_BIC))
    else:
        BIC_new = np.copy(BIC_old)

    N_ijk, q_i = count_combinations(graph, X)

    N_ij = {}
    for node in N_ijk:
        nodeN_ij = []
        index = 0

        for i in range(q_i[node]):
            nodeN_ij.append(sum(N_ijk[node][index:index + graph.categories[node]]))
            index += graph.categories[node]

        N_ij.update({node:nodeN_ij})

    LL = {}
    for node in graph.nodes_w_changed_BIC:
        index = 0
        LL_node = 0
        for combinationCount in range(len(N_ij[node])):
            for categoryCount in range(graph.categories[node]):
                if N_ij[node][combinationCount] == 0:
                    N_ij[node][combinationCount] = 1

                if N_ijk[node][categoryCount + index] == 0:
                    N_ijk[node][categoryCount + index] = 1

                LL_node += N_ijk[node][categoryCount + index] \
                        * math.log(N_ijk[node][categoryCount + index]\
                        / N_ij[node][combinationCount])

            index += graph.categories[node]
        LL.update({node:LL_node})

    graph.categories = np.array(graph.categories)
    q_i = np.array(q_i.items())[:,1]    

    complexityPenalty = 0.5 * math.log(len(X)) * abs((graph.categories[graph.nodes_w_changed_BIC]-1) * q_i)

    for index, node in enumerate(graph.nodes_w_changed_BIC):
        BIC_new[node] = LL[node] - complexityPenalty[index]

    graph.nodes_w_changed_BIC = []

    return BIC_new

def make_samples(X, bin_size):
    '''Create samples from data by discretizing over specified bins. 
    Return the samples in form of a dict and the bins for each node'''
    samples = np.empty(X.shape)
    bins = []
    for i, node in enumerate(X.T):
        samples[:,i] = np.digitize(node, np.histogram(node, bins=bin_size)[1], right=True) - 1
        bins.append(np.histogram(node, bins=bin_size)[1])
    return samples, bins

def find_connected_components(graph):
    adj_list = []

    for i in range(graph.size):
        for j in range(graph.size):
            if graph.adj_matrix[i, j] == 1:
                adj_list.append([i, j])
    i_indices, j_indices = zip(*adj_list)

    sparse_matrix = sp.csr_matrix(
        (np.ones(len(j_indices)), (i_indices, j_indices)),
        shape=(graph.size, graph.size))

    return sp.csgraph.connected_components(sparse_matrix)


class StructureLearner:

    def __init__(self, graph, data, learn_method, max_iter):

        self.data = data
        self.graph = graph
        self.learn_method = learn_method
        self.max_iter = max_iter
        self.score = 0

        self.X, self.bins = make_samples(data, self.graph.bin_size)

    def learn(self):
        iteration = 0
        converged = False
        patience = 10
        self.learn_method._setup(self.graph, self.X)

        while iteration < self.max_iter and not converged:
            iteration += 1
            self.learn_method.step()
            if sum(self.learn_method.delta) == 0 and iteration > patience: 
                converged = True
                # print 'Converged after %d iterations.'%iteration
                # sys.stdout.flush()


        self.graph = self.learn_method.best_graph
        self.score = sum(self.learn_method.best_score)

        # print('BIC: %f'%self.score)

    def separate_components(self):
        """Sometimes, the learned structure will consist of multiple graphs.
        In this case, we want to separate the components into separate graphs.
        We also remove components with single nodes, as they do not hold any
        valuable information"""
        components = find_connected_components(self.graph)
        graphs = []
        X_list = []
        bins_list = []
        all_nodes = np.arange(self.graph.size)

        values, indices = np.unique(components[1], return_index=True)
        for val in values:
            len_component = list(components[1]).count(val)
            if len_component > 1:
                i_new_graph = [k for k, v
                               in enumerate(components[1])
                               if v == val]
                i_rem = [k for j, k
                         in enumerate(all_nodes)
                         if j not in i_new_graph]

                new_graph = cp.deepcopy(self.graph)
                new_graph.remove_nodes(i_rem)
                graphs.append(new_graph)

                X_list.append(self.X.T[np.ix_(i_new_graph)].T)
                bins_list.append(list(self.bins[k] for k in i_new_graph))

        return graphs, X_list, bins_list

class HillClimber():
    """Add, delete or reverse an edge.
    If this improves the BIC save the graph. 
    Stop after maximum number of iterations."""

    def __init__(self):
        self.n_delete_fails = 0
        self.delta = [0]

    def _setup(self, graph, X):
        self.X = X
        self.best_graph = graph
        self.new_graph = cp.deepcopy(self.best_graph)
        self.best_score = graph.BIC(X)

    def step(self):
        self.new_graph.adj_matrix = cp.copy(self.best_graph.adj_matrix)
        mode = rnd.randint(1, 3)

        if mode == 1:
            self.new_graph.add_edge()
        elif mode == 2:
            try:
                self.new_graph.delete_edge()
            except:
                self.n_delete_fails += 1
                self.new_graph.add_edge()
        elif mode == 3:
            self.new_graph.reverse_edge()

        # If there is a bi-directional connection, delete one at random

        self.new_graph.clear_bi_dir

        new_score = self.new_graph.BIC(self.X, self.best_score)

        if sum(new_score) > sum(self.best_score):
            self.delta = abs(self.best_score - new_score)
            self.best_score = new_score
            self.best_graph.adj_matrix = cp.copy(self.new_graph.adj_matrix)

class SimulatedAnnealing():
    """ Add, delete or reverse an edge.
        Accept the change if the BIC is better,
        or if the BIC is worse with a decreasing
        changing acceptance probability"""

    def __init__(self, starting_temperature=1, alpha=0.95):
        self.temperature = starting_temperature
        self.alpha = alpha
        self.n_delete_fails = 0
        self.delta = [0]

    def _setup(self, graph, X):
        self.X = X
        self.best_graph = graph
        self.new_graph = cp.deepcopy(self.best_graph)
        self.best_score = graph.BIC(X)

    def step(self):
        assert self.n_delete_fails < 20, 'There is no causal relation between the variables.'
        self.new_graph.adj_matrix = cp.copy(self.best_graph.adj_matrix)
        mode = rnd.randint(1, 3)
        if mode == 1:
            self.new_graph.add_edge()
        elif mode == 2:
            try:
                self.new_graph.delete_edge()
            except:
                self.n_delete_fails += 1
                self.new_graph.add_edge()
        elif mode == 3:
            self.new_graph.reverse_edge()

        new_score = self.new_graph.BIC(self.X, self.best_score)

        if sum(new_score) < sum(self.best_score):
            delta = abs(sum(new_score)) - abs(sum(self.best_score))
            p_accept = np.exp(-(delta) / self.temperature)

            if rnd.uniform(0, 1) < p_accept:
                self.delta = abs(self.best_score - new_score)
                self.best_score = new_score
                self.best_graph.adj_matrix = cp.copy(self.new_graph.adj_matrix)  
        else:
            self.delta = abs(self.best_score - new_score)
            self.best_score = new_score
            self.best_graph.adj_matrix = cp.copy(self.new_graph.adj_matrix)

        self.temperature *= self.alpha        

class GeneticAlgorithm():
    """Create an initial population of random graphs.
    In each generation, let the best graphs remain 
    in the population (=survivors). Mutate and crossover 
    individual graphs from the survivors."""


    def __init__(self, size_initial_population=40, nr_survivors=10,
                 nr_crossovers=10, nr_mutations=20):
        self.size_population = size_initial_population
        self.nr_survivors = nr_survivors
        self.nr_crossovers = nr_crossovers
        self.mutation_start = nr_survivors + nr_crossovers
        self.nr_mutations = nr_mutations
        self.n_delete_fails = 0
        self.delta = [0]

    def _setup(self, graph, X):
        self.X = X
        self.graph_size = graph.size
        self.population = {}
        self.population_BICs = np.empty((self.size_population, graph.size))

        for i in range(self.size_population):
            new_graph = cp.deepcopy(graph)
            new_graph.new_rnd()
            self.population.update({i:new_graph})
            self.population_BICs[i, :] = self.population[i].BIC(X)

        index_best = np.argmax(np.sum(self.population_BICs, 1))
        self.best_score = self.population_BICs[index_best]
        self.best_graph = self.population[index_best]

    def step(self):
        assert self.n_delete_fails < 20*self.size_population,\
            'There is no causal relation between the variables.'
        self.select_fittest()

        for c in range(self.nr_survivors, self.nr_survivors + self.nr_crossovers):
            self.crossover(c)

        for m in range(self.mutation_start, self.mutation_start + self.nr_mutations):
            self.mutate(m)

        index_best = np.argmax(np.sum(self.population_BICs, 1))
        new_score = self.population_BICs[index_best]
        self.delta = abs(self.best_score - new_score)
        self.best_score = new_score
        self.best_graph = self.population[index_best]

    def select_fittest(self):
        """Selecting the fittest individuals means keeping the graphs
        with the highest BIC. Those graphs are called survivors."""
        fittest_index = np.argpartition(
            np.sum(self.population_BICs, 1), -self.nr_survivors)\
                [-self.nr_survivors:]
        
        self.population.update({i: self.population[s] for i, s in enumerate(fittest_index)})

        for i, s in enumerate(fittest_index):
            self.population_BICs[i, :] = self.population_BICs[s, :]

    def mutate(self, m):
        """We mutate randomly chosen survivors by adding, deleting or
        reversing a randomly chosen edge"""
        parent_index = rnd.randint(0, self.nr_survivors - 1)

        offspring = cp.deepcopy(self.population[parent_index])

        mode = rnd.randint(1, 3)
        if mode == 1:
            offspring.add_edge()
        elif mode == 2:
            try:
                offspring.delete_edge()
            except:
                self.n_delete_fails += 1
                offspring.add_edge()
        elif mode == 3:
            try:
                offspring.reverse_edge()
            except:                
                offspring.add_edge()        

        self.population.update({m: offspring})

        self.population_BICs[m, :] = offspring.BIC(self.X,
            self.population_BICs[parent_index, :])

    def crossover(self, c):
        """We crossover by choosing the parents of N nodes and merging
        in random order."""
        graph_selection = rnd.sample(xrange(self.size_population), self.graph_size)
        for n, graph in enumerate(graph_selection):
            self.population[c].adj_matrix[n, :] = \
                self.population[graph].adj_matrix[n, :]

            self.population_BICs[c, n] = self.population_BICs[c, n]

        # After crossing over, make sure the resulting graph is acyclic.
        # Should that not be the case, try to make the graph acyclic.
        # If that fails create a new random graph.
        if not self.population[c].is_DAG():
            try:
                self.population[c].make_acyclic()
            except:
                self.population[c].new_rnd()

            self.population_BICs[c, :] = self.population[c].BIC(self.X)