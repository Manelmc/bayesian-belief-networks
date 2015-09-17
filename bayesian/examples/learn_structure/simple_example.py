import bayesian.bbn as bg
import bayesian.structure_learning as sl
import numpy as np
from scipy import random, linalg



if __name__ == '__main__':

    n_samples = 400
    n_variables = 5

    var_names = []
    for n in range(n_variables):
        var_names.append('var_%d'%n)

    # Generate random mean vector and covariance matrix
    mu = random.sample(n_variables)
    A = random.rand(n_variables, n_variables)
    R = np.dot(A, A.transpose())

    # Generate random samples of correlated variables
    data = np.random.multivariate_normal(mu, R, size=n_samples)


    bin_size = 3
    max_parents = 2
    p_link = 0.7

    max_iter = 100

    size_initial_population = 20
    nr_survivors = 5
    nr_crossovers = 8
    nr_mutations = 7

    initial_graph = sl.Graph(var_names, p_link, max_parents, bin_size)

    learn_method = sl.GeneticAlgorithm(
                        size_initial_population,
                        nr_survivors,
                        nr_crossovers,
                        nr_mutations)

    model = sl.StructureLearner(initial_graph, data, learn_method, max_iter)
    model.learn()

    graphs, X_list, bins_list = model.separate_components()

    for n, g in enumerate(graphs):
        print g.adj_matrix
        bd = bg.build_bbn_from_data(g, X_list[n])
        bd.export('graph_%d.gv'%n)
        bd.q()
