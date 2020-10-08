import random

import numpy as np
import networkx as nx
import matplotlib.pyplot as plt


def get_degree_num_dict(node_degree_dict):
    d_list = [nd[1] for nd in node_degree_dict]
    d_num_dict = {}
    for d in d_list:
        if d not in d_num_dict.keys():
            d_num_dict[d] = 1
        else:
            d_num_dict[d] += 1
    return d_num_dict


def plot_degree_bar(degree_num_dict):
    x, y = degree_num_dict.keys(), degree_num_dict.values()
    y = np.array(list(y))
    print(f'{max(y)} {min(y)}')
    y = y / num_nodes
    x = np.array(list(x))
    print(f'{max(x)} {min(x)}')
    plt.xlim(-1, 50)
    plt.ylim(0, 0.5)
    ax = plt.gca()
    ax.yaxis.set_major_locator(plt.MultipleLocator(0.05))
    ax.xaxis.set_major_locator(plt.MultipleLocator(5))
    plt.bar(x, y)
    # plt.show()


def uniform_edge_remove(remove_num):
    for _ in range(int(remove_num)):
        target_index = np.random.randint(low=0, high=node_num, size=1)[0]
        target_node = node_list[target_index]
        node_edges = list(base_graph.edges(target_node))
        while len(node_edges) == 0:
            target_index = np.random.randint(low=0, high=node_num, size=1)[0]
            target_node = node_list[target_index]
            node_edges = list(base_graph.edges(target_node))
        target_edge_index = np.random.randint(low=0, high=len(node_edges), size=1)[0]
        target_edge = node_edges[target_edge_index]
        base_graph.remove_edge(target_edge[0], target_edge[1])


def preferential_edge_remove():
    removed_indices = np.random.randint(low=0, high=edge_num, size=int(frac_removed * edge_num))
    removed_edges = edge_list[removed_indices].tolist()
    removed_edges = [(k[0], k[1]) for k in removed_edges]
    base_graph.remove_edges_from(removed_edges)


def save_graph(g, save_name):
    nodes, edges = g.nodes, g.edges
    param_dict = {
        'nodes': nodes,
        'edges': edges
    }
    # pickle.dump(param_dict, file=f'network/{save_name}')
    np.save(f'network/{save_name}.npy', param_dict, allow_pickle=True)


def load_graph(name):
    p_dict = np.load(f'network/{name}.npy', allow_pickle=True).item()
    nodes, edges = p_dict['nodes'], p_dict['edges']
    G = nx.Graph()
    G.add_nodes_from(nodes)
    G.add_edges_from(edges)
    return G


if __name__ == "__main__":
    # manual_seed = 0
    frac_removed = 2.6
    # random.seed(manual_seed)
    # np.random.seed(manual_seed)
    num_nodes = 10000
    # base_graph = nx.barabasi_albert_graph(n=num_nodes, m=9, seed=manual_seed)
    base_graph = load_graph(name='baseline')
    close_graph = load_graph(name='close_contact')
    node_list, edge_list = np.array(base_graph.nodes), np.array(base_graph.edges)
    node_num, edge_num = len(node_list), len(edge_list)
    print(f'Initial node:{node_num}\tedge:{edge_num}')
    # preferential_edge_remove()
    # uniform_edge_remove(remove_num=frac_removed * edge_num)
    # print(f'After Removal: edge:{len(base_graph.edges)}')
    nd_dict = base_graph.degree(base_graph.nodes)
    deg_num_dict = get_degree_num_dict(nd_dict)
    plot_degree_bar(deg_num_dict)
    nd_dict = close_graph.degree(close_graph.nodes)
    deg_num_dict = get_degree_num_dict(nd_dict)
    plot_degree_bar(deg_num_dict)
    # save_graph(g=base_graph, save_name='close_contact')
    plt.xlabel('network degree')
    plt.ylabel('% of all nodes')
    plt.legend(['baseline network G', 'close contact network Q'])
    # plt.show()
    plt.savefig('figure/network_degree_distribution.png')
