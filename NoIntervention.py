import argparse
import math
import random
from network import load_graph
import networkx as nx
import matplotlib.pyplot as plt
import numpy as np
import pickle


def get_inf_n_ratio(node_i):
    neighbors = np.array([i for i in base_graph.neighbors(node_i)])
    n_all_adj = len(neighbors)
    if not n_all_adj:
        return .0
    neighbors_states = state_list[neighbors]
    I_neighbors = neighbors[np.where(neighbors_states == 2)]
    n_I_adj = len(I_neighbors)
    return n_I_adj / n_all_adj


def trans_S2E(A, S_node, I_num):
    ip_global = p * beta * I_num / A
    inf_n_ratio_list = [get_inf_n_ratio(node) for node in S_node]
    ip_local = (1 - p) * beta
    return ip_global, ip_local, inf_n_ratio_list


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--base_graph', default=r'baseline', help='name of base graph model')
    parser.add_argument('--log_name', default=r'NoIntervention', help='name of log of current simulation')
    parser.add_argument('--manual_seed', type=int, help='manual random seed', default=0)
    parser.add_argument('--num_node', type=int, default=10000, help='number of nodes in both graphs')
    parser.add_argument('--num_init_inf', type=int, default=150, help='number of initial infectious')
    parser.add_argument('--max_time', type=int, default=200, help='maximum duration of simulation')
    parser.add_argument('--removal_rate', type=float, default=0.2, help='[recovery rate + fatality rate]')
    parser.add_argument('--fatal_ratio', type=float, default=0.0037, help='ratio of fatality rate in the removal rate')
    parser.add_argument('--infect_rate', type=float, default=0.5, help='infected rate')
    parser.add_argument('--incubation', type=float, default=5.2, help='E take [incubation] days to convert to I')
    parser.add_argument('--global_ratio', type=float, default=0.5, help='ratio of global contact')
    config = parser.parse_args()

    manual_seed = config.manual_seed
    random.seed(manual_seed)
    np.random.seed(manual_seed)
    N = config.num_node
    I0 = config.num_init_inf
    t_max = config.max_time
    r = config.removal_rate
    mor_rate_portion = config.fatal_ratio
    mu = mor_rate_portion * r
    gamma = r - mu
    beta = config.infect_rate
    sigma = 1 / config.incubation
    p = config.global_ratio

    state_icon = ['S', 'E', 'I', 'R', 'F']
    base_graph = load_graph(name=config.base_graph)
    edge_list, node_list = list(base_graph.edges), np.array(list(base_graph.nodes))
    node_num, edge_num = len(node_list), len(edge_list)
    snnl = len(str(node_num))
    print(f'Graph stats: \nnode_num:{node_num}\tedge_num: {edge_num}')
    # init states
    state_list = np.zeros(shape=node_list.shape, dtype=np.int)
    init_infected_indices = np.random.randint(low=0, high=len(state_list), size=I0)
    state_list[init_infected_indices] = 2

    log_name = f'{config.log_name}_p{config.global_ratio}.txt'
    log_lines = []
    t_cur = 0
    P_num = 0
    while t_cur < t_max:
        rand_tuple = np.random.uniform(low=0, high=1, size=2)
        r1, r2 = rand_tuple[0], rand_tuple[1]
        S_node = np.where(state_list == 0)[0]
        E_node = np.where(state_list == 1)[0]
        I_node = np.where(state_list == 2)[0]
        R_node = np.where(state_list == 3)[0]
        F_node = np.where(state_list == 4)[0]
        S_num, E_num, I_num, R_num, F_num = len(S_node), len(E_node), len(I_node), len(R_node), len(F_node)
        A = N - F_num

        ip_global, ip_local, n_infect_ratio_list = trans_S2E(A, S_node, I_num)
        p_S2E = np.array([ip_global + ip_local * ni_ratio for ni_ratio in n_infect_ratio_list])
        p_E2I = np.array([sigma] * E_num)
        p_I2RF = np.array([gamma + mu] * I_num)
        prob_arr = np.concatenate([p_S2E, p_E2I, p_I2RF], axis=0)

        # Timing
        pi = np.sum(prob_arr)
        if pi == 0: break
        tau = np.log(1 / r1) / pi
        t_cur += tau

        # Sampling transition
        temp_node_list = np.concatenate([S_node, E_node, I_node], axis=0)
        # assert len(prob_arr) == len(node_list) == len(state_list) == len(temp_node_list)
        prob_arr = prob_arr / np.sum(prob_arr)
        chosen_node = np.random.choice(temp_node_list, replace=True, size=1, p=prob_arr)
        chosen_node_state = state_list[chosen_node]
        if chosen_node_state == 0:
            target_state = 1
        elif chosen_node_state == 1:
            target_state = 2
        else:
            target_state = 3 if r2 > mor_rate_portion else 4
        state_list[chosen_node] = target_state

        t_log = f'T:{t_cur:.3f}\tnS:{S_num}\tnE:{E_num}\tnI:{I_num}\tnP:{P_num}\tnR:{R_num}\t' \
                f'nF:{F_num}\tnode:{chosen_node}\t{state_icon[chosen_node_state[0]]} -> {state_icon[target_state]}'
        log_lines.append(f'{t_log}\n')
        print(f'\r{t_log}', end='')

    with open(f'log/{log_name}', 'w') as f:
        f.writelines(log_lines)
