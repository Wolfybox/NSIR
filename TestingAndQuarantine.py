import argparse
import math
import random
from network import load_graph
import networkx as nx
import matplotlib.pyplot as plt
import numpy as np
import pickle


def get_inf_pos_n_ratio(node_i):
    base_neighbors = np.array([i for i in base_graph.neighbors(node_i)])
    close_neighbors = np.array([i for i in close_graph.neighbors(node_i)])
    n_all_adj = len(base_neighbors)
    if not n_all_adj:
        return .0
    base_neighbors_states = state_list[base_neighbors]
    I_neighbors = base_neighbors[np.where(base_neighbors_states == 2)]
    if len(close_neighbors):
        close_neighbors_states = state_list[close_neighbors]
        P_neighbors = close_neighbors[np.where(close_neighbors_states == 3)]
    else:
        P_neighbors = []
    n_I_adj, n_P_adj = len(I_neighbors), len(P_neighbors)
    return (n_I_adj + n_P_adj) / n_all_adj


def trans_S2E(A, S_node, I_num):
    ip_global = p * beta * I_num / A
    inf_n_ratio_list = [get_inf_pos_n_ratio(node) for node in S_node]
    ip_local = (1 - p) * beta
    return ip_global, ip_local, inf_n_ratio_list


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--base_graph', default=r'baseline', help='name of base graph model')
    parser.add_argument('--close_graph', default=r'close_contact', help='name of close contact graph model')
    parser.add_argument('--log_name', default=r'Testing_Quarantine', help='name of log of current simulation')
    parser.add_argument('--manual_seed', type=int, help='manual random seed', default=0)
    parser.add_argument('--num_node', type=int, default=10000, help='number of nodes in both graphs')
    parser.add_argument('--num_init_inf', type=int, default=150, help='number of initial infectious')
    parser.add_argument('--max_time', type=int, default=200, help='maximum duration of simulation')
    parser.add_argument('--removal_rate', type=float, default=0.2, help='[recovery rate + fatality rate]')
    parser.add_argument('--fatal_ratio', type=float, default=0.0037, help='ratio of fatality rate in the removal rate')
    parser.add_argument('--infect_rate', type=float, default=0.5, help='infected rate')
    parser.add_argument('--incubation', type=float, default=5.2, help='E take [incubation] days to convert to I')
    parser.add_argument('--mass_testing_rate', type=float, default=0.05, help='ratio of massive testing')
    parser.add_argument('--global_ratio', type=float, default=0.5, help='ratio of global contact')
    config = parser.parse_args()

    manual_seed = config.manual_seed
    random.seed(manual_seed)
    np.random.seed(manual_seed)

    N = config.num_node
    I0 = config.num_init_inf
    t_max = config.max_time

    r = config.removal_rate
    gamma_mu_rate = config.fatal_ratio
    mu = gamma_mu_rate * r
    gamma = r - mu

    gamma_mu_rate_p = config.fatal_ratio
    mu_p = gamma_mu_rate_p * r
    gamma_p = r - mu_p

    beta = config.infect_rate
    sigma = 1 / config.incubation
    theta = config.mass_testing_rate
    p = config.global_ratio

    baseline_name = config.base_graph
    close_contact_name = config.close_graph

    state_icon = ['S', 'E', 'I', 'P', 'R', 'F']
    base_graph = load_graph(name=baseline_name)
    close_graph = load_graph(name=close_contact_name)
    edge_list, node_list = list(base_graph.edges), np.array(list(base_graph.nodes))
    cc_edge_list, cc_node_list = list(close_graph.edges), np.array(list(close_graph.nodes))

    node_num, edge_num = len(node_list), len(edge_list)
    cc_node_num, cc_edge_num = len(cc_node_list), len(cc_edge_list)
    print(f'Base Graph stats: \nnode_num:{node_num}\tedge_num: {edge_num}')
    print(f'Close Graph stats: \nnode_num:{cc_node_num}\tedge_num: {cc_edge_num}')

    # init states
    state_list = np.zeros(shape=node_list.shape, dtype=np.int)
    init_infected_indices = np.random.randint(low=0, high=len(state_list), size=I0)
    state_list[init_infected_indices] = 2

    # Init prob arrays
    I_trans_prob = np.array([theta, gamma, mu])
    I_trans_prob = I_trans_prob / np.sum(I_trans_prob)
    P_trans_prob = np.array([gamma_p, mu_p])
    P_trans_prob = P_trans_prob / np.sum(P_trans_prob)

    log_name = f'{config.log_name}_p{config.global_ratio}_theta{config.mass_testing_rate}.txt'
    log_lines = []
    t_cur = 0
    while t_cur < t_max:
        r1 = np.random.uniform(low=0, high=1, size=1)
        S_node = np.where(state_list == 0)[0]
        E_node = np.where(state_list == 1)[0]
        I_node = np.where(state_list == 2)[0]
        P_node = np.where(state_list == 3)[0]
        R_node = np.where(state_list == 4)[0]
        F_node = np.where(state_list == 5)[0]
        S_num, E_num, I_num, P_num, R_num, F_num = len(S_node), len(E_node), len(I_node), len(P_node), len(R_node), len(
            F_node)
        A = N - F_num - P_num

        ip_global, ip_local, n_infect_pos_ratio_list = trans_S2E(A, S_node, I_num)
        p_S2E = np.array([ip_global + ip_local * nip_ratio for nip_ratio in n_infect_pos_ratio_list])
        p_E2I = np.array([sigma] * E_num)
        p_I2PRF = np.array([gamma + theta + mu] * I_num)
        p_P2RF = np.array([gamma_p + mu_p] * P_num)
        prob_arr = np.concatenate([p_S2E, p_E2I, p_I2PRF, p_P2RF], axis=0)

        # Timing
        pi = np.sum(prob_arr)
        if pi == 0: break  # Transition Stop
        tau = np.log(1 / r1) / pi
        t_cur += tau[0]

        # Sampling transition
        temp_node_list = np.concatenate([S_node, E_node, I_node, P_node], axis=0)
        prob_arr = prob_arr / np.sum(prob_arr)
        chosen_node = np.random.choice(temp_node_list, replace=True, size=1, p=prob_arr)
        chosen_node_state = state_list[chosen_node][0]
        if chosen_node_state == 0:
            target_state = 1
        elif chosen_node_state == 1:
            target_state = 2
        elif chosen_node_state == 2:
            target_state = np.random.choice([3, 4, 5], replace=True, size=1, p=I_trans_prob)[0]
        else:
            target_state = np.random.choice([4, 5], replace=True, size=1, p=P_trans_prob)[0]
        state_list[chosen_node] = target_state

        t_log = f'T:{t_cur:.3f}\tnS:{S_num}\tnE:{E_num}\tnI:{I_num}\tnP:{P_num}\tnR:{R_num}\t' \
                f'nF:{F_num}\tnode:{chosen_node}\t{state_icon[chosen_node_state]} -> {state_icon[target_state]}'
        log_lines.append(f'{t_log}\n')
        print(f'\r{t_log}', end='')

    with open(f'log/{log_name}', 'w') as f:
        f.writelines(log_lines)
