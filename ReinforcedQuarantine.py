import argparse
import math
import random

import numpy as np

from network import load_graph


def get_non_q_S_ngr_stats(node_i):
    base_neighbors = np.array([i for i in base_graph.neighbors(node_i)])
    close_neighbors = np.array([i for i in close_graph.neighbors(node_i)])
    # remove quarantine neighbors (does not affect close neighbors)
    q_nodes = np.where(q_flags == 1)[0]
    base_neighbors = np.setdiff1d(base_neighbors, q_nodes)
    n_all_adj = len(base_neighbors)
    if not n_all_adj: return .0
    base_neighbors_states = state_list[base_neighbors]
    I_neighbors = base_neighbors[np.where(base_neighbors_states == state_icon['I'])[0]]
    if len(close_neighbors):
        close_neighbors_states = state_list[close_neighbors]
        P_neighbors = close_neighbors[np.where(close_neighbors_states == state_icon['P'])[0]]
    else:
        P_neighbors = []
    n_I_adj, n_P_adj = len(I_neighbors), len(P_neighbors)
    return (n_I_adj + n_P_adj) / n_all_adj


def get_q_S_ngr_stats(node_i):
    close_neighbors = np.array([i for i in close_graph.neighbors(node_i)])
    if not len(close_neighbors): return .0
    neighbors_states = state_list[close_neighbors]
    local_P_neighbors = len(np.where(neighbors_states == state_icon['P'])[0])
    local_I_neighbors = len(np.where(neighbors_states == state_icon['I'])[0])
    local_all_neighbors = len(close_neighbors)
    return (local_P_neighbors + local_I_neighbors) / local_all_neighbors


def get_ct_ngr(node_i):
    base_neighbors = np.array([i for i in base_graph.neighbors(node_i)])
    close_neighbors = np.array([i for i in close_graph.neighbors(node_i)])
    if q_flags[node_i]:
        if len(close_neighbors) == 0: return 0
        neighbors_states = state_list[close_neighbors]
        P_neighbors = len(np.where(neighbors_states == state_icon['P'])[0])
    else:
        if len(base_neighbors) == 0: return 0
        neighbors_states = state_list[base_neighbors]
        P_neighbors = len(np.where(neighbors_states == state_icon['P'])[0])
    return P_neighbors


def trans_S2E():
    I_non_q_num = len(np.where(q_flags[I_node] == 0)[0])
    global_trans_prob = p * beta * I_non_q_num / A
    coeff_global_ct = (1 - p) * beta
    coeff_local_ct = beta
    p_S2E = []
    for node in S_node:
        if q_flags[node]:
            n_stats = get_q_S_ngr_stats(node)
            p_S2E.append(coeff_local_ct * n_stats)
        else:
            n_stats = get_non_q_S_ngr_stats(node)
            p_S2E.append(global_trans_prob + coeff_global_ct * n_stats)
    return np.array(p_S2E)


def in_quarantine(infect_node):
    base_neighbors = np.array([i for i in base_graph.neighbors(infect_node[0])])
    base_neighbors = np.setdiff1d(base_neighbors, P_node)
    ngr_q_states = q_flags[base_neighbors]
    non_q_ngr = base_neighbors[np.where(ngr_q_states == 0)[0]]
    # for non-quarantine node, initialize duration
    q_counters[non_q_ngr] = q_days
    # set all target neighbors into quarantine states
    q_flags[base_neighbors] = 1


def out_quarantine():
    Q_node = np.where(q_flags == 1)[0]
    q_expired = np.where(q_counters <= 0)[0]
    release_Q = np.intersect1d(Q_node, q_expired)
    q_flags[release_Q] = 0


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--base_graph', default=r'baseline', help='name of base graph model')
    parser.add_argument('--close_graph', default=r'close_contact', help='name of close contact graph model')
    parser.add_argument('--log_name', default=r'ReinforcedQuarantine', help='name of log of current simulation')
    parser.add_argument('--manual_seed', type=int, help='manual random seed', default=0)
    parser.add_argument('--num_node', type=int, default=10000, help='number of nodes in both graphs')
    parser.add_argument('--num_init_inf', type=int, default=150, help='number of initial infectious')
    parser.add_argument('--max_time', type=int, default=500, help='maximum duration of simulation')
    parser.add_argument('--removal_rate', type=float, default=0.2, help='[recovery rate + fatality rate]')
    parser.add_argument('--fatal_ratio', type=float, default=0.0037, help='ratio of fatality rate in the removal rate')
    parser.add_argument('--infect_rate', type=float, default=0.5, help='infected rate')
    parser.add_argument('--incubation', type=float, default=5.2, help='E take [incubation] days to convert to I')
    parser.add_argument('--mass_testing_rate', type=float, default=0.05, help='ratio of massive testing')
    parser.add_argument('--contact_trace_rate', type=float, default=0.0, help='ratio of contact tracing')
    parser.add_argument('--global_ratio', type=float, default=0.5, help='ratio of global contact')
    parser.add_argument('--quarantine_days', type=int, default=14, help='duration of quarantine')
    config = parser.parse_args()

    manual_seed = config.manual_seed
    random.seed(manual_seed)
    np.random.seed(manual_seed)
    N = config.num_node
    I0 = config.num_init_inf
    t_max = config.max_time
    q_days = config.quarantine_days

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
    ph = config.contact_trace_rate
    p = config.global_ratio

    state_icon = {'S': 0, 'E': 1, 'I': 2, 'P': 3, 'R': 4, 'F': 5}
    base_graph = load_graph(name=config.base_graph)
    close_graph = load_graph(name=config.close_graph)
    edge_list, node_list = list(base_graph.edges), np.array(list(base_graph.nodes))
    cc_edge_list, cc_node_list = list(close_graph.edges), np.array(list(close_graph.nodes))

    node_num, edge_num = len(node_list), len(edge_list)
    cc_node_num, cc_edge_num = len(cc_node_list), len(cc_edge_list)
    print(f'Base Graph stats: \nnode_num:{node_num}\tedge_num: {edge_num}')
    print(f'Close Graph stats: \nnode_num:{cc_node_num}\tedge_num: {cc_edge_num}')

    # init states
    state_list = np.zeros(shape=node_list.shape, dtype=np.int)
    init_infected_indices = np.random.randint(low=0, high=len(state_list), size=I0)
    state_list[init_infected_indices] = state_icon['I']

    # quarantine days counter
    q_counters = np.zeros(shape=node_list.shape, dtype=np.int)
    # quarantine flags
    # q_init_indices = np.random.randint(low=0, high=len(node_list), size=8000)
    q_flags = np.zeros(shape=node_list.shape, dtype=np.int)
    # q_flags[q_init_indices] = 1
    # q_counters[q_init_indices] = q_days

    # static trans prob arrays
    P_trans_prob = np.array([gamma_p, mu_p])
    P_trans_prob = P_trans_prob / np.sum(P_trans_prob)

    log_name = f'{config.log_name}_p{config.global_ratio}_days{q_days}.txt'
    log_lines = []
    t_cur, t_prev = 0, 0
    while t_cur < t_max:
        # quarantine days countdown
        if int(t_cur) - int(t_prev) == 1:
            q_counters = q_counters - 1
        t_prev = t_cur
        # sketch statistics
        r1 = np.random.uniform(low=0, high=1, size=1)
        S_node = np.where(state_list == state_icon['S'])[0]
        Q_node = np.where(q_flags == 1)[0]
        E_node = np.where(state_list == state_icon['E'])[0]
        I_node = np.where(state_list == state_icon['I'])[0]
        P_node = np.where(state_list == state_icon['P'])[0]
        R_node = np.where(state_list == state_icon['R'])[0]
        F_node = np.where(state_list == state_icon['F'])[0]
        S_num, Q_num, E_num, I_num, P_num, R_num, F_num = len(S_node), len(Q_node), len(E_node), len(I_node), len(
            P_node), len(R_node), len(
            F_node)
        A = N - F_num - P_num - Q_num

        # Susceptible transition
        p_S2E = trans_S2E()

        # Exposed transition
        p_E2I = np.array([sigma] * E_num)

        # Infectious Transition
        p_I2P = np.array([min([1.0, theta + ph * get_ct_ngr(node)]) for node in I_node])
        p_I2RF = np.array([gamma + mu] * I_num)
        p_I2PRF = p_I2RF + p_I2P

        # Testing transition
        p_P2RF = np.array([gamma_p + mu_p] * P_num)

        prob_arr = np.concatenate([p_S2E, p_E2I, p_I2PRF, p_P2RF], axis=0)

        # probability for transition into specific states
        I_trans_prob = np.concatenate([p_I2P, np.array([gamma] * I_num), np.array([mu] * I_num)], axis=0)
        I_trans_prob = I_trans_prob / np.sum(I_trans_prob)

        # Timing
        pi = np.sum(prob_arr)
        if pi == 0: break
        tau = np.log(1 / r1) / pi
        t_cur += tau[0]

        # Sampling transition
        temp_node_list = np.concatenate([S_node, E_node, I_node, P_node], axis=0)
        prob_arr = prob_arr / np.sum(prob_arr)
        chosen_node = np.random.choice(temp_node_list, replace=True, size=1, p=prob_arr)

        chosen_node_state = state_list[chosen_node][0]
        if chosen_node_state == state_icon['S']:
            target_state = state_icon['E']
        elif chosen_node_state == state_icon['E']:
            target_state = state_icon['I']
        elif chosen_node_state == state_icon['I']:
            target_state = \
                np.random.choice([state_icon['P']] * I_num + [state_icon['R']] * I_num + [state_icon['F']] * I_num,
                                 replace=True, size=1, p=I_trans_prob)[0]
            # trans all of its valid neighbors into quarantine
            # if one infectious case is found
            if target_state == state_icon['P']:
                in_quarantine(chosen_node)
        else:
            target_state = np.random.choice([state_icon['R'], state_icon['F']], replace=True, size=1, p=P_trans_prob)[0]
        out_quarantine()  # check remained quarantine days and released node reaching deadline
        state_list[chosen_node] = target_state

        t_log = f'T:{t_cur:.3f}\tnS:{S_num}\tnQ:{Q_num}\tnE:{E_num}\tnI:{I_num}\tnP:{P_num}\tnR:{R_num}\t' \
                f'nF:{F_num}\tnode:{chosen_node}\t{list(state_icon.keys())[chosen_node_state]} -> {list(state_icon.keys())[target_state]}'
        log_lines.append(f'{t_log}\n')
        print(f'\r{t_log}', end='')

    with open(f'log/{log_name}', 'w') as f:
        f.writelines(log_lines)
