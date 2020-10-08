import math
from scipy.interpolate import make_interp_spline
import matplotlib.pyplot as plt
import numpy as np


def get_xy(log_name, include_pos=False, include_qua=False):
    N = 10000
    with open(f'log/{log_name}.txt', 'r') as f:
        lines = f.readlines()
    x, y = [], []
    for line in lines:
        line = line.strip()
        line_comps = line.split()
        if include_qua:
            T_str, nS_str, nQ_str, nE_str, nI_str, nP_str, nR_str, nF_str, node_str, src_s, _, dst_s = line_comps
            T, nS, nQ, nE, nI, nP, nR, nF, node = float(T_str[2:]), int(nS_str[3:]), int(nQ_str[3:]), int(
                nE_str[3:]), int(nI_str[3:]), int(
                nP_str[3:]), int(
                nR_str[3:]), int(nF_str[3:]), node_str[
                                              6:10]
        else:
            if include_pos:
                T_str, nS_str, nE_str, nI_str, nP_str, nR_str, nF_str, node_str, src_s, _, dst_s = line_comps
            else:
                T_str, nS_str, nE_str, nI_str, nR_str, nF_str, node_str, src_s, _, dst_s = line_comps
                nP_str = 'nP:0'
            T, nS, nE, nI, nP, nR, nF, node = float(T_str[2:]), int(nS_str[3:]), int(nE_str[3:]), int(nI_str[3:]), int(
                nP_str[3:]), int(
                nR_str[3:]), int(nF_str[3:]), node_str[
                                              6:10]
        x.append(T)
        y.append(nI)
    # final_tau = x[-1] - x[-2]
    # cpl_num = math.floor((200 - x[-1]) / final_tau)
    # x.extend([x[-1] + final_tau * k for k in range(cpl_num)])
    # y.extend([y[-1]] * cpl_num)
    y = np.array(y) / N * 100
    return x, y


def draw_curve(x, y, label='', linestyle='', linewidth=3, marker='', marker_size=3):
    # remove duplicate in x
    non_dup_x, non_dup_y = [], []
    for i in range(len(x)):
        cur_x, cur_y = x[i], y[i]
        if cur_x not in non_dup_x:
            non_dup_x.append(cur_x)
            non_dup_y.append(cur_y)
    # interpolate curve
    x_num = int(max(non_dup_x))
    # x_smooth = np.linspace(start=0, stop=x_num, num=x_num)
    x_smooth = [k for k in range(x_num)]
    y_smooth = make_interp_spline(non_dup_x, non_dup_y)(x_smooth)
    plt.plot(x_smooth, y_smooth, label=label, linestyle=linestyle, linewidth=linewidth, marker=marker,
             markersize=marker_size, markevery=10)


def ni_tq_tqct(choice=0, mode='show'):
    policy_set1 = ['global transmission only', 'both global and network transmission', 'network transmission only']
    p_set = [1.0, 0.5, 0.0]
    p = p_set[choice]
    model_name = policy_set1[choice]
    save_dir = 'figure/NoIntervention VS TestingAndContactTracing'
    NI_x, NI_y = get_xy(log_name=f'NoIntervention_p{p}', include_pos=True)
    TQ_x, TQ_y = get_xy(log_name=f'Testing_Quarantine_p{p}_theta0.05', include_pos=True)
    TQC_x, TQC_y = get_xy(log_name=f'Testing_Quarantine_ContactTracing_p{p}_theta0.05_ph0.1', include_pos=True)
    plt.xlim(0, 200)
    plt.ylim(0, 12.5)
    plt.grid()
    ax = plt.gca()
    ax.yaxis.set_major_locator(plt.MultipleLocator(2))
    ax.xaxis.set_major_locator(plt.MultipleLocator(20))
    draw_curve(NI_x, NI_y, label='no intervention', linestyle='-')
    draw_curve(TQ_x, TQ_y, label='testing and quarantine', linestyle='--')
    draw_curve(TQC_x, TQC_y, label='testing, quarantine and contact tracing', linestyle=':')
    plt.legend()
    plt.xlabel('Time(days)')
    plt.ylabel('Infectious(%)')
    title_txt = f'NSIR-{model_name}, p={p}'
    plt.title(title_txt)
    if mode == 'save':
        plt.savefig(f'{save_dir}/{title_txt}.jpg')
    elif mode == 'show':
        plt.show()
    else:
        return


def distancing(timing='early', duration='short', mode='show'):
    timing_dict = {'early': 0, 'delayed': 30}
    duration_dict = {'short': 30, 'medium': 60, 'long': 120}
    start, dur = timing_dict[timing], duration_dict[duration]
    model_name = f'{timing} {duration} distancing'
    save_dir = 'figure/Distancing Policies'
    base_x, base_y = get_xy(log_name=f'Distancing_start{0}_duration{0}', include_pos=False)
    dist_x, dist_y = get_xy(log_name=f'Distancing_start{start}_duration{dur}', include_pos=False)
    plt.xlim(0, 300)
    plt.ylim(0, 12.5)
    plt.grid()
    ax = plt.gca()
    ax.yaxis.set_major_locator(plt.MultipleLocator(2))
    ax.xaxis.set_major_locator(plt.MultipleLocator(20))
    ax.add_patch(plt.Rectangle(xy=(start, 0), width=dur, height=12.5, color='maroon', alpha=0.5))
    draw_curve(base_x, base_y, label='no distancing', linestyle='--')
    draw_curve(dist_x, dist_y, label='distancing', linestyle='-')
    plt.legend()
    plt.xlabel('Time(days)')
    plt.ylabel('Infectious(%)')
    title_txt = f'{model_name}'
    plt.title(title_txt)
    if mode == 'save':
        plt.savefig(f'{save_dir}/{title_txt}.jpg')
    elif mode == 'show':
        plt.show()
    else:
        return


def lockdown(p=0.0, days=-1, mode='show'):
    model_name = 'SIR model (p=1)' if p == 1 else 'NSIR model (p=0)'
    save_dir = 'figure/Lockdown'
    base_x, base_y = get_xy(log_name=f'Lockdown_p{p}_lambda0.0_days{days}', include_pos=False)
    lock30_x, lock30_y = get_xy(log_name=f'Lockdown_p{p}_lambda0.3_days{days}', include_pos=False)
    lock70_x, lock70_y = get_xy(log_name=f'Lockdown_p{p}_lambda0.7_days{days}', include_pos=False)
    lock90_x, lock90_y = get_xy(log_name=f'Lockdown_p{p}_lambda0.9_days{days}', include_pos=False)
    plt.xlim(0, 200)
    plt.ylim(0, 12.5)
    plt.grid()
    ax = plt.gca()
    ax.yaxis.set_major_locator(plt.MultipleLocator(2))
    ax.xaxis.set_major_locator(plt.MultipleLocator(20))
    draw_curve(base_x, base_y, label='no lockdown', linestyle='-')
    draw_curve(lock30_x, lock30_y, label='30% lockdown', linestyle='--')
    draw_curve(lock70_x, lock70_y, label='70% lockdown', linestyle='-.')
    draw_curve(lock90_x, lock90_y, label='90% lockdown', linestyle=':')
    plt.legend()
    plt.xlabel('Time(days)')
    plt.ylabel('Infectious(%)')
    title_txt = f'{model_name}'
    plt.title(title_txt)
    if mode == 'save':
        plt.savefig(f'{save_dir}/{title_txt}.jpg')
    elif mode == 'show':
        plt.show()
    else:
        return


def reinforced_quarantine(p=.5, q_days=2, mode='show'):
    model_name = 'testing and quarantine'
    save_dir = 'figure/Our VS Original'
    ori_x, ori_y = get_xy(log_name=f'Testing_Quarantine_ContactTracing_p{p}_theta0.05_ph0.1', include_pos=True)
    strict_x, strict_y = get_xy(log_name=f'ReinforcedQuarantine_p{p}_days{q_days}', include_pos=True, include_qua=True)
    plt.xlim(0, 200)
    plt.ylim(0, 4.0)
    plt.grid()
    ax = plt.gca()
    ax.yaxis.set_major_locator(plt.MultipleLocator(2))
    ax.xaxis.set_major_locator(plt.MultipleLocator(20))
    draw_curve(ori_x, ori_y, label='original', linestyle=':')
    draw_curve(strict_x, strict_y, label='reinforced quarantine', linestyle='-')
    plt.legend()
    plt.xlabel('Time(days)')
    plt.ylabel('Infectious(%)')
    title_txt = f'NSIR-{model_name}, p={p}'
    plt.title(title_txt)
    if mode == 'save':
        plt.savefig(f'{save_dir}/{title_txt}.jpg')
    elif mode == 'show':
        plt.show()
    else:
        return


def sweep_q_days(p=0.5, mode='show'):
    model_name = 'reinforced quarantine'
    save_dir = 'figure/Our VS Original'
    x_day0, y_day0 = get_xy(log_name=f'ReinforcedQuarantine_p{p}_days0', include_pos=True, include_qua=True)
    x_day2, y_day2 = get_xy(log_name=f'ReinforcedQuarantine_p{p}_days2', include_pos=True, include_qua=True)
    x_day8, y_day8 = get_xy(log_name=f'ReinforcedQuarantine_p{p}_days8', include_pos=True, include_qua=True)
    x_day14, y_day14 = get_xy(log_name=f'ReinforcedQuarantine_p{p}_days14', include_pos=True, include_qua=True)
    plt.xlim(0, 200)
    plt.ylim(0, 6.0)
    plt.grid()
    ax = plt.gca()
    ax.yaxis.set_major_locator(plt.MultipleLocator(2))
    ax.xaxis.set_major_locator(plt.MultipleLocator(20))
    draw_curve(x_day0, y_day0, label='0 day quarantine', linestyle='-', marker='.')
    draw_curve(x_day2, y_day2, label='2 days quarantine', linestyle='-', marker='v')
    draw_curve(x_day8, y_day8, label='8 days quarantine', linestyle='-', marker='1')
    draw_curve(x_day14, y_day14, label='14 days quarantine', linestyle='-', marker='s')
    plt.legend()
    plt.xlabel('Time(days)')
    plt.ylabel('Infectious(%)')
    title_txt = f'NSIR-{model_name}, p={p}'
    plt.title(title_txt)
    if mode == 'save':
        plt.savefig(f'{save_dir}/{title_txt}.jpg')
    elif mode == 'show':
        plt.show()
    else:
        return


def sweep_contact_trace_rate(mode='save'):
    p = 0.5
    model_name = 'impact of contact tracing rate'
    save_dir = 'figure/Our VS Original'
    x01, y01 = get_xy(log_name=f'Testing_Quarantine_ContactTracing_p{p}_theta0.05_ph0.1', include_pos=True)
    x03, y03 = get_xy(log_name=f'Testing_Quarantine_ContactTracing_p{p}_theta0.05_ph0.3', include_pos=True)
    x06, y06 = get_xy(log_name=f'Testing_Quarantine_ContactTracing_p{p}_theta0.05_ph0.6', include_pos=True)
    x08, y08 = get_xy(log_name=f'Testing_Quarantine_ContactTracing_p{p}_theta0.05_ph0.8', include_pos=True)
    x10, y10 = get_xy(log_name=f'Testing_Quarantine_ContactTracing_p{p}_theta0.05_ph1.0', include_pos=True)
    plt.xlim(0, 200)
    plt.ylim(0, 3.5)
    plt.grid()
    ax = plt.gca()
    ax.yaxis.set_major_locator(plt.MultipleLocator(2))
    ax.xaxis.set_major_locator(plt.MultipleLocator(20))
    draw_curve(x01, y01, label='ph0.1', linestyle='-', marker='', linewidth=1)
    draw_curve(x03, y03, label='ph0.3', linestyle='-', marker='o', linewidth=1)
    draw_curve(x06, y06, label='ph0.6', linestyle='-', marker='*', linewidth=1)
    draw_curve(x08, y08, label='ph0.8', linestyle='-', marker='>', linewidth=1)
    draw_curve(x10, y10, label='ph1.0', linestyle='-', marker='+', linewidth=1)
    plt.legend()
    plt.xlabel('Time(days)')
    plt.ylabel('Infectious(%)')
    title_txt = f'NSIR-{model_name}, p={p}'
    plt.title(title_txt)
    if mode == 'save':
        plt.savefig(f'{save_dir}/{title_txt}.jpg')
    elif mode == 'show':
        plt.show()
    else:
        return


if __name__ == "__main__":
    # reinforced_quarantine(p=0.5, mode='show')
    sweep_q_days(p=0.0, mode='save')
    # tqct_and_strict_tqct(choice='1')
    # sweep_contact_trace_rate(mode='save')
    # tqct_and_strict_tqct(mode='save', choice=0)
    # tqct_and_strict_tqct(mode='save', choice=1)
    # tqct_and_strict_tqct(mode='save', choice=2)
    # distancing(timing='early', duration='short', mode='save')
    # distancing(timing='early', duration='medium', mode='save')
    # distancing(timing='early', duration='long', mode='save')
    # distancing(timing='delayed', duration='short', mode='save')
    # distancing(timing='delayed', duration='medium', mode='save')
    # distancing(timing='delayed', duration='long', mode='save')
    # lockdown(p=1.0, days=-1, mode='show')
