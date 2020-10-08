import os

if __name__ == '__main__':
    # policies = ['NoIntervention', 'TestingAndQuarantine', 'TQAndContactTracing', 'Distancing', 'Lockdown']
    policies = ['ReinforcedQuarantine']
    dash_len = 16
    print(f'{"=" * dash_len * 2} Simulation Start {"=" * dash_len * 2}')
    for policy in policies:
        if policy == 'NoIntervention':
            for p in [0.0, 0.5, 1.0]:
                script = f'python {policy}.py --global_ratio {p}'
                print(f'\n{"-" * dash_len} Running {policy} with p={p} {"-" * dash_len}')
                os.system(script)
        if policy == 'TestingAndQuarantine':
            theta = 0.05
            for p in [0.0, 0.5, 1.0]:
                script = f'python {policy}.py --global_ratio {p} --mass_testing_rate {theta}'
                print(f'\n{"-" * dash_len} Running {policy} with p={p} theta={theta} {"-" * dash_len}')
                os.system(script)
        if policy == 'TQAndContactTracing':
            theta, ph = 0.05, 0.1
            for p in [0.0, 0.5, 1.0]:
                script = f'python {policy}.py --global_ratio {p} --mass_testing_rate {theta} --contact_trace_rate {ph}'
                print(f'\n{"-" * dash_len} Running {policy} with p={p} theta={theta} ph={ph} {"-" * dash_len}')
                os.system(script)
        if policy == 'Distancing':
            for dtc_start in [0, 30]:
                for dtc_duration in [0, 30, 60, 120]:
                    script = f'python {policy}.py --dtc_start {dtc_start} --dtc_duration {dtc_duration}'
                    print(f'\n{"-" * dash_len} Running {policy} with s={dtc_start} dur={dtc_duration} {"-" * dash_len}')
                    os.system(script)
        if policy == 'Lockdown':
            lockdown_days = -1
            for p in [0.0, 1.0]:
                for lockdown_rate in [0.0, 0.3, 0.7, 0.9]:
                    script = f'python {policy}.py --global_ratio {p} --lockdown_rate {lockdown_rate} --lockdown_days {lockdown_days}'
                    print(
                        f'\n{"-" * dash_len} Running {policy} with p={p} lock_rate={lockdown_rate} lock_days={lockdown_days} {"-" * dash_len}')
                    os.system(script)
        if policy == 'ReinforcedQuarantine':
            theta, ph = 0.05, 0.1
            for p in [0.0, 0.5]:
                for q_days in [0]:
                    script = f'python {policy}.py --global_ratio {p} --mass_testing_rate {theta} --contact_trace_rate {ph} --quarantine_days {q_days}'
                    print(
                        f'\n{"-" * dash_len} Running {policy} with p={p} theta={theta} ph={ph} q_days={q_days} {"-" * dash_len}')
                    os.system(script)
