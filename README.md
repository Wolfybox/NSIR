# NSIR  
An unofficial implementation of paper "Karaivanov, Alexander. (2020). A Social Network Model of COVID-19"  
Paper Link: (https://papers.ssrn.com/sol3/papers.cfm?abstract_id=3584895)  

# Note  
This is a course project for HITSZ Modeling of Complex Networks Semester A 2020.   
# Features  
We have reproduced part of the results (or scenarios simulation) from the paper, which includes:  
1. Herd Immunity  
2. Testing and Quarantine  
3. Contact Tracing  
4. Distancing  
5. Lockdown  

Additionally, we proposed two model variations by considering two real-world pandemic policies in China:  
1. Strict Contact Tracing, where all neighbors of an infectious agents will have to go through a test  
2. Quarantine with contact tracing at certain duration, where neighbors of an infectious agents being tested positive will be sent into finite-long quarantine. (usually 14 days)  

The above two policies are also named as Reinforced Contact tracing and Quarantine. For implementation details, please see codes.  
# Dependencies   
|type|name|version|  
|-----|-----|-----|  
|os|Ubuntu or Windows|16.04(Ubuntu); Win10(Windows)|  
|language|python|3.8|  
|networktools|networkx|2.5|  
|calculation|numpy|1.19.1|  
# Usage  
Run command in terminal:  
```shell  
$ cd project_dir  
$ python  NoIntervention.py [--max_time 200] [--removal_rate 0.2] [--fatal_ratio 0.0037] [--infect_rate 0.5] [--incubation 5.2]  [--global_ratio 0.5] [--num_init_inf 150]  
```    
