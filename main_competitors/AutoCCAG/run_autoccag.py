#!/usr/bin/env python3

import os
import sys
from autoccag.autoccag import get_time_pid_random_string, get_num_way, get_num_options, get_num_constraints, generate_list_schedule, generate_conf_file_from_list_schedule

if __name__ == '__main__':
    
    if len(sys.argv) != 6:
        print('c Usage: python3 %s <model_file> <constraint_file> <seed> <cutoff_time> <ca_file>' % (sys.argv[0]))
        sys.exit(-1)
    
    model_file_path = os.path.abspath(sys.argv[1])
    constraint_file_path = os.path.abspath(sys.argv[2])
    seed = int(sys.argv[3])
    cutoff_time = int(sys.argv[4])
    ca_file = sys.argv[5]
    
    time_ratio = 1.0
    cutoff_time_pre = int(cutoff_time * time_ratio) + 1
    
    num_way = get_num_way(model_file_path)
    num_options = get_num_options(model_file_path)
    num_constraints = get_num_constraints(constraint_file_path)
    
    tmp_dir = 'tmp/'
    if not os.path.isdir(tmp_dir):
        os.system('mkdir -p ' + tmp_dir)
    os.system("rm " + tmp_dir + "*.txt")
    
    key_str = get_time_pid_random_string()
    conf_file_name = 'conf_' + key_str + '.txt'
    conf_file_path = tmp_dir + conf_file_name

    pre_schedule_file_path = 'data/pre_solver_schedule.txt'
    list_schedule = generate_list_schedule(pre_schedule_file_path, cutoff_time_pre)
    generate_conf_file_from_list_schedule(list_schedule, conf_file_path, num_way, num_options, num_constraints, cutoff_time_pre)
    
    command = 'cd MetaTCA/; ./TCA -model_file %s -constraints_file %s -configuration_file ../%s -seed %d -ca_file %s -final_check 1 2> /dev/null' % (model_file_path, constraint_file_path, conf_file_path, seed, ca_file)
    os.system(command)
