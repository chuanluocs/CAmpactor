import os
import sys
import time
import random

def get_time_pid_random_string():
    my_time_str = time.strftime('%Y-%m-%d-%H:%M:%S', time.localtime(time.time()))
    my_pid = os.getpid()
    my_pid_str = str(my_pid)
    my_random = random.randint(1, sys.maxsize)
    my_random_str = str(my_random)
    my_time_pid_random_str = my_time_str + '_' + my_pid_str + '_' + my_random_str
    return my_time_pid_random_str

def get_num_way(model_file_path):
    fin = open(model_file_path)
    num_way = int(fin.readline().strip())
    fin.close()
    return num_way

def get_num_options(model_file_path):
    fin = open(model_file_path)
    fin.readline()
    num_options = int(fin.readline().strip())
    fin.close()
    return num_options

def get_num_constraints(constraint_file_path):
    fin = open(constraint_file_path)
    num_constraints = int(fin.readline().strip())
    fin.close()
    return num_constraints

def generate_list_schedule(pre_schedule_file_path, cutoff_time_pre):
    fin = open(pre_schedule_file_path, 'r')
    mylines = fin.readlines()
    myline = mylines[-1].strip()
    fin.close()
    
    list_schedule = eval(myline)
    list_schedule = [[item[0], item[1]*cutoff_time_pre+1] for item in list_schedule]
    return list_schedule

def generate_conf_file_adaptive(conf_file_path, num_way, num_options, num_constraints, cutoff_time_pre):
    if (num_way == 5 and num_options == 9 and num_constraints == 20):
        planned_time = cutoff_time_pre
        fout = open(conf_file_path, 'w+')
        myline = 'Conf_5 %d 1 0 0 4 1000' % (planned_time)
        fout.write(myline + '\n')
        fout.close()
        return
    elif (num_way == 5 and num_options == 15 and num_constraints == 13) or \
        (num_way == 5 and num_options == 15 and num_constraints == 48) or \
        (num_way == 5 and num_options == 18 and num_constraints == 13):
        planned_time = cutoff_time_pre
        fout = open(conf_file_path, 'w+')
        myline = 'Conf_6 %d 1 0 0 4 10000' % (planned_time)
        fout.write(myline + '\n')
        fout.close()
        return
    else:
        return

def generate_conf_file_from_list_schedule(list_schedule, conf_file_path, num_way, num_options, num_constraints, cutoff_time_pre):
    fout = open(conf_file_path, 'w+')
    for item in list_schedule:
        solver_name = item[0]
        planned_time = int(item[1]) + 1

        if solver_name == '<NEW_CHECK><LIGHTWEIGHT_SCORING><USE_GRADDEC><TABU_5><MODE_BALANCE_100000>':
            myline = 'Conf_1 %d 1 1 1 5 100000' % (planned_time)
        elif solver_name == '<NEW_CHECK><TRADITIONAL_SCORING><USE_GRADDEC><TABU_8><MODE_BALANCE_100000>':
            myline = 'Conf_2 %d 1 0 1 8 100000' % (planned_time)
        elif solver_name == '<NEW_CHECK><TRADITIONAL_SCORING><NO_GRADDEC><TABU_4><PURE_GREEDY>':
            myline = 'Conf_3 %d 1 0 0 4 -1' % (planned_time)
        elif solver_name == '<NEW_CHECK><LIGHTWEIGHT_SCORING><NO_GRADDEC><TABU_8><MODE_BALANCE_100000>':
            myline = 'Conf_4 %d 1 1 0 8 100000' % (planned_time)
        else:
            myline = 'Conf_1 %d 1 1 1 5 100000' % (planned_time)
    
        fout.write(myline + '\n')
    fout.close()
    generate_conf_file_adaptive(conf_file_path, num_way, num_options, num_constraints, cutoff_time_pre)
    return