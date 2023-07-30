#ifndef CLIHELPER_INCLUDE_H
#define CLIHELPER_INCLUDE_H

#include <cstdio>
#include <string>
#include <unordered_map>
#include <utility>
#include <functional>

using std::string;
using std::unordered_map;
using std::pair;

enum Opt_type { Opt_toggle, Opt_int, Opt_string, Opt_double };

static void simple_int_input(int* dst, const char* src){
    sscanf(src, "%d", dst);
}

static void simple_toggle_input(bool* dst){
    (*dst) = !(*dst);
}

static void simple_string_input(string& dst, const char* src){
    dst = string(src);
}

struct Argument {
    bool flag_input_cnf_path;
    string input_cnf_path;
    bool flag_reduced_cnf_file_path;
    string reduced_cnf_file_path;
    bool flag_output_testcase_path;
    string output_testcase_path;
    bool flag_use_existing_ca;
    string init_ca_path;

    int seed;
    bool flag_use_cnf_reduction;
    int given_limit;
    int testcase_taboo;
    int verbosity;
    int forced_greedy_percent;
    bool flag_use_random_step;
    bool flag_use_greedy_step;
    bool flag_use_testcase_taboo;
    int stop_length;
    bool flag_use_cell_tabu;

    int twise;

    Argument():
        seed(1), 
        flag_input_cnf_path(false), 
        flag_reduced_cnf_file_path(false), 
        flag_output_testcase_path(false), 
        flag_use_existing_ca(false),
        flag_use_cnf_reduction(true),
        given_limit(10), 
        testcase_taboo(10), 
        verbosity(0), 
        forced_greedy_percent(10), 
        flag_use_random_step(true), 
        flag_use_greedy_step(true), 
        flag_use_testcase_taboo(true), 
        stop_length(10000), 
        flag_use_cell_tabu(false), 
        twise(2)
    {}
};

typedef unordered_map<string, pair<Opt_type, std::function<void(Argument*, const char*)> > > action_table;

static action_table build_default_actions(){
    action_table res;

    res["-input_cnf_path"] = std::make_pair(Opt_string, [](Argument* argu, const char* v){ 
        simple_string_input(argu->input_cnf_path, v); 
        argu->flag_input_cnf_path = true;
    });
    res["-i"] = res["-input_cnf_path"];
    res["-output_testcase_path"] = std::make_pair(Opt_string, [](Argument* argu, const char* v){ 
        simple_string_input(argu->output_testcase_path, v); 
        argu->flag_output_testcase_path = true;
    });
    res["-o"] = res["-output_testcase_path"];
    res["-reduced_cnf_path"] = std::make_pair(Opt_string, [](Argument* argu, const char* v){ 
        simple_string_input(argu->reduced_cnf_file_path, v); 
        argu->flag_reduced_cnf_file_path = true;
    });
    res["-initwith"] = std::make_pair(Opt_string, [](Argument* argu, const char* v){ 
        simple_string_input(argu->init_ca_path, v); 
        argu->flag_use_existing_ca = true;
    });
    res["--init_PCA_path"] = res["-initwith"];

    res["-seed"] = std::make_pair(Opt_int, [](Argument* argu, const char* v){ 
        simple_int_input(&(argu->seed), v); 
    });
    res["--seed"] = res["-seed"];
    res["-tc_taboo"] = std::make_pair(Opt_int, [](Argument* argu, const char* v){ 
        simple_int_input(&(argu->testcase_taboo), v); 
    });
    res["--delta"] = res["-tc_taboo"];
    res["-verbose"] = std::make_pair(Opt_int, [](Argument* argu, const char* v){ 
        simple_int_input(&(argu->verbosity), v); 
    });
    res["-v"] = res["-verbose"];
    res["-forced_greedy"] = std::make_pair(Opt_int, [](Argument* argu, const char* v){ 
        double tmp;
        sscanf(v, "%lf", &tmp);
        argu->forced_greedy_percent = tmp * 100;
    });
    res["--psi"] = res["-forced_greedy"];
    res["-stop_length"] = std::make_pair(Opt_int, [](Argument* argu, const char* v){ 
        simple_int_input(&(argu->stop_length), v); 
    });
    res["--gamma"] = res["-stop_length"];

    res["--twise"] = std::make_pair(Opt_int, [](Argument* argu, const char* v){ 
        simple_int_input(&(argu->twise), v); 
    });
    
    res["-simplcnf"] = std::make_pair(Opt_toggle, [](Argument* argu, const char* v){ 
        simple_toggle_input(&(argu->flag_use_cnf_reduction)); 
    });
    res["--nosimplcnf"] = res["-simplcnf"];
    res["-no_random_step"] = std::make_pair(Opt_toggle, [](Argument* argu, const char* v){ 
        simple_toggle_input(&(argu->flag_use_random_step)); 
    });
    res["-no_greedy_step"] = std::make_pair(Opt_toggle, [](Argument* argu, const char* v){ 
        simple_toggle_input(&(argu->flag_use_greedy_step)); 
    });
    res["-no_testcase_taboo"] = std::make_pair(Opt_toggle, [](Argument* argu, const char* v){ 
        simple_toggle_input(&(argu->flag_use_testcase_taboo)); 
    });
    res["-use_cell_tabu"] = std::make_pair(Opt_toggle, [](Argument* argu, const char* v){ 
        simple_toggle_input(&(argu->flag_use_cell_tabu)); 
    });
    res["--use_cell_tabu"] = res["-use_cell_tabu"];

    return res;
}

static bool match_by_action_table(const action_table& act, int argc, char **argv, Argument* argu){
    if (argc < 2) return false;
    for (int i = 1; i < argc; ++i){
        auto iter = act.find(string(argv[i]));
        if (iter == act.end()){
            return false;
        }

        if ((iter->second).first != Opt_toggle){
            ++i;
            if (i >= argc) return false;
        }
        ((iter->second).second)(argu, argv[i]);
    }

    return argu->flag_input_cnf_path && argu->flag_use_existing_ca && argu->flag_output_testcase_path;
}

#endif
