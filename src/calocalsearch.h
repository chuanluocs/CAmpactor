#ifndef LOCALSEARCH_OPTIMIZER_INCLUDE_H
#define LOCALSEARCH_OPTIMIZER_INCLUDE_H

#include "core/Solver.h"
#include "clihelper.h"
#include <vector>
#include <numeric>
#include <random>
#include <utility>
#include <algorithm>
#include <set>
#include <map>

using std::vector;
using std::pair;
using std::set;
using std::map;

class LSOptimizer {
public:
    LSOptimizer(const Argument& params);
    ~LSOptimizer();

    void search();
    vector<vector<int> > get_testcase_set();
private:
    int verbosity;
    
    int nvar;
    int nclauses;
    int testcase_size;
    long long num_combination_all_possible;
    long long num_tuple_all_possible;
   
    vector<int> tuple_cov_cnt;

    vector<vector<int> > clauses;
    vector<vector<int> > pos_in_cls;
    vector<vector<int> > neg_in_cls;

    vector<vector<int> > clauses_cov;

    vector<vector<int> > ca;
    vector<vector<int> > last_ca;

    vector<pair<int, int> > uncovered_tuples;
    vector<int> uncovered_tuple_positions;

    std::mt19937_64 gen;
    int limit;
    
    int testcase_taboo;
    int greedy_limit;
    vector<int> last_greedy_time;

    vector<pair<int, int> > tmp_break, tmp_gain;
    
    vector<int> tmp_flip_sequence;

    int (LSOptimizer::*remove_testcase_strategy)();

    int remove_testcase_randomly();
    int remove_testcase_greedily();

    bool random_greedy_step();
    void greedy_step_forced(int comp1, int comp2);
    void random_step();

    Minisat::Solver *cdcl_solver;

    long long get2tupleid(int i, int vi, int j, int vj);

    long long new_uncovered_tuples_after_remove_testcase(const vector<int>& tc, bool dump);

    pair<bool, pair<int, int> > get_gain_for_flip(int tcid, int vid, bool dump);
    void flip_bit(int tcid, int vid);

    pair<bool, pair<int, int> > get_gain_for_forcetuple(int tcid, int comp1, int comp2, bool dump);
    void forcetuple(int tcid, int comp1, int comp2);

    void break_and_gain();
    void update_clause_cov_by_flip(int tcid, int vid, int curbit, bool trailing);

    void validate();
    void final_validate();
    bool validate_testcase(const vector<int>& tc);

    pair<int, int> get_gain_for_forcetestcase(int tcid, const vector<int>& tc2, bool dump);
    void forcetestcase(int tcid, const vector<int>& tc2);

    bool __enable_random_step;
    bool __enable_greedy_step;
    int __forced_greedy_percent;
    int __stop_length;
    int __last_success_step;
    
    vector<vector<int> > last_flipped_time;
    bool __use_cell_tabu;
};

#endif
