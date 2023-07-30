#include "calocalsearch.h"
#include "cnfinfo.h"
#include <iostream>
#include <unistd.h>
#include <chrono>

LSOptimizer::LSOptimizer(const Argument& params):
    verbosity(params.verbosity),
    gen(params.seed), 
    limit(params.given_limit), 
    testcase_taboo(params.testcase_taboo), 
    twise(params.twise),
    __forced_greedy_percent(params.forced_greedy_percent),
    __enable_random_step(params.flag_use_random_step), 
    __enable_greedy_step(params.flag_use_greedy_step),
    __stop_length(params.stop_length), 
    __use_cell_tabu(params.flag_use_cell_tabu)
{
    std::string real_input_cnf_path = params.input_cnf_path;
    if (params.flag_use_cnf_reduction){
        std::mt19937_64 tmp_gen(params.seed);
        std::string reduced_cnf_path = "/tmp/" + std::to_string(getpid()) + std::to_string(tmp_gen()) + "_reduced.cnf";
        std::string cmd = "./bin/coprocessor -enabled_cp3 -up -subsimp -no-bve -no-bce -no-dense -dimacs=" + 
            reduced_cnf_path + " " + params.input_cnf_path;
        int _ = system(cmd.c_str());
        real_input_cnf_path = reduced_cnf_path;
    }

    CNFInfo cnf(real_input_cnf_path);
    cnf.dump(nvar, nclauses, clauses, pos_in_cls, neg_in_cls);

    combnum = new long long *[twise + 1];
    for (int i = 2; i <= twise; ++i){
        combnum[i] = new long long[nvar + 1]{0};
    }
    for (int i = 2; i <= nvar; ++i){
        combnum[2][i] = 1ll * i * (i - 1) / 2;
        for (int j = 3; j <= twise; ++j){
            combnum[j][i] = combnum[j][i - 1] + combnum[j - 1][i - 1];
        }
    }

    tuple_iter_base = new int[twise + 1];
    get_gain_for_flip_tuple_tmp = vector<pair<int, int> >(twise, {0, 0});
    get_gain_for_forcetuple_tuple_tmp = vector<pair<int, int> >(twise, {0, 0});

    ca.clear();
    last_ca.clear();
    FILE *in_ca = fopen(params.init_ca_path.c_str(), "r");
    vector<int> tmp(nvar, 0);
    for (int c, p = 0; (c = fgetc(in_ca)) != EOF; ){
        if (c == '\n'){
            ca.emplace_back(tmp);
            p = 0;
        } else if (isdigit(c)){
            tmp[p++] = c - '0';
        }
    }
    fclose(in_ca);

    std::cout << "c read PCA from " << params.init_ca_path << "; initial size = " << ca.size() << std::endl;
    last_ca = ca;

    testcase_size = ca.size();
    num_combination_all_possible = combnum[twise][nvar];
    num_tuple_all_possible = num_combination_all_possible * (1ll << twise);
    tuple_cov_cnt = vector<int>(num_tuple_all_possible, 0);
    uncovered_tuples.clear();
    uncovered_tuple_positions = vector<int>(num_tuple_all_possible, -1);
    cdcl_solver = new Minisat::Solver;

    Minisat::vec<Minisat::Lit> lits;
    for (const vector<int>& cl: clauses){
        lits.clear();
        for (int rawvar: cl){
            int var = abs(rawvar) - 1;
            while (var >= cdcl_solver->nVars()) cdcl_solver->newVar();
            lits.push(Minisat::mkLit(var, rawvar < 0));
        }
        cdcl_solver->addClause_(lits);
    }

    while (nvar > cdcl_solver->nVars()){
        cdcl_solver->newVar();
    }

    for (const vector<int>& testcase: ca){
        target_twise = twise;
        tuple_iter_action_mode = 1;
        iterate_all_tuples(0, 0, testcase);
    }

    for (int i = 0; i < num_tuple_all_possible; ++i){
        if (tuple_cov_cnt[i] == 0) tuple_cov_cnt[i] = -1;
    }

    nclauses = clauses.size();
    clauses_cov = vector<vector<int> >(testcase_size, vector<int>(nclauses, 0));
    for (int i = 0; i < testcase_size; ++i){
        const vector<int>& testcase = ca[i];
        vector<int>& cur_clauses_cov = clauses_cov[i];
        for (int j = 0; j < nvar; ++j){
            const vector<int>& vec = (testcase[j] ? pos_in_cls[j + 1]: neg_in_cls[j + 1]);
            for (int x: vec) ++cur_clauses_cov[x];
        }        
    }

    if (params.flag_use_cell_tabu) last_flipped_time = vector<vector<int> >(testcase_size, vector<int>(nvar, -testcase_taboo - 1));

    if (!params.flag_use_testcase_taboo) testcase_taboo = -1;

    greedy_limit = 0;
    last_greedy_time = vector<int>(testcase_size, -testcase_taboo - 1);

    tmp_flip_sequence = vector<int>(testcase_size * nvar, 0);
    std::iota(tmp_flip_sequence.begin(), tmp_flip_sequence.end(), 0);

    remove_testcase_strategy = &LSOptimizer::remove_testcase_greedily;

    if (testcase_size <= testcase_taboo + 1) testcase_taboo = std::min(testcase_size - 2, testcase_size >> 1);

    __last_success_step = 0;
}

LSOptimizer::~LSOptimizer(){
    delete cdcl_solver;

    for (int i = 2; i <= twise; ++i){
        delete [] combnum[i];
    }
    delete [] combnum;

    delete [] tuple_iter_base;
}

int LSOptimizer::remove_testcase_greedily(){
    int besttc = -1;
    long long mini = num_combination_all_possible + 1;
    for (int i = 0; i < testcase_size; ++i){
        long long res = new_uncovered_tuples_after_remove_testcase(ca[i], false);
        if (res < mini){
            mini = res;
            besttc = i;
        }
    }
    return besttc;
}

int LSOptimizer::remove_testcase_randomly(){
    return gen() % testcase_size;
}

void LSOptimizer::search(){
    int cur_step = 0;
    for (; testcase_size > 1 && __last_success_step < __stop_length; ){
        if (uncovered_tuples.empty()){
            std::cout << "\033[;32mc current PCA size: " << testcase_size << ", step #" << cur_step << " \033[0m" << std::endl;
            last_ca = ca;

            int rmid = (this->*remove_testcase_strategy)();
            new_uncovered_tuples_after_remove_testcase(ca[rmid], true);

            if (verbosity >= 2){
                std::cout << "c     testcase removed: #" << rmid << std::endl;
            }

            int newsz = 0;
            for (const vector<pair<int, int> >& p: tmp_break){
                uncovered_tuples.push_back(p);
                uncovered_tuple_positions[getXtupleid(p)] = newsz;
                ++newsz;
            }

            if (rmid != testcase_size - 1){
                std::swap(clauses_cov[rmid], clauses_cov[testcase_size - 1]);
                std::swap(ca[rmid], ca[testcase_size - 1]);
                if (__use_cell_tabu) std::swap(last_flipped_time[rmid], last_flipped_time[testcase_size - 1]);
                std::swap(last_greedy_time[rmid], last_greedy_time[testcase_size - 1]);
            }

            for (int j = 0; j < nvar; ++j) tmp_flip_sequence.pop_back();
            std::iota(tmp_flip_sequence.begin(), tmp_flip_sequence.end(), 0);
            clauses_cov.pop_back();
            ca.pop_back();
            if (__use_cell_tabu) last_flipped_time.pop_back();
            last_greedy_time.pop_back();
            --testcase_size;
            if (testcase_size <= testcase_taboo + 1) testcase_taboo = std::min(testcase_size - 2, testcase_size >> 1);

            __last_success_step = 0;

            continue;
        }

        ++cur_step;
        ++__last_success_step;

        int cyc = gen() % 100;
        if (cyc < 1 || !random_greedy_step()){
            random_step();
        }

        if (verbosity >= 1) std::cout << "c     #uncovered tuples: " << uncovered_tuples.size() << std::endl;
    }

    final_validate();
}

void LSOptimizer::flip_bit(int tcid, int vid){
    auto res = get_gain_for_flip(tcid, vid, true);
    if (verbosity >= 1) std::cout << "c   cnt_break = " << res.second.first << ", cnt_gain = " << res.second.second << std::endl;
    
    break_and_gain();

    int curbit = ca[tcid][vid];

    update_clause_cov_by_flip(tcid, vid, curbit, false);
    ca[tcid][vid] ^= 1;
}

bool LSOptimizer::random_greedy_step(){
    int uncovered_cnt = uncovered_tuples.size();
    int picked_tuple = gen() % uncovered_cnt;

    vector<pair<int, int> > tp = uncovered_tuples[picked_tuple];

    int besttcid = -1;
    long long maxi = -num_combination_all_possible - 1; 
    for (int i = 0; i < testcase_size; ++i){
        if (!__use_cell_tabu){
            if (greedy_limit - last_greedy_time[i] <= testcase_taboo){
                continue;
            }
        } else {
        }

        auto res = get_gain_for_forcetuple(i, tp, false);
        if (res.first){
            int net_gain = res.second.second - res.second.first;
            if (net_gain > maxi){
                besttcid = i;
                maxi = net_gain;
            }
        }
    }

    if (besttcid != -1){
        forcetuple(besttcid, tp);

        ++greedy_limit;
        if (!__use_cell_tabu){
            last_greedy_time[besttcid] = greedy_limit;
        }

        return true;
    }

    if (verbosity >= 1) std::cout << "c    random greedy step failed!" << std::endl;

    if ((gen() % 100) < __forced_greedy_percent){
        greedy_step_forced(tp);
        return true;
    }

    return false;
}

void LSOptimizer::greedy_step_forced(const vector<pair<int, int> >& chosen_tp){
    Minisat::vec<Minisat::Lit> assu;

    for (auto& pp: chosen_tp){
        assu.push(Minisat::mkLit(pp.first, pp.second == 0));
    }

    int besttcid = -1;
    long long maxi = -num_combination_all_possible - 1; 
    vector<int> besttc2;

    for (int i = 0; i < testcase_size; ++i){
        if (!__use_cell_tabu){
            if (greedy_limit - last_greedy_time[i] <= testcase_taboo){
                continue;
            }
        } else {
        }

        for (int j = 0; j < nvar; ++j){
            cdcl_solver->setPolarity(j, ca[i][j] == 0);
        }

        if (!cdcl_solver->simplify() || !cdcl_solver->solve(assu)){
            std::cout << "c \033[1;31mError: SAT solve failing!\033[0m" << std::endl;
            return;
        }

        vector<int> tc2;
        for (int j = 0; j < nvar; ++j){
            if (cdcl_solver->model[j] == (Minisat::lbool((uint8_t)0))) tc2.emplace_back(1);
            else tc2.emplace_back(0);
        }
        auto res = get_gain_for_forcetestcase(i, tc2, false);
        int net_gain = res.second - res.first;
        if (net_gain > maxi){
            besttcid = i;
            besttc2 = tc2;
            maxi = net_gain;
        }
    }

    forcetestcase(besttcid, besttc2);

    ++greedy_limit;
    if (!__use_cell_tabu){
        last_greedy_time[besttcid] = greedy_limit;
    }

}

void LSOptimizer::random_step(){
    if (!__enable_random_step) return ;

    std::shuffle(tmp_flip_sequence.begin(), tmp_flip_sequence.end(), gen);
    for (int pp: tmp_flip_sequence){
        int tcid = pp / nvar, vid = pp % nvar;
        auto res = get_gain_for_flip(tcid, vid, false);
        if (res.first){
            flip_bit(tcid, vid);
            if (verbosity >= 1) std::cout << "c    ramdom step for cell (" << tcid << ", " << vid << ")!" << std::endl;
            break;
        }
    }
}

long long LSOptimizer::new_uncovered_tuples_after_remove_testcase(const vector<int>& tc, bool dump){
    if (dump) tmp_break.clear();

    target_twise = twise;
    tuple_iter_action_mode = 2;
    new_uncovered_tuples_after_remove_testcase_dump = dump;
    new_uncovered_tuples_after_remove_testcase_break_cnt = 0;
    iterate_all_tuples(0, 0, tc);

    return new_uncovered_tuples_after_remove_testcase_break_cnt;
}

vector<vector<int> > LSOptimizer::get_testcase_set(){
    return last_ca;
}

pair<bool, pair<int, int> > LSOptimizer::get_gain_for_flip(int tcid, int vid, bool dump){
    const vector<int>& tc = ca[tcid];
    int curbit = tc[vid];
    
    const vector<int>& var_cov_old = (curbit ? pos_in_cls[vid + 1]: neg_in_cls[vid + 1]);
    const vector<int>& var_cov_new = (curbit ? neg_in_cls[vid + 1]: pos_in_cls[vid + 1]);
    vector<int>& cur_clauses_cov = clauses_cov[tcid];
    
    bool has0 = false;
    for (int cid: var_cov_new) ++cur_clauses_cov[cid];
    for (int cid: var_cov_old){
        --cur_clauses_cov[cid];
        if (cur_clauses_cov[cid] == 0){
            has0 = true;
        }
    }

    for (int cid: var_cov_new) --cur_clauses_cov[cid];
    for (int cid: var_cov_old) ++cur_clauses_cov[cid];

    if (has0) return {false, {0, 0}};

    if (dump){
        tmp_break.clear();
        tmp_gain.clear();
    }

    target_twise = twise - 1;
    tuple_iter_action_mode = 4;
    get_gain_for_flip_dump = dump;
    get_gain_for_flip_vid = vid;
    get_gain_for_flip_curbit = curbit;
    get_gain_for_flip_break_cnt = 0;
    get_gain_for_flip_gain_cnt = 0;
    iterate_all_tuples(0, 0, tc);

    return {true, {get_gain_for_flip_break_cnt, get_gain_for_flip_gain_cnt}};
}

pair<bool, pair<int, int> > LSOptimizer::get_gain_for_forcetuple(int tcid, const vector<pair<int, int> >& chosen_tp, bool dump){
    const vector<int>& tc = ca[tcid];
    
    vector<int>& cur_clauses_cov = clauses_cov[tcid];

    for (auto& pp: chosen_tp){
        int vid = pp.first;
        int curbit = tc[vid], tt = pp.second;
        if (curbit == tt) continue;

        const vector<int>& var_cov_old = (curbit ? pos_in_cls[vid + 1]: neg_in_cls[vid + 1]);
        const vector<int>& var_cov_new = (tt ? pos_in_cls[vid + 1]: neg_in_cls[vid + 1]);
        for (int cid: var_cov_new) ++cur_clauses_cov[cid];   
        for (int cid: var_cov_old) --cur_clauses_cov[cid]; 
    }

    bool has0 = false;
    for (int i = 0; i < nclauses; ++i){
        if (cur_clauses_cov[i] == 0){
            has0 = true;
            break;
        }
    }

    for (auto& pp: chosen_tp){
        int vid = pp.first;
        int curbit = tc[vid], tt = pp.second;
        if (curbit == tt) continue;

        const vector<int>& var_cov_old = (curbit ? pos_in_cls[vid + 1]: neg_in_cls[vid + 1]);
        const vector<int>& var_cov_new = (tt ? pos_in_cls[vid + 1]: neg_in_cls[vid + 1]);
        for (int cid: var_cov_new) --cur_clauses_cov[cid];   
        for (int cid: var_cov_old) ++cur_clauses_cov[cid]; 
    }

    if (has0) return {false, {0, 0}};

    if (dump){
        tmp_break.clear();
        tmp_gain.clear();
    }

    tuple_iter_action_mode = 5;
    get_gain_for_forcetuple_dump = dump;
    get_gain_for_forcetuple_chosen_tp = chosen_tp;
    get_gain_for_forcetuple_break_cnt = 0;
    get_gain_for_forcetuple_gain_cnt = 0;

    for (target_twise = 1; target_twise < twise; ++target_twise){
        iterate_all_tuples(0, 0, tc);
    }
    for (int i = 0; i < twise; ++i){
        get_gain_for_forcetuple_tuple_tmp[i].first = chosen_tp[i].first;
        get_gain_for_forcetuple_tuple_tmp[i].second = tc[chosen_tp[i].first];
    }
    long long tupleid_oldtp = getXtupleid(get_gain_for_forcetuple_tuple_tmp);
    long long tupleid_newtp = getXtupleid(chosen_tp);

    if (dump){
        int& tcc_old = tuple_cov_cnt[tupleid_oldtp];
        --tcc_old;
        if (tcc_old == 0){
            tmp_break.push_back(get_gain_for_forcetuple_tuple_tmp);
        }
        int& tcc_new = tuple_cov_cnt[tupleid_newtp];
        ++tcc_new;
        if (tcc_new == 1){
            tmp_gain.push_back(chosen_tp);
        }
    } else {
        if (tuple_cov_cnt[tupleid_oldtp] == 1){
            ++get_gain_for_forcetuple_break_cnt;
        }
        if (tuple_cov_cnt[tupleid_newtp] == 0){
            ++get_gain_for_forcetuple_gain_cnt;
        }
    }

    return {true, {get_gain_for_forcetuple_break_cnt, get_gain_for_forcetuple_gain_cnt}};
}

void LSOptimizer::forcetuple(int tcid, const vector<pair<int, int> >& chosen_tp){
    auto res = get_gain_for_forcetuple(tcid, chosen_tp, true);
    if (verbosity >= 1) std::cout << "c   cnt_break = " << res.second.first << ", cnt_gain = " << res.second.second << std::endl;
    
    break_and_gain();

    const vector<int>& tc = ca[tcid];

    int diffcnt = 0;
    for (auto& pp: chosen_tp){
        if (pp.second != tc[pp.first]){
            ++diffcnt;
        }
    }

    for (auto& pp: chosen_tp){
        if (pp.second != tc[pp.first]){
            --diffcnt;
            update_clause_cov_by_flip(tcid, pp.first, tc[pp.first], diffcnt != 0);
        }
    }

    for (auto& pp: chosen_tp){
        ca[tcid][pp.first] = pp.second;
    }
}

void LSOptimizer::break_and_gain(){
    int uncovered_cnt = uncovered_tuples.size();
    for (const vector<pair<int, int> >& p: tmp_break){
        uncovered_tuples.push_back(p);
        uncovered_tuple_positions[getXtupleid(p)] = uncovered_cnt;
        ++uncovered_cnt;
    }
    for (const vector<pair<int, int> >& p: tmp_gain){
        long long this_tuple_id = getXtupleid(p);
        int& original_pos = uncovered_tuple_positions[this_tuple_id];
        if (original_pos != uncovered_cnt - 1){
            long long tail_tupleid = getXtupleid(uncovered_tuples[uncovered_cnt - 1]);
            uncovered_tuple_positions[tail_tupleid] = original_pos;
            swap(uncovered_tuples[original_pos], uncovered_tuples[uncovered_cnt - 1]);
        }
        uncovered_tuples.pop_back();
        original_pos = -1;
        --uncovered_cnt;
    }
}

void LSOptimizer::update_clause_cov_by_flip(int tcid, int vid, int curbit, bool trailing){
    const vector<int>& var_cov_old = (curbit ? pos_in_cls[vid + 1]: neg_in_cls[vid + 1]);
    const vector<int>& var_cov_new = (curbit ? neg_in_cls[vid + 1]: pos_in_cls[vid + 1]);
    vector<int>& cur_clauses_cov = clauses_cov[tcid];
    for (int cid: var_cov_new) ++cur_clauses_cov[cid];
    for (int cid: var_cov_old) --cur_clauses_cov[cid];  
}

bool LSOptimizer::validate_testcase(const vector<int>& tc){
    vector<int> cl(nclauses, 0);
    for (int j = 0; j < nvar; ++j){
        const vector<int>& var_cov = (tc[j] ? pos_in_cls[j + 1]: neg_in_cls[j + 1]);
        for (int cid: var_cov) ++cl[cid];
    }
    for (int j = 0; j < nclauses; ++j){
        if (cl[j] == 0) return false;
    }
    return true;
}

void LSOptimizer::validate(){
    for (int i = 0; i < testcase_size; ++i){
        vector<int> cl(nclauses, 0);
        for (int j = 0; j < nvar; ++j){
            const vector<int>& var_cov = (ca[i][j] ? pos_in_cls[j + 1]: neg_in_cls[j + 1]);
            for (int cid: var_cov) ++cl[cid];
        }
        for (int j = 0; j < nclauses; ++j){
            if (cl[j] == 0){
                std::cout << "c \033[1;31mError: SAT broken in validation!\033[0m" << std::endl;
                return ;
            } else if (cl[j] != clauses_cov[i][j]){
                std::cout << "c \033[1;31mError: clauses_cov mismatches at testcase #" << i << "!\033[0m" << std::endl;
                return ;
            }
        }
    }
}

void LSOptimizer::final_validate(){
    for (const vector<int>& testcase: ca){
        target_twise = twise;
        tuple_iter_action_mode = 0;
        iterate_all_tuples(0, 0, testcase);
    }

    long long tuple2_cnt = 0;
    for (int i = 0; i < num_tuple_all_possible; ++i){
        if (tuple_cov_cnt[i] == 0){
            ++tuple2_cnt;
        } else if (tuple_cov_cnt[i] != -1){
            std::cout << "c \033[1;31mError: tuple_cov not accurate!\033[0m" << std::endl;
            return ;
        }
    }

    std::cout << "c total tuple num = " << tuple2_cnt << std::endl;
}

pair<int, int> LSOptimizer::get_gain_for_forcetestcase(int tcid, const vector<int>& tc2, bool dump){

    const vector<int>& tc = ca[tcid];

    if (dump){
        tmp_break.clear();
        tmp_gain.clear();
    }

    target_twise = twise;
    tuple_iter_action_mode = 3;
    get_gain_for_forcetestcase_dump = dump;
    get_gain_for_forcetestcase_break_cnt = 0;
    get_gain_for_forcetestcase_gain_cnt = 0;
    iterate_all_tuples_double(0, 0, tc, 0, tc2);

    return {get_gain_for_forcetestcase_break_cnt, get_gain_for_forcetestcase_gain_cnt};
}

void LSOptimizer::forcetestcase(int tcid, const vector<int>& tc2){
    auto res = get_gain_for_forcetestcase(tcid, tc2, true);
    if (verbosity >= 1) std::cout << "c   cnt_break = " << res.first << ", cnt_gain = " << res.second << std::endl;
    
    break_and_gain();

    vector<int>& cur_clauses_cov = clauses_cov[tcid];
    for (int i = 0; i < nclauses; ++i) cur_clauses_cov[i] = 0;
    for (int i = 0; i < nvar; ++i){
        const vector<int>& var_cov_new = (tc2[i] ? pos_in_cls[i + 1]: neg_in_cls[i + 1]);
        for (int x: var_cov_new) ++cur_clauses_cov[x];
    }

    vector<int>& tc = ca[tcid];

    tc = tc2;
}

void LSOptimizer::iterate_all_tuples(int cur, long long tupleid, const vector<int>& tc){
    if (cur == target_twise){
        long long tuplebase = getXbase();
        
        switch (tuple_iter_action_mode){
        case 0: 
            tupleid += tuplebase;
            tuple_iter_local_action_0(tupleid);
            break;
        case 1: 
            tupleid += tuplebase;
            tuple_iter_local_action_1(tupleid);
            break;
        case 2: 
            tupleid += tuplebase;
            tuple_iter_local_action_2(tupleid, tc);
            break;
        case 4: 
            tuple_iter_local_action_4(tc);
            break;
        case 5: 
            tuple_iter_local_action_5(tc);
            break;
        }
        return ;
    }

    int lst = -1; 
    if (cur > 0) lst = tuple_iter_base[cur - 1];
    int ubd = nvar - target_twise + cur;
    for (int i = lst + 1; i <= ubd; ++i){
        tuple_iter_base[cur] = i;
        iterate_all_tuples(cur + 1, (tupleid << 1) + (tc[i] ? num_combination_all_possible: 0), tc);
    }
}

void LSOptimizer::iterate_all_tuples_double(int cur, long long tupleid_old, const vector<int>& tc_old, long long tupleid_new, const vector<int>& tc_new){
    if (cur == target_twise){
        long long tuplebase = getXbase();
        tupleid_old += tuplebase;
        tupleid_new += tuplebase;
        
        switch (tuple_iter_action_mode){
        case 3: 
            tuple_iter_local_action_3(tupleid_old, tc_old, tupleid_new, tc_new);
            break;
        }
        return ;
    }

    int lst = -1; 
    if (cur > 0) lst = tuple_iter_base[cur - 1];
    int ubd = nvar - target_twise + cur;
    for (int i = lst + 1; i <= ubd; ++i){
        tuple_iter_base[cur] = i;
        iterate_all_tuples_double(cur + 1, (tupleid_old << 1) + (tc_old[i] ? num_combination_all_possible: 0), tc_old, 
            (tupleid_new << 1) + (tc_new[i] ? num_combination_all_possible: 0), tc_new);
    }
}