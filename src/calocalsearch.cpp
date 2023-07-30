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
    num_combination_all_possible = 1ll * nvar * (nvar - 1) / 2ll;
    num_tuple_all_possible = 4ll * num_combination_all_possible;
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
        long long tuple_pos = 0;
        for (int j = 0; j < nvar; ++j){
            long long tuple_id_base0 = testcase[j] ? num_combination_all_possible << 1: 0;
            for (int k = j + 1; k < nvar; ++k, ++tuple_pos){
                long long tuple_id_base1 = testcase[k] ? num_combination_all_possible: 0;
                long long tuple_id = tuple_id_base0 + tuple_id_base1 + tuple_pos;
                ++tuple_cov_cnt[tuple_id];
            }
        }
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
}

long long LSOptimizer::get2tupleid(int i, int vi, int j, int vj){
    long long base = (vi << 1 | vj) * num_combination_all_possible;
    long long pos = (2ll * nvar - i - 1) * i / 2 + j - i - 1;
    return base + pos;
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

            int break_cnt = 0;
            long long tuple_pos = 0;
            for (int j = 0; j < nvar; ++j){
                long long tuple_id_base0 = ca[rmid][j] ? num_combination_all_possible << 1: 0;
                for (int k = j + 1; k < nvar; ++k, ++tuple_pos){
                    long long tuple_id_base1 = ca[rmid][k] ? num_combination_all_possible: 0;
                    long long tuple_id = tuple_id_base0 + tuple_id_base1 + tuple_pos;
                    --tuple_cov_cnt[tuple_id];
                    if (tuple_cov_cnt[tuple_id] == 0){
                        ++break_cnt;
                    }
                }
            }

            int newsz = 0;
            for (const pair<int, int>& p: tmp_break){
                uncovered_tuples.push_back(p);
                int comp1 = p.first, comp2 = p.second;
                uncovered_tuple_positions[get2tupleid(abs(comp1) - 1, comp1 > 0 ? 1: 0, abs(comp2) - 1, comp2 > 0 ? 1: 0)] = newsz;
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

        if (verbosity >= 3){
            for (const pair<int, int>& p: uncovered_tuples){
                std::cout << "c        uncovered is: (" << p.first << ", " << p.second << ")" << std::endl;
            }
        }

        int cyc = gen() % 100;
        if (cyc < 1 || !random_greedy_step()){
            random_step();
        }

        if (verbosity >= 1) std::cout << "c     #uncovered tuples: " << uncovered_tuples.size() << std::endl;

        validate();
    }

    final_validate();
}

void LSOptimizer::flip_bit(int tcid, int vid){
    auto res = get_gain_for_flip(tcid, vid, true);
    if (verbosity >= 1) std::cout << "c   cnt_break = " << res.second.first << ", cnt_gain = " << res.second.second << std::endl;
    
    break_and_gain();

    int curbit = ca[tcid][vid];

    update_clause_cov_by_flip(tcid, vid, curbit, false);

    for (int i = 0; i < nvar; ++i){
        if (i == vid) continue;
        int thatbit = ca[tcid][i];
        long long tupleid_old = (i < vid ? get2tupleid(i, thatbit, vid, curbit): get2tupleid(vid, curbit, i, thatbit));
        long long tupleid_new = (i < vid ? get2tupleid(i, thatbit, vid, curbit ^ 1): get2tupleid(vid, curbit ^ 1, i, thatbit));
        --tuple_cov_cnt[tupleid_old];
        ++tuple_cov_cnt[tupleid_new];
    }

    ca[tcid][vid] ^= 1;
}

bool LSOptimizer::random_greedy_step(){
    int uncovered_cnt = uncovered_tuples.size();
    int picked_tuple = gen() % uncovered_cnt;
    int comp1 = uncovered_tuples[picked_tuple].first, comp2 = uncovered_tuples[picked_tuple].second;

    int besttcid = -1;
    long long maxi = -num_combination_all_possible - 1; 
    for (int i = 0; i < testcase_size; ++i){
        if (!__use_cell_tabu){
            if (greedy_limit - last_greedy_time[i] <= testcase_taboo){
                continue;
            }
        } else {
            int vid1 = abs(comp1) - 1, vid2 = abs(comp2) - 1;
            bool check1 = ca[i][vid1] == (comp1 > 0 ? 1: 0);
            bool check2 = ca[i][vid2] == (comp2 > 0 ? 1: 0);
            if (greedy_limit - std::max(check1 ? (-testcase_taboo - 1): last_flipped_time[i][vid1], 
                check2 ? (-testcase_taboo - 1): last_flipped_time[i][vid2]) <= testcase_taboo){
                continue;
            }
        }

        auto res = get_gain_for_forcetuple(i, comp1, comp2, false);
        if (res.first){
            int net_gain = res.second.second - res.second.first;
            if (net_gain > maxi){
                besttcid = i;
                maxi = net_gain;
            }
        }
    }

    if (besttcid != -1){
        forcetuple(besttcid, comp1, comp2);

        ++greedy_limit;
        if (!__use_cell_tabu){
            last_greedy_time[besttcid] = greedy_limit;
        }

        if (verbosity >= 1) std::cout << "c    random greedy step for (" << comp1 << ", " << comp2 << ") and testcase #" << besttcid << "!" << std::endl;
        return true;
    }

    if (verbosity >= 1) std::cout << "c    random greedy step failed!" << std::endl;

    if ((gen() % 100) < __forced_greedy_percent){
        greedy_step_forced(comp1, comp2);
        return true;
    }

    return false;
}

void LSOptimizer::greedy_step_forced(int comp1, int comp2){
    Minisat::vec<Minisat::Lit> assu;

    assu.push(Minisat::mkLit(abs(comp1) - 1, comp1 < 0));
    assu.push(Minisat::mkLit(abs(comp2) - 1, comp2 < 0));

    int besttcid = -1;
    long long maxi = -num_combination_all_possible - 1; 
    vector<int> besttc2;

    for (int i = 0; i < testcase_size; ++i){
        if (!__use_cell_tabu){
            if (greedy_limit - last_greedy_time[i] <= testcase_taboo){
                continue;
            }
        } else {
            int vid1 = abs(comp1) - 1, vid2 = abs(comp2) - 1;
            bool check1 = ca[i][vid1] == (comp1 > 0 ? 1: 0);
            bool check2 = ca[i][vid2] == (comp2 > 0 ? 1: 0);
            if (greedy_limit - std::max(check1 ? (-testcase_taboo - 1): last_flipped_time[i][vid1], 
                check2 ? (-testcase_taboo - 1): last_flipped_time[i][vid2]) <= testcase_taboo){
                continue;
            }
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

    if (verbosity >= 1) std::cout << "c    forced greedy step for (" << comp1 << ", " << comp2 << ") and testcase #" << besttcid << "!" << std::endl;

    if (verbosity >= 2){
        std::cout << "c    forced testcase: ";
        for (int x: besttc2){
            std::cout << x;
        }
        std::cout << std::endl;
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
    int break_cnt = 0;
    if (dump) tmp_break.clear();

    long long tuple_pos = 0;
    for (int j = 0; j < nvar; ++j){
        long long tuple_id_base0 = tc[j] ? num_combination_all_possible << 1: 0;
        for (int k = j + 1; k < nvar; ++k, ++tuple_pos){
            long long tuple_id_base1 = tc[k] ? num_combination_all_possible: 0;
            long long tuple_id = tuple_id_base0 + tuple_id_base1 + tuple_pos;
            if (tuple_cov_cnt[tuple_id] == 1){
                ++break_cnt;
                if (dump) tmp_break.emplace_back((tc[j] ? (j+1): (-j-1)), (tc[k] ? (k+1): (-k-1)));
            }
        }
    }

    return break_cnt;
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

    int break_cnt = 0, gain_cnt = 0;

    if (dump){
        tmp_break.clear();
        tmp_gain.clear();
    }
    
    for (int i = 0; i < nvar; ++i){
        if (i == vid) continue;
        int thatbit = tc[i];
        long long tupleid = (i < vid ? get2tupleid(i, thatbit, vid, curbit): get2tupleid(vid, curbit, i, thatbit));
        if (tuple_cov_cnt[tupleid] == 1){
            ++break_cnt;
            if (dump){
                if (i < vid) tmp_break.emplace_back((thatbit ? (i+1): (-i-1)), (curbit ? (vid+1): (-vid-1)));
                else tmp_break.emplace_back((curbit ? (vid+1): (-vid-1)), (thatbit ? (i+1): (-i-1)));
            }
        }
    }

    curbit ^= 1;
    for (int i = 0; i < nvar; ++i){
        if (i == vid) continue;
        int thatbit = tc[i];
        long long tupleid = (i < vid ? get2tupleid(i, thatbit, vid, curbit): get2tupleid(vid, curbit, i, thatbit));
        if (tuple_cov_cnt[tupleid] == 0){
            ++gain_cnt;
            if (dump){
                if (i < vid) tmp_gain.emplace_back((thatbit ? (i+1): (-i-1)), (curbit ? (vid+1): (-vid-1)));
                else tmp_gain.emplace_back((curbit ? (vid+1): (-vid-1)), (thatbit ? (i+1): (-i-1)));
            }
        }
    }

    return {true, {break_cnt, gain_cnt}};
}

pair<bool, pair<int, int> > LSOptimizer::get_gain_for_forcetuple(int tcid, int comp1, int comp2, bool dump){
    int vid1 = abs(comp1) - 1, vid2 = abs(comp2) - 1;
    int t1 = (comp1 > 0 ? 1: 0), t2 = (comp2 > 0 ? 1: 0);

    const vector<int>& tc = ca[tcid];
    int curbit1 = tc[vid1], curbit2 = tc[vid2];
    int comp1_old = (curbit1 ? (vid1+1): (-vid1-1)), comp2_old = (curbit2 ? (vid2+1): (-vid2-1));

    if (curbit1 == t1 && curbit2 == t2){
        std::cout << "c \033[1;31mError: this \"uncovered tuple\" is covered!\033[0m" << std::endl;
        return {false, {0, 0}};
    }
    
    const vector<int>& var_cov_old1 = (curbit1 ? pos_in_cls[vid1 + 1]: neg_in_cls[vid1 + 1]);
    const vector<int>& var_cov_new1 = (t1 ? pos_in_cls[vid1 + 1]: neg_in_cls[vid1 + 1]);
    const vector<int>& var_cov_old2 = (curbit2 ? pos_in_cls[vid2 + 1]: neg_in_cls[vid2 + 1]);
    const vector<int>& var_cov_new2 = (t2 ? pos_in_cls[vid2 + 1]: neg_in_cls[vid2 + 1]);

    vector<int>& cur_clauses_cov = clauses_cov[tcid];
    
    if (t1 != curbit1){
        for (int cid: var_cov_new1) ++cur_clauses_cov[cid];   
        for (int cid: var_cov_old1) --cur_clauses_cov[cid]; 
    }
    if (t2 != curbit2){
        for (int cid: var_cov_new2) ++cur_clauses_cov[cid];   
        for (int cid: var_cov_old2) --cur_clauses_cov[cid]; 
    }

    bool has0 = false;
    for (int i = 0; i < nclauses; ++i){
        if (cur_clauses_cov[i] == 0){
            has0 = true;
            break;
        }
    }

    if (t1 != curbit1){
        for (int cid: var_cov_new1) --cur_clauses_cov[cid];   
        for (int cid: var_cov_old1) ++cur_clauses_cov[cid]; 
    }
    if (t2 != curbit2){
        for (int cid: var_cov_new2) --cur_clauses_cov[cid];   
        for (int cid: var_cov_old2) ++cur_clauses_cov[cid]; 
    }

    if (has0) return {false, {0, 0}};

    int break_cnt = 0, gain_cnt = 0;

    if (dump){
        tmp_break.clear();
        tmp_gain.clear();
    }
    
    if (t1 != curbit1){
        for (int i = 0; i < nvar; ++i){
            if (i == vid1 || i == vid2) continue;
            int thatbit = tc[i];
            long long tupleid = (i < vid1 ? get2tupleid(i, thatbit, vid1, curbit1): get2tupleid(vid1, curbit1, i, thatbit));
            if (tuple_cov_cnt[tupleid] == 1){
                ++break_cnt;
                if (dump){
                    if (i < vid1) tmp_break.emplace_back((thatbit ? (i+1): (-i-1)), comp1_old);
                    else tmp_break.emplace_back(comp1_old, (thatbit ? (i+1): (-i-1)));
                }
            }
        }

        for (int i = 0; i < nvar; ++i){
            if (i == vid1 || i == vid2) continue;
            int thatbit = tc[i];
            long long tupleid = (i < vid1 ? get2tupleid(i, thatbit, vid1, t1): get2tupleid(vid1, t1, i, thatbit));
            if (tuple_cov_cnt[tupleid] == 0){
                ++gain_cnt;
                if (dump){
                    if (i < vid1) tmp_gain.emplace_back((thatbit ? (i+1): (-i-1)), comp1);
                    else tmp_gain.emplace_back(comp1, (thatbit ? (i+1): (-i-1)));
                }
            }
        }
    }

    if (t2 != curbit2){
        for (int i = 0; i < nvar; ++i){
            if (i == vid1 || i == vid2) continue;
            int thatbit = tc[i];
            long long tupleid = (i < vid2 ? get2tupleid(i, thatbit, vid2, curbit2): get2tupleid(vid2, curbit2, i, thatbit));
            if (tuple_cov_cnt[tupleid] == 1){
                ++break_cnt;
                if (dump){
                    if (i < vid2) tmp_break.emplace_back((thatbit ? (i+1): (-i-1)), comp2_old);
                    else tmp_break.emplace_back(comp2_old, (thatbit ? (i+1): (-i-1)));
                }
            }
        }

        for (int i = 0; i < nvar; ++i){
            if (i == vid1 || i == vid2) continue;
            int thatbit = tc[i];
            long long tupleid = (i < vid2 ? get2tupleid(i, thatbit, vid2, t2): get2tupleid(vid2, t2, i, thatbit));
            if (tuple_cov_cnt[tupleid] == 0){
                ++gain_cnt;
                if (dump){
                    if (i < vid2) tmp_gain.emplace_back((thatbit ? (i+1): (-i-1)), comp2);
                    else tmp_gain.emplace_back(comp2, (thatbit ? (i+1): (-i-1)));
                }
            }
        }
    }

    long long tupleid_oldtp = get2tupleid(vid1, curbit1, vid2, curbit2);
    if (tuple_cov_cnt[tupleid_oldtp] == 1){
        ++break_cnt;
        if (dump) tmp_break.emplace_back(comp1_old, comp2_old);
    }
    long long tupleid_newtp = get2tupleid(vid1, t1, vid2, t2);
    ++gain_cnt;
    if (dump) tmp_gain.emplace_back(comp1, comp2);

    return {true, {break_cnt, gain_cnt}};
}

void LSOptimizer::forcetuple(int tcid, int comp1, int comp2){
    auto res = get_gain_for_forcetuple(tcid, comp1, comp2, true);
    if (verbosity >= 1) std::cout << "c   cnt_break = " << res.second.first << ", cnt_gain = " << res.second.second << std::endl;
    
    break_and_gain();

    int vid1 = abs(comp1) - 1, vid2 = abs(comp2) - 1;
    int t1 = (comp1 > 0 ? 1: 0), t2 = (comp2 > 0 ? 1: 0);

    const vector<int>& tc = ca[tcid];
    int curbit1 = tc[vid1], curbit2 = tc[vid2];

    if (t1 != curbit1) update_clause_cov_by_flip(tcid, vid1, curbit1, t2 != curbit2);
    if (t2 != curbit2) update_clause_cov_by_flip(tcid, vid2, curbit2, false);

    if (t1 != curbit1){
        for (int i = 0; i < nvar; ++i){
            if (i == vid1 || i == vid2) continue;
            int thatbit = tc[i];
            long long tupleid_old = (i < vid1 ? get2tupleid(i, thatbit, vid1, curbit1): get2tupleid(vid1, curbit1, i, thatbit));
            long long tupleid_new = (i < vid1 ? get2tupleid(i, thatbit, vid1, t1): get2tupleid(vid1, t1, i, thatbit));
            --tuple_cov_cnt[tupleid_old];
            ++tuple_cov_cnt[tupleid_new];
        }
    }
    if (t2 != curbit2){
        for (int i = 0; i < nvar; ++i){
            if (i == vid1 || i == vid2) continue;
            int thatbit = tc[i];
            long long tupleid_old = (i < vid2 ? get2tupleid(i, thatbit, vid2, curbit2): get2tupleid(vid2, curbit2, i, thatbit));
            long long tupleid_new = (i < vid2 ? get2tupleid(i, thatbit, vid2, t2): get2tupleid(vid2, t2, i, thatbit));
            --tuple_cov_cnt[tupleid_old];
            ++tuple_cov_cnt[tupleid_new];
        }
    }

    long long tupleid_oldtp = get2tupleid(vid1, curbit1, vid2, curbit2);
    long long tupleid_newtp = get2tupleid(vid1, t1, vid2, t2);
    --tuple_cov_cnt[tupleid_oldtp];
    ++tuple_cov_cnt[tupleid_newtp];

    if (t1 != curbit1){
        if (__use_cell_tabu) last_flipped_time[tcid][vid1] = greedy_limit + 1;
        ca[tcid][vid1] = t1;
    }
    if (t2 != curbit2){
        if (__use_cell_tabu) last_flipped_time[tcid][vid2] = greedy_limit + 1;
        ca[tcid][vid2] = t2;
    }
}

void LSOptimizer::break_and_gain(){
    int uncovered_cnt = uncovered_tuples.size();
    for (const pair<int, int>& p: tmp_break){
        uncovered_tuples.push_back(p);
        int comp1 = p.first, comp2 = p.second;
        uncovered_tuple_positions[get2tupleid(abs(comp1) - 1, comp1 > 0 ? 1: 0, abs(comp2) - 1, comp2 > 0 ? 1: 0)] = uncovered_cnt;
        ++uncovered_cnt;
    }
    for (const pair<int, int>& p: tmp_gain){
        int comp1 = p.first, comp2 = p.second;
        long long this_tuple_id = get2tupleid(abs(comp1) - 1, comp1 > 0 ? 1: 0, abs(comp2) - 1, comp2 > 0 ? 1: 0);
        int& original_pos = uncovered_tuple_positions[this_tuple_id];
        if (original_pos != uncovered_cnt - 1){
            int tail_comp1 = uncovered_tuples[uncovered_cnt - 1].first, 
                tail_comp2 = uncovered_tuples[uncovered_cnt - 1].second;
            uncovered_tuple_positions[get2tupleid(abs(tail_comp1) - 1, tail_comp1 > 0 ? 1: 0, abs(tail_comp2) - 1, tail_comp2 > 0 ? 1: 0)] = original_pos;
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

    if (!trailing){
        for (int x: cur_clauses_cov){
            if (x == 0){
                std::cout << "c \033[1;31mError: SAT broken!\033[0m" << std::endl;
            }
        } 
    }
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
        long long tuple_pos = 0;
        for (int j = 0; j < nvar; ++j){
            long long tuple_id_base0 = testcase[j] ? num_combination_all_possible << 1: 0;
            for (int k = j + 1; k < nvar; ++k, ++tuple_pos){
                long long tuple_id_base1 = testcase[k] ? num_combination_all_possible: 0;
                long long tuple_id = tuple_id_base0 + tuple_id_base1 + tuple_pos;
                --tuple_cov_cnt[tuple_id];
            }
        }
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
    if (!validate_testcase(tc2)){
        std::cout << "c \033[1;31mError: SAT broken in pre testcase validation!\033[0m" << std::endl;
        return {-num_combination_all_possible - 1, 0};
    }

    const vector<int>& tc = ca[tcid];

    int break_cnt = 0, gain_cnt = 0;

    if (dump){
        tmp_break.clear();
        tmp_gain.clear();
    }

    long long tuple_pos = 0;
    for (int j = 0; j < nvar; ++j){
        long long tuple_id_base0_old = tc[j] ? num_combination_all_possible << 1: 0;
        long long tuple_id_base0_new = tc2[j] ? num_combination_all_possible << 1: 0;
        for (int k = j + 1; k < nvar; ++k, ++tuple_pos){
            long long tuple_id_base1_old = tc[k] ? num_combination_all_possible: 0;
            long long tuple_id_base1_new = tc2[k] ? num_combination_all_possible: 0;
            long long tuple_id_old = tuple_id_base0_old + tuple_id_base1_old + tuple_pos;
            long long tuple_id_new = tuple_id_base0_new + tuple_id_base1_new + tuple_pos;
            if (tuple_id_old != tuple_id_new){
                if (tuple_cov_cnt[tuple_id_old] == 1){
                    ++break_cnt;
                    if (dump){
                        tmp_break.emplace_back((tc[j] ? (j+1): (-j-1)), (tc[k] ? (k+1): (-k-1)));
                    }
                }
                if (tuple_cov_cnt[tuple_id_new] == 0){
                    ++gain_cnt;
                    if (dump){
                        tmp_gain.emplace_back((tc2[j] ? (j+1): (-j-1)), (tc2[k] ? (k+1): (-k-1)));
                    }
                }
            }
        }
    }

    return {break_cnt, gain_cnt};
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
    long long tuple_pos = 0;
    for (int j = 0; j < nvar; ++j){
        long long tuple_id_base0_old = tc[j] ? num_combination_all_possible << 1: 0;
        long long tuple_id_base0_new = tc2[j] ? num_combination_all_possible << 1: 0;
        for (int k = j + 1; k < nvar; ++k, ++tuple_pos){
            long long tuple_id_base1_old = tc[k] ? num_combination_all_possible: 0;
            long long tuple_id_base1_new = tc2[k] ? num_combination_all_possible: 0;
            long long tuple_id_old = tuple_id_base0_old + tuple_id_base1_old + tuple_pos;
            long long tuple_id_new = tuple_id_base0_new + tuple_id_base1_new + tuple_pos;
            if (tuple_id_old != tuple_id_new){
                --tuple_cov_cnt[tuple_id_old];
                ++tuple_cov_cnt[tuple_id_new];
            }
        }
    }

    if (__use_cell_tabu){
        for (int i = 0; i < nvar; ++i){
            if (tc[i] != tc2[i]){
                last_flipped_time[tcid][i] = greedy_limit + 1;
            }
        }
    }

    tc = tc2;
}
