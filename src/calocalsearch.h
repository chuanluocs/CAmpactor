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
#include <iostream>

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

    int twise;
    
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

    vector<vector<pair<int, int> > > uncovered_tuples;
    vector<pair<int, int> > uncovered_tuple_tmp;
    vector<int> uncovered_tuple_positions;

    std::mt19937_64 gen;
    int limit;
    
    int testcase_taboo;
    int greedy_limit;
    vector<int> last_greedy_time;

    vector<vector<pair<int, int> > > tmp_break, tmp_gain;
    
    vector<int> tmp_flip_sequence;

    int (LSOptimizer::*remove_testcase_strategy)();

    int remove_testcase_randomly();
    int remove_testcase_greedily();

    bool random_greedy_step();
    void greedy_step_forced(const vector<pair<int, int> >& chosen_tp);
    void random_step();

    Minisat::Solver *cdcl_solver;

    long long new_uncovered_tuples_after_remove_testcase(const vector<int>& tc, bool dump);

    pair<bool, pair<int, int> > get_gain_for_flip(int tcid, int vid, bool dump);
    void flip_bit(int tcid, int vid);

    pair<bool, pair<int, int> > get_gain_for_forcetuple(int tcid, const vector<pair<int, int> >& chosen_tp, bool dump);
    void forcetuple(int tcid, const vector<pair<int, int> >& chosen_tp);

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

    long long **combnum;

    long long get2base(int n, int i, int j){
        return (2ll * n - i - 1) * i / 2 + j - i - 1;
    }

    long long get3base(int n, int i1, int i2, int i3){
        long long res = combnum[3][n] - combnum[3][n - i1];
        return res + get2base(n - i1 - 1, i2 - i1 - 1, i3 - i1 - 1);
    }

    long long get4base(int n, int i1, int i2, int i3, int i4){
        long long res = combnum[4][n] - combnum[4][n - i1];
        return res + get3base(n - i1 - 1, i2 - i1 - 1, i3 - i1 - 1, i4 - i1 - 1);
    }

    long long get5base(int n, int i1, int i2, int i3, int i4, int i5){
        long long res = combnum[5][n] - combnum[5][n - i1];
        return res + get4base(n - i1 - 1, i2 - i1 - 1, i3 - i1 - 1, i4 - i1 - 1, i5 - i1 - 1);
    }

    long long get6base(int n, int i1, int i2, int i3, int i4, int i5, int i6){
        long long res = combnum[6][n] - combnum[6][n - i1];
        return res + get5base(n - i1 - 1, i2 - i1 - 1, i3 - i1 - 1, i4 - i1 - 1, i5 - i1 - 1, i6 - i1 - 1);
    }

    long long get2tupleid(int i, int vi, int j, int vj){
        long long base = (vi << 1 | vj) * num_combination_all_possible;
        return base + get2base(nvar, i, j);
    }

    long long get3tupleid(int i1, int v1, int i2, int v2, int i3, int v3){
        long long base = (v1 << 2 | v2 << 1 | v3) * num_combination_all_possible;
        return base + get3base(nvar, i1, i2, i3);
    }

    long long get4tupleid(int i1, int v1, int i2, int v2, int i3, int v3, int i4, int v4){
        long long base = (v1 << 3 | v2 << 2 | v3 << 1 | v4) * num_combination_all_possible;
        return base + get4base(nvar, i1, i2, i3, i4);
    }

    long long get5tupleid(int i1, int v1, int i2, int v2, int i3, int v3, int i4, int v4, int i5, int v5){
        long long base = (v1 << 4 | v2 << 3 | v3 << 2 | v4 << 1 | v5) * num_combination_all_possible;
        return base + get5base(nvar, i1, i2, i3, i4, i5);
    }

    long long get6tupleid(int i1, int v1, int i2, int v2, int i3, int v3, int i4, int v4, int i5, int v5, int i6, int v6){
        long long base = (v1 << 5 | v2 << 4 | v3 << 3 | v4 << 2 | v5 << 1 | v6) * num_combination_all_possible;
        return base + get6base(nvar, i1, i2, i3, i4, i5, i6);
    }

    int *tuple_iter_base;
    int target_twise;
    int tuple_iter_action_mode;

    vector<pair<int, int> > build_tuple_from_tuple_iter_base(const vector<int>& tc){
        vector<pair<int, int> > tp;
        tp.resize(target_twise);
        for (int i = 0; i < target_twise; ++i){
            tp[i].first = tuple_iter_base[i];
            tp[i].second = tc[tp[i].first];
        }
        return tp;
    }

    void iterate_all_tuples(int cur, long long tupleid, const vector<int>& tc);
    void iterate_all_tuples_double(int cur, long long tupleid_old, const vector<int>& tc_old, 
        long long tupleid_new, const vector<int>& tc_new);

    void tuple_iter_local_action_0(long long tuple_id){
        --tuple_cov_cnt[tuple_id];
    }

    void tuple_iter_local_action_1(long long tuple_id){
        ++tuple_cov_cnt[tuple_id];
    }

    bool new_uncovered_tuples_after_remove_testcase_dump;
    long long new_uncovered_tuples_after_remove_testcase_break_cnt;

    void tuple_iter_local_action_2(long long tuple_id, const vector<int>& tc){
        if (new_uncovered_tuples_after_remove_testcase_dump){
            int& tcc = tuple_cov_cnt[tuple_id];
            --tcc;
            if (tcc == 0){
                tmp_break.push_back(build_tuple_from_tuple_iter_base(tc));
            }
        } else {
            if (tuple_cov_cnt[tuple_id] == 1){
                ++new_uncovered_tuples_after_remove_testcase_break_cnt;
            }
        }
    }

    bool get_gain_for_forcetestcase_dump;
    long long get_gain_for_forcetestcase_break_cnt;
    long long get_gain_for_forcetestcase_gain_cnt;

    void tuple_iter_local_action_3(long long tuple_id_old, const vector<int>& tc_old, long long tuple_id_new, const vector<int>& tc_new){
        if (tuple_id_old == tuple_id_new) return ;
        if (get_gain_for_forcetestcase_dump){
            int& tcc_old = tuple_cov_cnt[tuple_id_old];
            --tcc_old;
            if (tcc_old == 0){
                tmp_break.push_back(build_tuple_from_tuple_iter_base(tc_old));
            }
            int& tcc_new = tuple_cov_cnt[tuple_id_new];
            ++tcc_new;
            if (tcc_new == 1){
                tmp_gain.push_back(build_tuple_from_tuple_iter_base(tc_new));
            }
        } else {
            if (tuple_cov_cnt[tuple_id_old] == 1){
                ++get_gain_for_forcetestcase_break_cnt;
            }
            if (tuple_cov_cnt[tuple_id_new] == 0){
                ++get_gain_for_forcetestcase_gain_cnt;
            }
        }
    }

    bool get_gain_for_flip_dump;
    int get_gain_for_flip_vid;
    int get_gain_for_flip_curbit;
    long long get_gain_for_flip_break_cnt;
    long long get_gain_for_flip_gain_cnt;

    vector<pair<int, int> > get_gain_for_flip_tuple_tmp;

    void tuple_iter_local_action_4(const vector<int>& tc){
        int cur = 0;
        while (cur < target_twise && tuple_iter_base[cur] <= get_gain_for_flip_vid){
            get_gain_for_flip_tuple_tmp[cur].first = tuple_iter_base[cur];
            get_gain_for_flip_tuple_tmp[cur].second = tc[tuple_iter_base[cur]];
            ++cur;
        }
        if (cur > 0 && get_gain_for_flip_tuple_tmp[cur - 1].first == get_gain_for_flip_vid){
            return ;
        }
        int lastcur = cur;
        get_gain_for_flip_tuple_tmp[cur].first = get_gain_for_flip_vid;
        get_gain_for_flip_tuple_tmp[cur].second = get_gain_for_flip_curbit;
        ++cur;
        while (cur < twise){
            get_gain_for_flip_tuple_tmp[cur].first = tuple_iter_base[cur - 1];
            get_gain_for_flip_tuple_tmp[cur].second = tc[tuple_iter_base[cur - 1]];
            ++cur;
        }

        if (get_gain_for_flip_dump){
            long long tupleid_old = getXtupleid(get_gain_for_flip_tuple_tmp);
            int& tcc_old = tuple_cov_cnt[tupleid_old];
            --tcc_old;
            if (tcc_old == 0){
                tmp_break.push_back(get_gain_for_flip_tuple_tmp);
            }

            get_gain_for_flip_tuple_tmp[lastcur].second ^= 1;

            long long tupleid_new = getXtupleid(get_gain_for_flip_tuple_tmp);
            int& tcc_new = tuple_cov_cnt[tupleid_new];
            ++tcc_new;
            if (tcc_new == 1){
                tmp_gain.push_back(get_gain_for_flip_tuple_tmp);
            }
        } else {
            long long tupleid_old = getXtupleid(get_gain_for_flip_tuple_tmp);
            if (tuple_cov_cnt[tupleid_old] == 1){
                ++get_gain_for_flip_break_cnt;
            }

            get_gain_for_flip_tuple_tmp[lastcur].second ^= 1;

            long long tupleid_new = getXtupleid(get_gain_for_flip_tuple_tmp);
            if (tuple_cov_cnt[tupleid_new] == 0){
                ++get_gain_for_flip_gain_cnt;
            }
        }
    }

    bool get_gain_for_forcetuple_dump;
    int get_gain_for_forcetuple_break_cnt;
    int get_gain_for_forcetuple_gain_cnt;

    vector<pair<int, int> > get_gain_for_forcetuple_chosen_tp;
    vector<pair<int, int> > get_gain_for_forcetuple_tuple_tmp;

    int get_gain_for_forcetuple_chosen_tp_ids[10];

    void tuple_iter_local_action_5(const vector<int>& tc){
        for (auto& pp: get_gain_for_forcetuple_chosen_tp){
            bool ok = true;
            for (int i = 0; i < target_twise; ++i){
                if (tuple_iter_base[i] == pp.first){
                    ok = false;
                    break;
                }
            }
            if (!ok){
                return ;
            }
        }

        for (int i = 0; i < (1 << twise); ++i){
            if (__builtin_popcount(i) + target_twise != twise){
                continue;
            }
            
            bool ok = false;
            int ntp = 0;
            for (int j = 0; j < twise; ++j){
                if ((i >> j) & 1){
                    get_gain_for_forcetuple_chosen_tp_ids[ntp] = j;
                    if (tc[get_gain_for_forcetuple_chosen_tp[j].first] != get_gain_for_forcetuple_chosen_tp[j].second){
                        ok = true;
                    }
                    ++ntp;
                }
            }
            
            if (!ok) continue;
            
            int qtp = 0, qfind = 0, cur = 0;
            while (qfind < target_twise && qtp < ntp){
                int rr = get_gain_for_forcetuple_chosen_tp_ids[qtp];
                if (get_gain_for_forcetuple_chosen_tp[rr].first < tuple_iter_base[qfind]){
                    get_gain_for_forcetuple_tuple_tmp[cur].first = get_gain_for_forcetuple_chosen_tp[rr].first;
                    ++qtp;
                } else {
                    get_gain_for_forcetuple_tuple_tmp[cur].first = tuple_iter_base[qfind];
                    ++qfind;
                }
                get_gain_for_forcetuple_tuple_tmp[cur].second = tc[get_gain_for_forcetuple_tuple_tmp[cur].first];
                ++cur;
            }
            while (qfind < target_twise){
                get_gain_for_forcetuple_tuple_tmp[cur].first = tuple_iter_base[qfind];
                ++qfind;
                get_gain_for_forcetuple_tuple_tmp[cur].second = tc[get_gain_for_forcetuple_tuple_tmp[cur].first];
                ++cur;
            }
            while (qtp < ntp){
                int rr = get_gain_for_forcetuple_chosen_tp_ids[qtp];
                get_gain_for_forcetuple_tuple_tmp[cur].first = get_gain_for_forcetuple_chosen_tp[rr].first;
                ++qtp;
                get_gain_for_forcetuple_tuple_tmp[cur].second = tc[get_gain_for_forcetuple_tuple_tmp[cur].first];
                ++cur;
            }
            long long tupleid_old = getXtupleid(get_gain_for_forcetuple_tuple_tmp);

            if (get_gain_for_forcetuple_dump){
                int& tcc_old = tuple_cov_cnt[tupleid_old];
                --tcc_old;
                if (tcc_old == 0){
                    tmp_break.push_back(get_gain_for_forcetuple_tuple_tmp);
                }
            } else {
                if (tuple_cov_cnt[tupleid_old] == 1){
                    ++get_gain_for_forcetuple_break_cnt;
                }
            }

            qtp = 0, qfind = 0, cur = 0;
            while (qfind < target_twise && qtp < ntp){
                int rr = get_gain_for_forcetuple_chosen_tp_ids[qtp];
                if (get_gain_for_forcetuple_chosen_tp[rr].first < tuple_iter_base[qfind]){
                    get_gain_for_forcetuple_tuple_tmp[cur].first = get_gain_for_forcetuple_chosen_tp[rr].first;
                    get_gain_for_forcetuple_tuple_tmp[cur].second = get_gain_for_forcetuple_chosen_tp[rr].second;
                    ++qtp;
                } else {
                    get_gain_for_forcetuple_tuple_tmp[cur].first = tuple_iter_base[qfind];
                    get_gain_for_forcetuple_tuple_tmp[cur].second = tc[get_gain_for_forcetuple_tuple_tmp[cur].first];
                    ++qfind;
                }
                ++cur;
            }
            while (qfind < target_twise){
                get_gain_for_forcetuple_tuple_tmp[cur].first = tuple_iter_base[qfind];
                get_gain_for_forcetuple_tuple_tmp[cur].second = tc[get_gain_for_forcetuple_tuple_tmp[cur].first];
                ++qfind;
                ++cur;
            }
            while (qtp < ntp){
                int rr = get_gain_for_forcetuple_chosen_tp_ids[qtp];
                get_gain_for_forcetuple_tuple_tmp[cur].first = get_gain_for_forcetuple_chosen_tp[rr].first;
                get_gain_for_forcetuple_tuple_tmp[cur].second = get_gain_for_forcetuple_chosen_tp[rr].second;
                ++qtp;
                ++cur;
            }
            long long tupleid_new = getXtupleid(get_gain_for_forcetuple_tuple_tmp);

            if (get_gain_for_forcetuple_dump){
                int& tcc_new = tuple_cov_cnt[tupleid_new];
                ++tcc_new;
                if (tcc_new == 1){
                    tmp_gain.push_back(get_gain_for_forcetuple_tuple_tmp);
                }
            } else {
                if (tuple_cov_cnt[tupleid_new] == 0){
                    ++get_gain_for_forcetuple_gain_cnt;
                }
            }
        }
    }
    
    long long getXbase(){
        switch (target_twise){
        case 1:
            return tuple_iter_base[0];
        case 2: 
            return get2base(nvar, tuple_iter_base[0], tuple_iter_base[1]);
        case 3: 
            return get3base(nvar, tuple_iter_base[0], tuple_iter_base[1], 
                            tuple_iter_base[2]);
        case 4: 
            return get4base(nvar, tuple_iter_base[0], tuple_iter_base[1], 
                            tuple_iter_base[2], tuple_iter_base[3]);
        case 5: 
            return get5base(nvar, tuple_iter_base[0], tuple_iter_base[1], 
                            tuple_iter_base[2], tuple_iter_base[3], 
                            tuple_iter_base[4]);
        case 6: 
            return get6base(nvar, tuple_iter_base[0], tuple_iter_base[1], 
                            tuple_iter_base[2], tuple_iter_base[3], 
                            tuple_iter_base[4], tuple_iter_base[5]);
        }
        return -1;
    }

    long long getXtupleid(const vector<pair<int, int> >& vec){
        switch (twise){
        case 2: 
            return get2tupleid(vec[0].first, vec[0].second,
                                            vec[1].first, vec[1].second);
        case 3: 
            return get3tupleid(vec[0].first, vec[0].second,
                                            vec[1].first, vec[1].second,
                                            vec[2].first, vec[2].second);
        case 4: 
            return get4tupleid(vec[0].first, vec[0].second,
                                            vec[1].first, vec[1].second,
                                            vec[2].first, vec[2].second,
                                            vec[3].first, vec[3].second);
        case 5: 
            return get5tupleid(vec[0].first, vec[0].second,
                                            vec[1].first, vec[1].second,
                                            vec[2].first, vec[2].second,
                                            vec[3].first, vec[3].second,
                                            vec[4].first, vec[4].second);
        case 6: 
            return get6tupleid(vec[0].first, vec[0].second,
                                            vec[1].first, vec[1].second,
                                            vec[2].first, vec[2].second,
                                            vec[3].first, vec[3].second,
                                            vec[4].first, vec[4].second,
                                            vec[5].first, vec[5].second);
        }
        return -1;
    }



};

#endif
