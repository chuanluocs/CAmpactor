#ifndef CNFINFO_INCLUDE_H
#define CNFINFO_INCLUDE_H

#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <iostream>

using std::vector;
using std::string;

class CNFInfo{
public:
    CNFInfo(): nvar(-1), nclauses(-1) {}
    
    CNFInfo(int nvar, int nclauses, 
        const vector<vector<int> >& clauses, 
        const vector<vector<int> >& pos_in_cls, 
        const vector<vector<int> >& neg_in_cls)
    : nvar(nvar), nclauses(nclauses), clauses(clauses), 
        pos_in_cls(pos_in_cls), neg_in_cls(neg_in_cls) {}
    
    CNFInfo(const string& cnf_input_path){
        if (!read_cnf(cnf_input_path)){
            nvar = nclauses = -1;
        }
    }

    void dump(int& _nvar, int& _nclauses, 
        vector<vector<int> >& _clauses, 
        vector<vector<int> >& _pos_in_cls, 
        vector<vector<int> >& _neg_in_cls
    ){
        _nvar = nvar;
        _nclauses = nclauses;
        _clauses = clauses;
        _pos_in_cls = pos_in_cls;
        _neg_in_cls = neg_in_cls;
    }

private:
    int nvar;
    int nclauses;
    vector<vector<int> > clauses;
    vector<vector<int> > pos_in_cls, neg_in_cls;

    bool read_cnf(const string& cnf_input_path){
        string line;
        std::istringstream iss;

        std::ifstream fin(cnf_input_path.c_str());
        if (!fin.is_open()) return false;

        while (getline(fin, line)){
            if (line.substr(0, 1) == "c")
                continue;
            else if (line.substr(0, 1) == "p"){
                string tempstr1, tempstr2;
                iss.clear();
                iss.str(line);
                iss.seekg(0, std::ios::beg);
                iss >> tempstr1 >> tempstr2 >> nvar >> nclauses;
                break;
            }
        }

        if (nvar < 0 || nclauses < 0){
            fin.close();
            return false;
        }

        pos_in_cls.resize(nvar + 1, vector<int>());
        neg_in_cls.resize(nvar + 1, vector<int>());
        clauses.resize(nclauses);

        for (int c = 0; c < nclauses; ++c){
            int cur_lit;
            fin >> cur_lit;
            while (cur_lit != 0){
                int v = abs(cur_lit);
                if (cur_lit > 0) pos_in_cls[v].push_back(c);
                else neg_in_cls[v].push_back(c);
                clauses[c].emplace_back(cur_lit);
                fin >> cur_lit;
            }
        }

        fin.close();

        return true;
    }
};

#endif
