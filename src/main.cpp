#include "calocalsearch.h"
#include "clihelper.h"
#include <iostream>
#include <signal.h>
#include <chrono>

using namespace std;

void HandleInterrupt(int sig){
    cout << "c" << endl;
    cout << "c caught signal... exiting" << endl;
    exit(-1);
}

void SetupSignalHandler(){
    signal(SIGTERM, HandleInterrupt);
    signal(SIGINT, HandleInterrupt);
    signal(SIGQUIT, HandleInterrupt);
    signal(SIGKILL, HandleInterrupt);
}

bool ParseArgument(int argc, char **argv, Argument *argu){
    action_table act = build_default_actions();
    return match_by_action_table(act, argc, argv, argu);
}

int main(int argc, char **argv){
    SetupSignalHandler();

    Argument params;

    if (!ParseArgument(argc, argv, &params)){
        cout << "c Argument Error!" << endl;
        return -1;
    }

    LSOptimizer q(params);

    auto start_time = chrono::system_clock::now().time_since_epoch();
    q.search();
    auto end_time = chrono::system_clock::now().time_since_epoch();
    auto elapsed_time = chrono::duration_cast<chrono::milliseconds>(end_time - start_time).count() / 1000.0;
    cout << "c CPU time cost by optimization: " << elapsed_time << " seconds" << endl;

    vector<vector<int> > ca = q.get_testcase_set();
    int nvar = (int) ca[0].size();
    FILE *out_ca = fopen(params.output_testcase_path.c_str(), "w");
    for (const vector<int>& testcase: ca){
        for (int i = 0; i < nvar; ++i)
            fprintf(out_ca, "%d%c", testcase[i], " \n"[i == nvar - 1]);
    }
    fclose(out_ca);
    cout << "c optimized PCA saved in " << params.output_testcase_path << endl;

    return 0;
}
