// =====================================================================================
//
//       Filename:  main.cc
//
//    Description:  
//
//        Version:  1.0
//        Created:  10/27/2014 10:26:54 AM
//       Revision:  none
//       Compiler:  g++
//
//         Author:  Jinkun Lin, jkunlin@gmail.com
//   Organization:  School of EECS, Peking University
//
// =====================================================================================

#include <iostream>
#include <string>

#include "SpecificationFile.h"
#include "ConstraintFile.H"
#include "LocalSearch.h"


using namespace std;

int main(int argc, char const *argv[]) {
	if (argc == 0) {
		return 1;
	}
	string modelFile(argv[1]);
	string constrFile;
  string initCAFile;
	unsigned long long maxTime;
	int seed;
	if (argc == 6) {
		constrFile = argv[2];
		maxTime = atoi(argv[3]);
		seed = atoi(argv[4]);
    initCAFile = argv[5];
	}
	else {
		maxTime = atoi(argv[2]);
		seed = atoi(argv[3]);
    initCAFile = argv[4];
	}
	SpecificationFile specificationFile(modelFile);
	ConstraintFile constraintFile(constrFile);
  localSearch(specificationFile, constrFile, maxTime, seed, initCAFile);
	return 0;
}
