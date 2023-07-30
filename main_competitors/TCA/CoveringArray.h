// =====================================================================================
//
//       Filename:  CoveringArray.h
//
//    Description:  
//
//        Version:  1.0
//        Created:  10/27/2014 09:18:19 PM
//       Revision:  none
//       Compiler:  g++
//
//         Author:  Jinkun Lin, jkunlin@gmail.com
//   Organization:  School of EECS, Peking University
//
// =====================================================================================

#include <vector>
#include <set>
#include <queue>
#include <cmath>
#include <string>
#include <fstream>
#include <sstream>

#include "ConstraintFile.H"
#include "Coverage.h"
#include "TupleSet.h"
#include "mersenne.h"
#include "Tabu.h"
#include "SAT.H"

class CoveringArray {
public:
	CoveringArray (const SpecificationFile &specificationFile,
			const ConstraintFile &constraintFile, unsigned long long maxT, int seed);
	void greedyConstraintInitialize();
  void arrayInitialize(const std::string file_name);
	void optimize();
private:
	SATSolver satSolver;
	Mersenne mersenne;
	const SpecificationFile &specificationFile;
	std::vector<std::vector<unsigned>> array;
	Coverage coverage;
	TupleSet uncoveredTuples;
	std::set<unsigned> varInUncovertuples;
	Tabu<Entry> entryTabu;

	unsigned long long maxTime;
	clock_t clock_start;
	
	void cover(const unsigned encode);
	void uncover(const unsigned encode);
	//produce one row at least cover one uncovered tuple.
	//Producing the row without update coverage
	void produceSatRow(std::vector<unsigned> &newLine, const unsigned encode);
	//greedily produce one row at least cover one uncovered tuple.
	//producing the row AND updating coverage
	void mostGreedySatRow(std::vector<unsigned> &newLine, const unsigned encode);
	void replaceRow(const unsigned lineIndex, const unsigned encode);
	void removeUselessRows();
	void removeOneRow();
	long long varScoreOfRow(const unsigned var, const unsigned lineIndex);
	void replace(const unsigned var, const unsigned lineIndex);

	long long multiVarRow(const std::vector<unsigned> &sortedMultiVars, const unsigned lineIndex,
			const bool change = false);
	long long multiVarScoreOfRow(const std::vector<unsigned> &sortedMultiVars,
			const unsigned lineIndex);
	void multiVarReplace(const std::vector<unsigned> &sortedMultiVars, const unsigned lineIndex);

	void tabuStep();
	void tmpPrint(std::vector<std::vector<unsigned>>& bestArray);
	bool verify(const std::vector<std::vector<unsigned>> &resultArray);
#ifndef NDEBUG 
	void print();
#endif
	void t();
};

