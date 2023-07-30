#include <cmath>
#include <fstream>
#include <queue>
#include <set>
#include <sstream>
#include <string>
#include <vector>

#include "ConfigurationFile.h"
#include "ConstraintFile.H"
#include "Coverage.h"
#include "LineVarTupleSet.h"
#include "SAT.H"
#include "Tabu.h"
#include "TupleSet.h"
#include "mersenne.h"
#include "Valid_check.h"

class CoveringArray {
public:
  CoveringArray(const SpecificationFile &specificationFile,
                const ConstraintFile &constraintFile, 
                const ConfigurationFile &configurationFile, 
                int seed, const int finalChk);
  void greedyConstraintInitialize();
  // void actsInitialize(const std::string file_name);
  void arrayInitialize(const std::string file_name);
  void optimize();

private:
  Valid::Validater validater;
  SATSolver satSolver;
  std::vector<bool> option_constrained_indicator;
  Mersenne mersenne;
  const SpecificationFile &specificationFile;
  const ConfigurationFile &configurationFile;
  std::vector<std::vector<unsigned>> bestArray; // = array;
  std::vector<std::vector<unsigned>> array;
  Coverage coverage;
  TupleSet uncoveredTuples;
  std::set<unsigned> varInUncovertuples;
  LineVarTupleSet oneCoveredTuples;
  Tabu<Entry> entryTabu;

  unsigned long long maxTime;
  clock_t clock_start;

  long long step;
  int checkFunction;
  int scoringFunction;
  int gradientDecent;
  long modeBalance;
  int finalCheck;
  long tabuSize;
  int specialRandomOn;

  void cover(const unsigned encode, unsigned lineIndex);
  void uncover(const unsigned encode, unsigned lineIndex);
  // produce one row at least cover one uncovered tuple.
  // Producing the row without update coverage
  void produceSatRow(std::vector<unsigned> &newLine, const unsigned encode);
  // greedily produce one row at least cover one uncovered tuple.
  // producing the row AND updating coverage
  void mostGreedySatRow(const unsigned lineIndex, const unsigned encode);
  void replaceRow(const unsigned lineIndex, const unsigned encode);
  void replaceRow2(const unsigned lineIndex, const unsigned encode);
  void removeUselessRows();
  void removeUselessRows2();
  void removeOneRow();
  void removeOneRowRandom();
  long long varScoreOfRow(const unsigned var, const unsigned lineIndex);
  long long varScoreOfRow3(const unsigned var, const unsigned lineIndex);
  void replace(const unsigned var, const unsigned lineIndex);

  long long multiVarRow(const std::vector<unsigned> &sortedMultiVars,
                        const unsigned lineIndex, const bool change = false);
  long long multiVarScoreOfRow(const std::vector<unsigned> &sortedMultiVars,
                               const unsigned lineIndex);
  long long multiVarScoreOfRow2(const std::vector<unsigned> &sortedMultiVars,
                                   const unsigned lineIndex);
  void multiVarReplace(const std::vector<unsigned> &sortedMultiVars,
                       const unsigned lineIndex);
  void multiVarReplace2(const std::vector<unsigned> &sortedMultiVars,
                       const unsigned lineIndex);                       

  void tabuStep();
  void tabuStepHybrid();
  void tabuStepHybridGD();
  void tabugw();
  void tmpPrint();
  bool verify(const std::vector<std::vector<unsigned>> &resultArray);
  int checkInstance(const ConstraintFile &constraintFile);
#ifndef NDEBUG
  void print();
#endif
  void t();
};
