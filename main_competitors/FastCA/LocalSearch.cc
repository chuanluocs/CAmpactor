#include "LocalSearch.h"
#include "TupleSet.h"

#include "CoveringArray.h"
#include <unistd.h>

void localSearch(const SpecificationFile &specificationFile,
                 const ConstraintFile &constraintFile,
                 const unsigned long long maxTime, int seed,
                 const std::string &initCAFile) {
  CoveringArray c(specificationFile, constraintFile, maxTime, seed);
  c.arrayInitialize(initCAFile);
  //	return ;
  c.optimize();
}
