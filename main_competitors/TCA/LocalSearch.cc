// =====================================================================================
//
//       Filename:  LocalSearch.cc
//
//    Description:  
//
//        Version:  1.0
//        Created:  10/27/2014 02:55:57 PM
//       Revision:  none
//       Compiler:  g++
//
//         Author:  Jinkun Lin, jkunlin@gmail.com
//   Organization:  School of EECS, Peking University
//
// =====================================================================================

#include "LocalSearch.h"
#include "TupleSet.h"

#include "CoveringArray.h"

void localSearch(const SpecificationFile &specificationFile,
                 const ConstraintFile &constraintFile,
                 const unsigned long long maxTime, int seed,
                 const std::string &initCAFile) {
	CoveringArray c(specificationFile, constraintFile, maxTime, seed);
  c.arrayInitialize(initCAFile);
	c.optimize();
}
