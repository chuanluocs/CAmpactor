#ifndef LOCALSEARCH_H
#define LOCALSEARCH_H

#include "ConstraintFile.H"
#include "SpecificationFile.h"

void localSearch(const SpecificationFile &specificationFile,
                 const ConstraintFile &constrFile,
                 const unsigned long long maxTime, int seed,
                 const std::string &initCAFile);

#endif /* end of include guard: LOCALSEARCH_H */
