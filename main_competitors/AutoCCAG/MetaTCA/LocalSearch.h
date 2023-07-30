#ifndef LOCALSEARCH_H
#define LOCALSEARCH_H

#include "ConstraintFile.H"
#include "SpecificationFile.h"
#include "ConfigurationFile.h"

void localSearch(const SpecificationFile &specificationFile,
                 const ConstraintFile &constrFile,
                 const ConfigurationFile &configurationFile,
                 int seed, const int finalCheck,
                 const std::string &initCAFile);

#endif /* end of include guard: LOCALSEARCH_H */
