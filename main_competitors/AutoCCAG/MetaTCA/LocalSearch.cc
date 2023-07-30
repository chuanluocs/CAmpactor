#include "LocalSearch.h"
#include "TupleSet.h"

#include "CoveringArray.h"
#include <unistd.h>

void localSearch(const SpecificationFile &specificationFile,
                 const ConstraintFile &constraintFile,
                 const ConfigurationFile &configurationFile,
                 const int seed, const int finalCheck,
                 const std::string &initCAFile) {

  CoveringArray c(specificationFile, constraintFile, configurationFile, seed, finalCheck); 
  // c.greedyConstraintInitialize();
  c.arrayInitialize(initCAFile);
  // ActsSolver ActsSolver;
  // char filename[L_tmpnam];
  // if (!tmpnam(filename)) {
  //   std::cerr << "tmp file name error" << std::endl;
  //   abort();
  // }
  // std::string acts_res_filename = filename;
  // acts_res_filename += std::to_string(getpid());
  // ActsSolver.solve(specificationFile, constraintFile, acts_res_filename);
  // c.actsInitialize(acts_res_filename);
  // std::string cmd = (std::string) "rm " + acts_res_filename;
  // if (system(cmd.c_str()) != 0) {
  //   std::cerr << "can't remove acts result file" << std::endl;
  //   exit(0);
  // };
  
  c.optimize();
}
