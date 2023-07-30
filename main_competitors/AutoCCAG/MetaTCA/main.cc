#include <iostream>
#include <string>

#include "ConfigurationFile.h"
#include "ConstraintFile.H"
#include "LocalSearch.h"
#include "SpecificationFile.h"

using namespace std;

int main(int argc, char const *argv[]) {
  if (argc == 0) {
    return 1;
  }
  string modelFile;
  string constrFile;
  string configFile;
  string initCAFile;
  int seed = 1;
  int finalCheck = 1;

  for (int i=0; i<argc-1; ++i){
    if (string(argv[i]).compare(string("-model_file"))==0) {
      modelFile = argv[i+1];
    }
    if (string(argv[i]).compare(string("-constraints_file"))==0) {
      constrFile = argv[i+1];
    }
    if (string(argv[i]).compare(string("-configuration_file"))==0) {
      configFile = argv[i+1];
    }
    if (string(argv[i]).compare(string("-seed"))==0) {
      seed = atoi(argv[i+1]);
    }

    if (string(argv[i]).compare(string("-final_check"))==0) {
      finalCheck = atoi(argv[i+1]);
      if (finalCheck != 0 && finalCheck != 1) {
        std::cerr<< "illegal parameter: -final_check"<<std::endl;
        return 1;
      }
    }

    if (string(argv[i]).compare(string("-ca_file"))==0) {
      initCAFile = argv[i + 1];
    }
  }

  SpecificationFile specificationFile(modelFile);
  ConstraintFile constraintFile(constrFile);
  ConfigurationFile configurationFile(configFile);
  localSearch(specificationFile, constrFile, configurationFile, seed, finalCheck, initCAFile);
  return 0;
}
