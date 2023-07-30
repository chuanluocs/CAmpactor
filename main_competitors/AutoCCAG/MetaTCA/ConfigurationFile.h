#ifndef CONFIGURATION_H_CP9ABXVN
#define CONFIGURATION_H_CP9ABXVN
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

class ConfigurationFile {
public:
  ConfigurationFile(std::string filename) {
    std::ifstream infile(filename);
    if (!infile) {
      std::cerr << "can't open conf file!" << std::endl;
      exit(1);
    }
    std::string line, tmp;
    while (getline(infile, line)) {
      conf.push_back(std::vector<std::string>());
      std::istringstream is(line);
      while (is >> tmp) {
        conf.rbegin()->push_back(tmp);
      }
    }
  }

  size_t size() const { return conf.size(); }

  std::vector<std::vector<std::string>>::const_iterator begin() const {
    return conf.begin();
  }
  std::vector<std::vector<std::string>>::const_iterator end() const {
    return conf.end();
  }

private:
  std::vector<std::vector<std::string>> conf;
};

#endif /* end of include guard: CONFIGURATION_H_CP9ABXVN */
