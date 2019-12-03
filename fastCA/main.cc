#include <iostream>
#include <map>
#include <string>

#include "ConstraintFile.H"
#include "LocalSearch.h"
#include "SpecificationFile.h"
#include "TestSetFile.H"

using namespace std;

int main(int argc, char const *argv[]) {
  if (argc < 3) {
    return 1;
  }
  if (argc % 2 != 1) {
    return 1;
  }

  string modelFile;
  string constrFile;
  unsigned long long maxTime;
  int seed;
  int threadsNum;
  string testFile("");

  map<string, string> parameters_map = {
      {"-f", ""}, {"-c", ""}, {"-t", "0"}, {"-s", "0"}, {"-p", "1"}};

  vector<string> parameterVec;
  for (int i = 1; i < argc - 1; i += 2) {
    string paraName = argv[i];
    string paraValue = argv[i + 1];
    if (parameters_map.find(paraName) == parameters_map.end()) {
      return 1;
    }
    parameters_map[paraName] = paraValue;
  }

  if (parameters_map["-f"] == "") {
    return 1;
  } else {
    modelFile = parameters_map["-f"];
  }

  constrFile = parameters_map["-c"];

  if (atoi(parameters_map["-t"].c_str()) < 0) {
    return 1;
  } else {
    maxTime = atoi(parameters_map["-t"].c_str());
  }

  seed = atoi(parameters_map["-s"].c_str());

  if (atoi(parameters_map["-p"].c_str()) < 1) {
    return 1;
  } else {
    threadsNum = atoi(parameters_map["-p"].c_str());
  }

  SpecificationFile specificationFile(modelFile);
  ConstraintFile constraintFile(constrFile);
  TestSetFile testSetFile(testFile);
  testSetFile.convert2acts(specificationFile);
  localSearch(specificationFile, constrFile, maxTime, seed, threadsNum,
              testSetFile);
  return 0;
}