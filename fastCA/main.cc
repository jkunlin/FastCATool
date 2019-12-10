#include <iostream>
#include <map>
#include <string>

#include "ConstraintFile.H"
#include "LocalSearch.h"
#include "SpecificationFile.h"
#include "TestSetFile.H"
#include "io.h"
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
  int minScoreTaskSize;
  int minReplaceTaskSize;

  string testFile("");

  map<string, string> parameters_map = {{"-f", ""},
                                        {"-c", ""},
                                        {"-t", "0"},
                                        {"-s", "0"},
                                        {"-p", "1"},
                                        {"-minScoreTaskSize", "100"},
                                        {"-minReplaceTaskSize", "120"}};

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

  if (atoi(parameters_map["-minScoreTaskSize"].c_str()) < 1) {
    return 1;
  } else {
    minScoreTaskSize = atoi(parameters_map["-minScoreTaskSize"].c_str());
  }

  if (atoi(parameters_map["-minReplaceTaskSize"].c_str()) < 1) {
    return 1;
  } else {
    minReplaceTaskSize = atoi(parameters_map["-minReplaceTaskSize"].c_str());
  }

  SpecificationFile specificationFile;
  ConstraintFile constraintFile;
  TestSetFile testSetFile;
  io.readInstance(modelFile, specificationFile, constraintFile, testSetFile);
  specificationFile.setStrenth(3);

  testSetFile.convert2acts(specificationFile);
  localSearch(specificationFile, constraintFile, testSetFile, maxTime, seed,
              threadsNum, minScoreTaskSize, minReplaceTaskSize);

  return 0;
}
