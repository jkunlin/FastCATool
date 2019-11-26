#include <iostream>
#include <string>

#include "ConstraintFile.H"
#include "LocalSearch.h"
#include "SpecificationFile.h"
#include "TestSetFile.H"

using namespace std;

int main(int argc, char const *argv[]) {
  if (argc == 0) {
    return 1;
  }
  string modelFile(argv[1]);
  string constrFile;
  string testFile("");
  unsigned long long maxTime;
  int seed;
  if (argc == 6) {
    constrFile = argv[2];
    testFile = argv[3];
    maxTime = atoi(argv[4]);
    seed = atoi(argv[5]);
  }
  else if (argc == 5) {
    constrFile = argv[2];
    maxTime = atoi(argv[3]);
    seed = atoi(argv[4]);
  } else {
    maxTime = atoi(argv[2]);
    seed = atoi(argv[3]);
  }
  SpecificationFile specificationFile(modelFile);
  ConstraintFile constraintFile(constrFile);
  TestSetFile testSetFile(testFile);
  testSetFile.convert2acts(specificationFile);
  localSearch(specificationFile, constrFile, maxTime, seed, testSetFile);
  return 0;
}
