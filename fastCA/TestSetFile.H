#ifndef TESTSET_FILE_H
#define TESTSET_FILE_H

#include "SpecificationFile.h"
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

// jkunlin, using vector instead of Array
class TestSetFile {
  std::vector<std::vector<int>> testSet;     // fastCA format
  std::vector<std::vector<int>> actsTestSet; // acts format
  int var_count;

public:
  TestSetFile(const std::string &filename);
  TestSetFile();
  unsigned int getTestSetCount();

  bool isEmpty() const { return testSet.size() == 0; }

  void convert2acts(const SpecificationFile &specificationFile);
  bool isExistedOption(const unsigned lineIndex, const unsigned option) const;
  int getSetSize() { return testSet.size(); };
  bool IsThisRow(unsigned inSetIndex, std::vector<unsigned> &rowInResult);
  void UpdateTestSetbyACTS(std::vector<std::vector<unsigned>> &array);

  std::string printInActsFormat() const;
};

#endif
