#include "SpecificationFile.h"
#include "TestSetFile.H"
#include <regex>

using namespace std;

TestSetFile::TestSetFile() {}

TestSetFile::TestSetFile(const string &filename) {
  if (!filename.size()) {
    return;
  }
  ifstream infile(filename.data());
  infile >> this->var_count;
  string line;
  regex ws_res(",");
  while (getline(infile, line)) {
    if (line.size() == 0)
      continue;
    vector<string> v(
        sregex_token_iterator(line.begin(), line.end(), ws_res, -1),
        sregex_token_iterator());
    vector<int> st;
    for (auto &&s : v) {
      if (s.compare("*") == 0) {
        st.push_back(-1);
      } else
        st.push_back(stoi(s));
    }
    this->testSet.push_back(st);
  }
  infile.close();
  // cout<<this->isEmpty()<<endl;
  // cout << "initial set size: "<<testSet.size() << endl;
}

void TestSetFile::convert2acts(const SpecificationFile &specificationFile) {
  const Options &options = specificationFile.getOptions();

  for (unsigned int i = 0; i < testSet.size(); i++) {
    vector<int> v;
    for (unsigned int j = 0; j < testSet[i].size(); j++) {
      int val = testSet[i][j];
      if (val < 0)
        v.push_back(-1);
      else {
        const unsigned option = options.option(val);
        v.push_back(val - options.firstSymbol(option));
      }
    }
    this->actsTestSet.push_back(v);
  }
}


string TestSetFile::printInActsFormat() const {
  ostringstream acts_infile;
  acts_infile << "p0";
  for (unsigned int i = 1; i < var_count; i++) {
    acts_infile << ",p" << i;
  }
  acts_infile << endl;

  //     for (unsigned int i = 0 ; i < contents.size(); i++) {
  //	     acts_infile << contents[i]<<endl;
  //
  //     }
  for (unsigned int i = 0; i < actsTestSet.size(); i++) {
    if (actsTestSet[i].size() == 0)
      continue;

    if (actsTestSet[i][0] == -1) {
      acts_infile << "*";
    } else {
      acts_infile << actsTestSet[i][0];
    }
    for (unsigned int j = 1; j < actsTestSet[i].size(); j++) {
      int val = actsTestSet[i][j];
      if (val == -1)
        acts_infile << ",*";
      else
        acts_infile << "," << val;
    }
    acts_infile << endl;
  }

  return acts_infile.str();
}

bool TestSetFile::isExistedOption(unsigned lineIndex, unsigned option) const {
  if (testSet.size() == 0 || lineIndex > testSet.size() - 1)
    return false;
  else if (testSet[lineIndex][option] < 0)
    return false;
  return true;
}

bool TestSetFile::isExistedRow(unsigned lineIndex) const {
  if (testSet.empty() || lineIndex > testSet.size() - 1)
    return false;
  return true;
}


bool TestSetFile::IsThisRow(unsigned inSetIndex,
                            vector<unsigned> &rowInResult) {
  for (unsigned i = 0; i < var_count; i++) {
    if (testSet[inSetIndex][i] > 0 && rowInResult[i] > 0 &&
        testSet[inSetIndex][i] != rowInResult[i])
      return false;
  }
  return true;
}

void TestSetFile::UpdateTestSetbyACTS(vector<vector<unsigned>> &array) {
  vector<vector<int>> validTestSet; // fastCA format
  unsigned scanArray = 0;

  for (unsigned i = 0; i < testSet.size() && scanArray < testSet.size(); i++) {
    if (IsThisRow(i, array[scanArray])) {
      validTestSet.push_back(testSet[i]);
      scanArray++;
    }
  }
  testSet.clear();
  testSet.assign(validTestSet.begin(), validTestSet.end());
  // cout<<"after clean:"<<testSet.size()<<endl;
}
