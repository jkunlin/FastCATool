#include "io.h"
#include "Valid_check.h"
#include <cassert>
#include <cctype>
#include <cstdio>
#include <sstream>
#include <unordered_set>
#include <utility>
#include <vector>

IO io;

void strip(std::string &str) {
  if (str.empty()) {
    return;
  }
  size_t beg = 0, end = str.size();
  while (std::isspace(str[beg])) {
    ++beg;
  }
  while (std::isspace(str[end - 1])) {
    --end;
  }

  if (beg < end) {
    str = str.substr(beg, end);
  } else {
    str.clear();
  }
}

void IO::readInstance(const std::string &filename,
                      SpecificationFile &specification,
                      ConstraintFile &constraint, TestSetFile &testSet) {
  infile.open(filename);
  if (!infile) {
    std::cout << "can't open instance file" << std::endl;
    abort();
  }

  std::string line;
  while (getline(infile, line)) {
    if (line.find("[System]") != line.npos) {
      readSystemName();
    } else if (line.find("[Parameter]") != line.npos) {
      readSpecification(specification);
    } else if (line.find("[Constraint]") != line.npos) {
      readConstraint(constraint);
    } else if (line.find("[Test Set]") != line.npos) {
      readTestSet(testSet);
    } else {
      assert(false);
    }
  }
  infile.close();
}

// void IO::readSystemName() {
//   std::string line;
//   getline(infile, line);
//   std::istringstream is(line);
//
//   std::string tmp;
//   is >> tmp >> systemName;
// }

void IO::readSystemName() {
  auto lead_pos = infile.tellg();
  std::istringstream is;
  std::string line;
  while (getline(infile, line)) {
    if (line.find("[") != line.npos) {
      infile.seekg(lead_pos);
      break;
    }
    lead_pos = infile.tellg();
    if (line.empty() || (line[0] == '-' && line[1] == '-')) {
      continue;
    }
    is.clear();
    is.str(line);
    std::string tmp;
    is >> tmp >> systemName;
  }
}

void IO::readSpecification(SpecificationFile &specification) {
  std::string line;

  std::istringstream is;
  std::string varName, value, tmp;

  auto lead_pos = infile.tellg();
  unsigned cumValueCount = 0;
  std::vector<unsigned> value_counts;
  while (getline(infile, line)) {
    if (line.find("[") != line.npos) {
      infile.seekg(lead_pos);
      break;
    }
    lead_pos = infile.tellg();
    if (line.empty() || (line[0] == '-' && line[1] == '-')) {
      continue;
    }
    is.clear();
    is.str(line);

    getline(is, varName, '('); // passing "(type):"
    getline(is, tmp, ':');
    strip(varName);

    if (cumulativeValueCounts.empty()) {
      cumulativeValueCounts.push_back(0);
    } else {
      cumulativeValueCounts.push_back(cumValueCount);
    }
    varName2index.insert(std::make_pair(varName, varNames.size()));
    varNames.push_back(varName);
    values.push_back(std::vector<std::string>());
    value2index.push_back(std::unordered_map<std::string, size_t>());

    value_counts.push_back(0);
    while (getline(is, value, ',')) {
      strip(value);
      value2index.back().insert(std::make_pair(value, values.back().size()));
      values.back().push_back(value);
      value_counts.back()++;
    }
    cumValueCount += value_counts.back();
  }
  specification.initialize(value_counts);
}

Valid::Literal IO::getLiteral(const std::string &term) {
  Valid::Literal lit;

  bool isNeg = true;
  size_t pos = term.find("!=");
  if (pos == term.npos) {
    pos = term.find('=');
    isNeg = false;
  }

  std::string varName = term.substr(0, pos);
  strip(varName);

  pos += (isNeg ? 2 : 1);
  std::string value = term.substr(pos, term.size());
  strip(value);

  unsigned var = getFastcaValue(varName, value);

  lit.assign(isNeg, var);
  return lit;
}

void IO::readConstraint(ConstraintFile &constraint) {
  std::string line;
  std::string term;

  auto lead_pos = infile.tellg();
  Valid::Formula formula;
  Valid::Clause clause;
  while (getline(infile, line)) {
    clause.clear();
    if (line.find("[") != line.npos) {
      infile.seekg(lead_pos);
      break;
    }
    lead_pos = infile.tellg();
    if (line.empty() || (line[0] == '-' && line[1] == '-')) {
      continue;
    }
    std::size_t prev = 0, pos;
    while ((pos = line.find("||", prev)) != std::string::npos) {
      if (pos > prev) {
        term = line.substr(prev, pos - prev);
        strip(term);
        clause.push_back(getLiteral(term));
      }
      prev = pos + 2;
    }
    if (prev < line.length()) {
      term = line.substr(prev, std::string::npos);
      strip(term);
      clause.push_back(getLiteral(term));
    }
    formula.push_back(clause);
  }
  constraint.initClause(formula);
}

unsigned IO::getFastcaValue(std::string varName, std::string value) {
  unsigned res = 0;
  size_t varIndex = varName2index[varName];
  res += cumulativeValueCounts[varIndex];
  res += value2index[varIndex][value];
  return res;
}

void IO::readTestSet(TestSetFile &testSet) {
  testSet.setVarCount(varSize());
  std::string line;
  std::istringstream is;

  // read variable names
  std::vector<std::string> names;
  while (getline(infile, line)) {
    if (line.empty() || (line[0] == '-' && line[1] == '-')) {
      continue;
    }
    std::string name;
    is.clear();
    is.str(line);
    while (getline(is, name, ',')) {
      strip(name);
      names.push_back(name);
    }
    break;
  }

  std::vector<int> test(varNames.size(), -1);
  std::vector<size_t> pos;
  for (auto &varName : names) {
    size_t index = varName2index[varName];
    pos.push_back(index);
    test[index] = 0;
  }

  // read test set
  auto lead_pos = infile.tellg();
  while (getline(infile, line)) {
    if (line.find("[") != line.npos) {
      infile.seekg(lead_pos);
      return;
    }
    lead_pos = infile.tellg();
    if (line.empty() || (line[0] == '-' && line[1] == '-')) {
      continue;
    }
    std::string value;
    is.clear();
    is.str(line);
    for (size_t i = 0; i < names.size(); i++) {
      getline(is, value, ',');
      strip(value);
      if (value == "*") {
        test[pos[i]] = -1;
        continue;
      }
      auto &varName = names[i];
      test[pos[i]] = getFastcaValue(varName, value);
      //      std::cout << pos[i] << '\t' << test[pos[i]] << std::endl;
    }
    testSet.addTest(test);
  }
}
