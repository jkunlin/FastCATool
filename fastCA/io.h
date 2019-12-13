#ifndef IO_H

#define IO_H

#include "ConstraintFile.H"
#include "SpecificationFile.h"
#include "TestSetFile.H"
#include <fstream>
#include <string>
#include <unordered_map>
#include <vector>
/*! \class io
 *  \brief Brief class description
 *
 *  Detailed description
 */
class IO {
public:
  IO(){};
  virtual ~IO(){};

  void readInstance(const std::string &filename,
                    SpecificationFile &specification,
                    ConstraintFile &constraint, TestSetFile &testSet);

  const std::string &getValue(size_t varIndex, size_t valueIndex) const {
    return values[varIndex][valueIndex];
  }
  const std::string &getVar(size_t varIndex) const {
    return varNames[varIndex];
  }
  const std::string &getSystemName() const{
    return systemName;
  }

protected:
  void readSystemName();
  void readSpecification(SpecificationFile &specification);
  void readConstraint(ConstraintFile &constraint);
  void readTestSet(TestSetFile &testSet);
  Valid::Literal getLiteral(const std::string &term);
  unsigned getFastcaValue(std::string varName, std::string value);

  size_t varSize() const { return varNames.size(); }

  std::ifstream infile;

  std::string  systemName;

  std::vector<std::string> varNames; /*!< variable names */
  std::unordered_map<std::string, size_t> varName2index;

  std::vector<std::vector<std::string>> values; /*!< variable values */
  std::vector<std::unordered_map<std::string, size_t>> value2index;

  std::vector<unsigned> cumulativeValueCounts;
};
extern IO io;

#endif /* end of include guard: IO_H */
