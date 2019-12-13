#ifndef LOCALSEARCH_H
#define LOCALSEARCH_H

#include "ConstraintFile.H"
#include "SpecificationFile.h"
#include "TestSetFile.H"

void localSearch(const SpecificationFile &specificationFile,
                 const ConstraintFile &constraintFile,
                 TestSetFile &testSetFile, const unsigned long long maxTime, int seed,
                 int threadsNum, int minScoreTaskSize, int minReplaceTaskSize, std::string outfile);

#endif /* end of include guard: LOCALSEARCH_H */
