#ifndef LOCALSEARCH_H
#define LOCALSEARCH_H

#include "ConstraintFile.H"
#include "SpecificationFile.h"
#include "TestSetFile.H"

void localSearch(const SpecificationFile &specificationFile,
                 const ConstraintFile &constrFile,
                 const unsigned long long maxTime, int seed,
		 TestSetFile &testSetFile);

#endif /* end of include guard: LOCALSEARCH_H */
