#include "LocalSearch.h"
#include "TupleSet.h"

#include "ActsSolver.h"
#include "CoveringArray.h"
#include "io.h"
#include <unistd.h>
#include <random>

void localSearch(const SpecificationFile &specificationFile,
                 const ConstraintFile &constraintFile, TestSetFile &testSetFile,
                 const unsigned long long maxTime, int seed, int threadsNum,
                 int minScoreTaskSize, int minReplaceTaskSize, std::string outfile) {
  CoveringArray c(specificationFile, constraintFile, testSetFile, maxTime, seed,
                  threadsNum, minScoreTaskSize, minReplaceTaskSize, outfile);
  //	c.greedyConstraintInitialize2();
  ActsSolver ActsSolver;
//  char filename[L_tmpnam];
//  if (!mktemp(filename)) {
//    std::cerr << "tmp file name error" << std::endl;
//    abort();
//
  std::string acts_res_filename = "acts_result_" + io.getSystemName();
  acts_res_filename += std::to_string(getpid());
  ActsSolver.solve(specificationFile, constraintFile, acts_res_filename,
                   testSetFile);
  c.actsInitialize(acts_res_filename);
  std::string cmd = (std::string) "rm " + acts_res_filename;
  if (system(cmd.c_str()) != 0) {
    std::cerr << "can't remove acts result file" << std::endl;
    exit(0);
  };
  // return ;
  c.optimize();
}
