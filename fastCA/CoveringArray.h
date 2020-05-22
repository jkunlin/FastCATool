#include <cmath>
#include <fstream>
#include <queue>
#include <set>
#include <sstream>
#include <string>
#include <vector>

#include <algorithm>
#include <atomic>
#include <condition_variable>
#include <fstream>
#include <functional>
#include <future>
#include <mutex>
#include <sys/time.h>
#include <thread>

#include "ConstraintFile.H"
#include "Coverage.h"
// #include "OptionTupleSet.h"
#include "LineVarTupleSet.h"
#include "SAT.H"
#include "Tabu.h"
#include "TestSetFile.H"
#include "ThreadTmpResult.h"
#include "TupleSet.h"
#include "Valid_check.h"
#include "mersenne.h"

class CoveringArray {
public:
  CoveringArray(const SpecificationFile &specificationFile,
                const ConstraintFile &constraintFile, TestSetFile &testSet,
                unsigned long long maxT, int seed, int threadsNum,
                int minScoreTaskSize, int minReplaceTaskSize,
                std::string outfile);
  ~CoveringArray();
  void actsInitialize(const std::string file_name);
  void optimize();

private:
  Valid::Validater validater;
  SATSolver satSolver;
  std::vector<bool> option_constrained_indicator;
  Mersenne mersenne;
  const SpecificationFile &specificationFile;
  TestSetFile &testSet;
  std::vector<std::vector<unsigned>> bestArray; // = array;
  std::vector<std::vector<unsigned>> array;
  Coverage coverage;
  TupleSet uncoveredTuples;
  std::set<unsigned> varInUncovertuples;
  LineVarTupleSet oneCoveredTuples;
  Tabu<Entry> entryTabu;
  std::string outfile;

  unsigned long long maxTime;
  unsigned long long step;
  struct timeval start_time;

  int minScoreTaskSize;
  int minReplaceTaskSize;
  std::vector<std::atomic<bool> *> taskReadyPtrs;
  std::vector<std::thread *> threadsPtr;
  std::vector<std::function<void()>> tasks;
  std::atomic<bool> programStop;
  int threadsNum;

  bool taskReady;
  int finishThreadNum;
  std::vector<std::mutex> taskMutex;
  std::vector<std::condition_variable> taskCv;

  std::mutex uncoveredTuplesMutex;
  std::mutex oneCoveredTuplesMutex;

  void cover(const unsigned encode, const unsigned oldLineIndex);
  void cover_with_lock(const unsigned encode, const unsigned oldLineIndex);
  void uncover(const unsigned encode, const unsigned oldLineIndex);
  void uncover_with_lock(const unsigned encode, const unsigned oldLineIndex);
  void updateTestSet();
  void replaceRow(const unsigned lineIndex, const unsigned encode);
  void replaceRowforTuple(const unsigned encode);
  void removeUselessRows();
  void removeOneRowRandom();
  long long varScoreOfRow(const unsigned var, const unsigned lineIndex);
  void replace(const unsigned var, const unsigned lineIndex);
  void replaceParallel(const unsigned int var, const unsigned int lineIndex);

  long long multiVarRow(const std::vector<unsigned> &sortedMultiVars,
                        const unsigned lineIndex, const bool change = false);
  long long multiVarScoreOfRow(const std::vector<unsigned> &sortedMultiVars,
                               const unsigned lineIndex);
  void multiVarReplace(const std::vector<unsigned> &sortedMultiVars,
                       const unsigned lineIndex);

  void tabugw();
  void tabugwParallel();
  void tabugwSubTask(const size_t start_index, const size_t end_index,
                     const unsigned &base, ThreadTmpResult &threadTmpResult);
  void tabugwReplaceSubTask(const unsigned &var, const unsigned &lineIndex,
                            size_t k, size_t count);
  void printBestArray() const;
  void outputBestArrayToFile() const;
  void tmpPrint();
  bool verify(const std::vector<std::vector<unsigned>> &resultArray);
  bool checkCovered(unsigned encode);
};
