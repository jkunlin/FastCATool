#include "CoveringArray.h"
#include "/home/jkunlin/cit/FastCATool/fastCA/Options.h"
#include "io.h"
CoveringArray::CoveringArray(const SpecificationFile &specificationFile,
                             const ConstraintFile &constraintFile,
                             TestSetFile &testSet, unsigned long long maxT,
                             int seed, int threadsNum, int minScoreTaskSize,
                             int minReplaceTaskSize)
    : validater(specificationFile), satSolver(constraintFile.isEmpty()),
      specificationFile(specificationFile), testSet(testSet),
      coverage(specificationFile), entryTabu(4), maxTime(maxT) {

  gettimeofday(&start_time, NULL);
  const Options &options = specificationFile.getOptions();
  // add constraint into satSolver
  const std::vector<InputClause> &clauses = constraintFile.getClauses();
  for (unsigned i = 0; i < clauses.size(); ++i) {
    satSolver.addClause(const_cast<InputClause &>(clauses[i]));
  }
  const Valid::Formula &formula = constraintFile.getFormula();
#ifndef NVISIBLE
  formula.print();
#endif
  for (auto &c : formula) {
    validater.addClause(c);
  }

  option_constrained_indicator.clear();
  option_constrained_indicator.resize(options.size(), false);
  for (auto &c : formula) {
    for (auto &lit : c) {
      option_constrained_indicator[options.option(lit.variable())] = true;
    }
  }
  for (unsigned option = 0; option < options.size(); ++option) {
    if (!option_constrained_indicator[option]) {
      continue;
    }
    InputClause atLeast;
    for (unsigned j = options.firstSymbol(option),
                  limit = options.lastSymbol(option);
         j <= limit; ++j) {
      atLeast.append(InputTerm(false, j));
    }
    satSolver.addClause(atLeast);
    for (unsigned j = options.firstSymbol(option),
                  limit = options.lastSymbol(option);
         j <= limit; ++j) {
      for (unsigned k = j + 1; k <= limit; ++k) {
        InputClause atMost;
        atMost.append(InputTerm(true, j));
        atMost.append(InputTerm(true, k));
        satSolver.addClause(atMost);
      }
    }
  }

  //		coverage.initialize(satSolver, option_constrained_indicator);
  //		uncoveredTuples.initialize(specificationFile, coverage, true);
  coverage.unconstrained_initialize();
  uncoveredTuples.initialize(specificationFile, coverage);

  mersenne.seed(seed);
  step = 0;

  this->threadsNum = threadsNum;
  this->minScoreTaskSize = minScoreTaskSize;
  this->minReplaceTaskSize = minReplaceTaskSize;

  // parallel
  if (this->threadsNum > 1) {
    tasks.resize(threadsNum);
    programStop.store(false);

    for (int t = 0; t < threadsNum; ++t) {
      taskReadyPtrs.push_back(new std::atomic<bool>(false));
    }

    for (int t = 1; t < threadsNum; ++t) {
      threadsPtr.push_back(new std::thread([t, this] {
        while (true) {
          while (true) {
            if (this->taskReadyPtrs[t]->load() || this->programStop.load())
              break;
          }
          if (this->taskReadyPtrs[t]->load())
            tasks[t]();
          if (this->programStop.load()) {
            break;
          }
          this->taskReadyPtrs[t]->store(false);
        }
      }));
    }
  }
}

CoveringArray::~CoveringArray() {
  // parallel
  if (this->threadsNum > 1) {
    this->programStop.store(true);

    for (auto threadptr : threadsPtr) {
      threadptr->join();
      delete threadptr;
    }

    for (auto taskReadyPtr : taskReadyPtrs) {
      delete taskReadyPtr;
    }
  }
}

void CoveringArray::actsInitialize(const std::string file_name) {
  const Options &options = specificationFile.getOptions();
  const unsigned &strenth = specificationFile.getStrenth();
  std::ifstream res_file(file_name);
  if (!res_file.is_open()) {
    std::cerr << "file open failed" << std::endl;
    exit(0);
  }
  std::string line;
  while (getline(res_file, line)) {
    if (line.find("Test Cases") != std::string::npos) {
      break;
    }
  }
  while (true) {
    bool begin = false;
    while (getline(res_file, line)) {
      if (line[0] == '1') {
        begin = true;
        break;
      }
    }
    if (!begin) {
      break;
    }
    array.push_back(std::vector<unsigned>(options.size()));
    std::vector<unsigned> &newRow = *array.rbegin();
    for (unsigned option = 0; option < options.size(); ++option) {
      unsigned value;
      std::string value_str(
          line.substr(line.find_last_of('=') + 1, line.size() - 1));
      value = atoi(value_str.c_str());
      newRow[option] = value + options.firstSymbol(option);
      getline(res_file, line);
    }
  }
  res_file.close();

  struct timeval start;
  gettimeofday(&start, NULL);
  std::vector<size_t> coverByLineIndex(coverage.tupleSize());
  std::vector<unsigned> tuple(strenth);
  for (unsigned lineIndex = 0; lineIndex < array.size(); ++lineIndex) {
    auto &line = array[lineIndex];
    for (std::vector<unsigned> columns = combinadic.begin(strenth);
         columns[strenth - 1] < line.size(); combinadic.next(columns)) {
      for (unsigned i = 0; i < strenth; ++i) {
        tuple[i] = line[columns[i]];
      }
      //			cover(coverage.encode(columns, tuple));
      unsigned encode = coverage.encode(columns, tuple);
      coverByLineIndex[encode] = lineIndex;
      coverage.cover(encode);
    }
  }
  coverage.set_zero_invalid();

  struct timeval end;
  gettimeofday(&end, NULL);
  double start_sec = start.tv_sec + start.tv_usec / 1000000.0;
  double end_sec = end.tv_sec + end.tv_usec / 1000000.0;
  std::cout << "actsInitialize: " << end_sec - start_sec << std::endl;

  entryTabu.initialize(Entry(array.size(), array.size()));
  validater.initialize(array);
  tmpPrint();
  oneCoveredTuples.initialize(specificationFile, array.size());
  oneCoveredTuples.pushOneCoveredTuple(coverage, coverByLineIndex);
#ifndef NDEBUG
  int i = 0;
  for (auto &line : array) {
    std::cout << i++ << ": ";
    for (auto v : line) {
      std::cout << v << ' ';
    }
    std::cout << std::endl;
  }
#endif
}

void CoveringArray::replaceRow(const unsigned lineIndex,
                               const unsigned encode) {
  std::vector<unsigned> &ranLine = array[lineIndex];
  const unsigned strength = specificationFile.getStrenth();
  std::vector<unsigned> tmpTuple(strength);
  // uncover the tuples
  for (std::vector<unsigned> columns = combinadic.begin(strength);
       columns[strength - 1] < ranLine.size(); combinadic.next(columns)) {
    for (unsigned i = 0; i < strength; ++i) {
      tmpTuple[i] = ranLine[columns[i]];
    }
    uncover(coverage.encode(columns, tmpTuple), lineIndex);
  }

  const std::vector<unsigned> &target_tuple = coverage.getTuple(encode);
  const std::vector<unsigned> &target_columns = coverage.getColumns(encode);
  for (auto &line : bestArray) {
    bool match = true;
    for (size_t i = 0; i < target_columns.size(); ++i) {
      if (line[target_columns[i]] != target_tuple[i]) {
        match = false;
        break;
      }
    }
    if (match) {
      ranLine = line;
      break;
    }
  }

  // cover the tuples
  for (std::vector<unsigned> columns = combinadic.begin(strength);
       columns[strength - 1] < ranLine.size(); combinadic.next(columns)) {
    for (unsigned i = 0; i < strength; ++i) {
      tmpTuple[i] = ranLine[columns[i]];
    }
    cover(coverage.encode(columns, tmpTuple), lineIndex);
  }
  entryTabu.initialize(
      Entry(array.size(), specificationFile.getOptions().size()));
  validater.change_row(lineIndex, ranLine);
}

void CoveringArray::removeUselessRows() {
  const Options &options = specificationFile.getOptions();
  const unsigned strength = specificationFile.getStrenth();
  std::vector<unsigned> tmpTuple(strength);

  for (size_t lineIndex = 0; lineIndex < array.size();) {
    if (oneCoveredTuples.oneCoveredCount(lineIndex) == 0) {
      const std::vector<unsigned> &line = array[lineIndex];
      for (std::vector<unsigned> columns = combinadic.begin(strength);
           columns[strength - 1] < options.size(); combinadic.next(columns)) {
        for (unsigned i = 0; i < strength; ++i) {
          tmpTuple[i] = line[columns[i]];
        }
        unsigned encode = coverage.encode(columns, tmpTuple);
        uncover(encode, lineIndex);
      }
      std::swap(array[lineIndex], array[array.size() - 1]);
      for (auto &entry : entryTabu) {
        if (entry.getRow() == array.size() - 1) {
          entry.setRow(lineIndex);
        }
        if (entry.getRow() == lineIndex) {
          entry.setRow(array.size() - 1);
        }
      }
      validater.exchange_row(lineIndex, array.size() - 1);
      validater.pop_back_row();
      oneCoveredTuples.exchange_row(lineIndex, array.size() - 1);
      oneCoveredTuples.pop_back_row();
      array.pop_back();
    } else {
      ++lineIndex;
    }
  }
}

void CoveringArray::removeOneRowRandom() {
  const Options &options = specificationFile.getOptions();
  const unsigned strength = specificationFile.getStrenth();

  unsigned rowToremoveIndex = mersenne.next(array.size());
  std::vector<unsigned> tmpTuple(strength);
  for (std::vector<unsigned> columns = combinadic.begin(strength);
       columns[strength - 1] < options.size(); combinadic.next(columns)) {
    for (unsigned j = 0; j < strength; ++j) {
      tmpTuple[j] = array[rowToremoveIndex][columns[j]];
    }
    unsigned encode = coverage.encode(columns, tmpTuple);
    uncover(encode, rowToremoveIndex);
  }

  std::swap(array[array.size() - 1], array[rowToremoveIndex]);
  validater.exchange_row(rowToremoveIndex, array.size() - 1);
  validater.pop_back_row();
  oneCoveredTuples.exchange_row(rowToremoveIndex, array.size() - 1);
  oneCoveredTuples.pop_back_row();
  for (auto &entry : entryTabu) {
    if (entry.getRow() == array.size() - 1) {
      entry.setRow(rowToremoveIndex);
    }
    if (entry.getRow() == rowToremoveIndex) {
      entry.setRow(array.size() - 1);
    }
  }
  array.pop_back();
}

void CoveringArray::updateTestSet() { testSet.UpdateTestSetbyACTS(array); }

void CoveringArray::optimize() {
  updateTestSet();
  while (true) {
    struct timeval end_time;
    gettimeofday(&end_time, NULL);
    double start_time_sec = start_time.tv_sec + start_time.tv_usec / 1000000.0;
    double end_time_sec = end_time.tv_sec + end_time.tv_usec / 1000000.0;
    if (end_time_sec - start_time_sec > maxTime) {
      break;
    }
    if (uncoveredTuples.size() == 0) {
      removeUselessRows();
      bestArray = array;
      tmpPrint();
      removeOneRowRandom();
    }

    size_t tasksNum = uncoveredTuples.size() * array.size();
    int neededThreadsNum = std::min(
        threadsNum, (int)std::ceil((double)tasksNum / minScoreTaskSize));
    if (threadsNum == 1 || neededThreadsNum < 2) {
      tabugw();
    } else {
      tabugwParallel();
    }

    struct timeval end;

    step++;
    continue;
  }

  if (uncoveredTuples.size() == 0) {
    removeUselessRows();
    bestArray = array;
    tmpPrint();
  }

  if (!verify(bestArray)) {
    std::cout << "wrong answer!!!!!" << std::endl;
    return;
  }
  std::cout << "total steps: " << step << std::endl;

  printBestArray();
}

void CoveringArray::printBestArray() const {
  const Options &options = specificationFile.getOptions();
  std::cout << "printing solution..." << std::endl;

  std::string sep;
  std::cout << "Parameters:" << std::endl;
  for (int opt = 0; opt < options.size(); opt++) {
    sep = "";
    std::cout << io.getVar(opt) << ": [";
    for (int symbol = options.firstSymbol(opt);
         symbol <= options.lastSymbol(opt); symbol++) {
      std::cout << sep << io.getValue(opt, symbol - options.firstSymbol(opt));
      sep = ", ";
    }
    std::cout << "]" << std::endl;
  }
  std::cout << std::endl;

  std::cout << "Configurations:" << std::endl;
  for (unsigned i = 0; i < bestArray.size(); ++i) {
    std::cout << i + 1 << "th  ";
    for (int j = 0; j < bestArray[i].size(); j++) {
      std::cout << ' '
                << io.getValue(j, bestArray[i][j] - options.firstSymbol(j));
    }
    std::cout << std::endl;
  }
  std::cout << "total size : " << bestArray.size() << std::endl;
}

void CoveringArray::tabugw() {
  unsigned base = mersenne.next(uncoveredTuples.size());
  std::vector<unsigned> firstBestRows;
  std::vector<unsigned> firstBestVars;
  std::vector<unsigned> bestRows;
  std::vector<unsigned> bestVars;
  long long bestScore = std::numeric_limits<long long>::min();
  for (size_t i = 0; i < uncoveredTuples.size(); ++i) {
    const unsigned tupleEncode =
        uncoveredTuples.encode((base + i) % uncoveredTuples.size());
    const std::vector<unsigned> &tuple = coverage.getTuple(tupleEncode);
    const std::vector<unsigned> &columns = coverage.getColumns(tupleEncode);
    for (unsigned lineIndex = 0; lineIndex < array.size(); ++lineIndex) {
      std::vector<unsigned> &line = array[lineIndex];
      unsigned diffCount = 0;
      unsigned diffVar;
      for (unsigned i = 0; i < tuple.size(); ++i) {
        if (line[columns[i]] != tuple[i]) {
          diffCount++;
          diffVar = tuple[i];
        }
      }
      if (diffCount > 1) {
        continue;
      }
      unsigned diffOption = specificationFile.getOptions().option(diffVar);
      // Tabu
      if (entryTabu.isTabu(Entry(lineIndex, diffOption))) {
        continue;
      }
      // my check
      if (!validater.valida_change(lineIndex, diffOption, line[diffOption],
                                   diffVar)) {
        continue;
      }
      long long tmpScore = varScoreOfRow(diffVar, lineIndex);
      if (bestScore < tmpScore) {
        bestScore = tmpScore;
        if (i == 0) {
          firstBestRows.clear();
          firstBestRows.push_back(lineIndex);
          firstBestVars.clear();
          firstBestVars.push_back(diffVar);
        }
        bestRows.clear();
        bestRows.push_back(lineIndex);
        bestVars.clear();
        bestVars.push_back(diffVar);
      } else if (bestScore == tmpScore) {
        bestRows.push_back(lineIndex);
        bestVars.push_back(diffVar);
        if (i == 0) {
          firstBestRows.push_back(lineIndex);
          firstBestVars.push_back(diffVar);
        }
      }
    }
  }
  if (bestScore > 0) {
    unsigned ran = mersenne.next(bestRows.size());
    replaceParallel(bestVars[ran], bestRows[ran]);
    return;
  }

  if (mersenne.nextClosed() < 0.001) {
    replaceRow(mersenne.next(array.size()), uncoveredTuples.encode(base));
    return;
  }
  if (firstBestRows.size() != 0) {
    unsigned ran = mersenne.next(firstBestRows.size());
    replaceParallel(firstBestVars[ran], firstBestRows[ran]);
    return;
  }

  const unsigned tupleEncode =
      uncoveredTuples.encode(base % uncoveredTuples.size());
  const std::vector<unsigned> &tuple = coverage.getTuple(tupleEncode);
  const std::vector<unsigned> &columns = coverage.getColumns(tupleEncode);

  bestRows.clear();

  std::vector<unsigned> changedVars;
  std::vector<unsigned> org_vars;
  for (unsigned lineIndex = 0; lineIndex < array.size(); ++lineIndex) {
    org_vars.clear();
    changedVars.clear();
    std::vector<unsigned> &line = array[lineIndex];
    for (unsigned i = 0; i < tuple.size(); ++i) {
      if (line[columns[i]] != tuple[i]) {
        org_vars.push_back(line[columns[i]]);
        changedVars.push_back(tuple[i]);
      }
    }
    if (changedVars.size() == 0) {
      continue;
    }
    // Tabu
    bool isTabu = false;
    for (auto v : changedVars) {
      unsigned diffOption = specificationFile.getOptions().option(v);
      if (entryTabu.isTabu(Entry(lineIndex, diffOption))) {
        isTabu = true;
        break;
      }
    }
    if (isTabu) {
      continue;
    }
    // check constraint, before tmpScore or after it?
    bool need_to_check = false;
    for (auto v : changedVars) {
      const Options &options = specificationFile.getOptions();
      if (option_constrained_indicator[options.option(v)]) {
        need_to_check = true;
        break;
      }
    }

    // my check
    if (need_to_check && !validater.valida_row(array[lineIndex], changedVars)) {
      continue;
    }
    // greedy
    long long tmpScore;
    if (changedVars.size() > 1) {
      tmpScore = multiVarScoreOfRow(changedVars, lineIndex);
    } else {
      tmpScore = varScoreOfRow(changedVars[0], lineIndex);
    }
    if (bestScore < tmpScore) {
      bestScore = tmpScore;
      bestRows.clear();
      bestRows.push_back(lineIndex);
    } else if (bestScore == tmpScore) {
      bestRows.push_back(lineIndex);
    }
  }
  // need to handle when bestRows.size() == 0
  if (bestRows.size() != 0) {
    // std::cout << "change_mutivar" << std::endl;
    unsigned lineIndex = bestRows[mersenne.next(bestRows.size())];
    changedVars.clear();
    for (unsigned i = 0; i < tuple.size(); ++i) {
      if (array[lineIndex][columns[i]] != tuple[i]) {
        changedVars.push_back(tuple[i]);
      }
    }
    if (changedVars.size() > 1) {
      // multiVarReplace(changedVars, lineIndex);
      for (auto v : changedVars) {
        replaceParallel(v, lineIndex);
      }
    } else {
      replaceParallel(changedVars[0], lineIndex);
    }
    return;
  }
  replaceRow(mersenne.next(array.size()), tupleEncode);
}

void CoveringArray::tabugwParallel() {
  unsigned base = mersenne.next(uncoveredTuples.size());
  std::vector<unsigned> firstBestRows;
  std::vector<unsigned> firstBestVars;
  std::vector<unsigned> bestRows;
  std::vector<unsigned> bestVars;
  long long bestScore = std::numeric_limits<long long>::min();

  size_t tasksNum = uncoveredTuples.size() * array.size();
  int neededThreadsNum =
      std::min(threadsNum, (int)std::ceil((double)tasksNum / minScoreTaskSize));
  std::vector<ThreadTmpResult> threadsTmpResult(neededThreadsNum);

  size_t taskSize = tasksNum / neededThreadsNum;
  int left = tasksNum % neededThreadsNum;
  size_t startIndex = 0;
  size_t endIndex;
  for (int t = 0; t < neededThreadsNum; ++t) {
    if (left > 0) {
      endIndex = startIndex + taskSize + 1;
      --left;
    } else {
      endIndex = startIndex + taskSize;
    }
    tasks[t] = std::packaged_task<void()>(
        [t, startIndex, endIndex, &base, &threadsTmpResult, this] {
          tabugwSubTask(startIndex, endIndex, base, threadsTmpResult[t]);
        });
    this->taskReadyPtrs[t]->store(true);
    startIndex = endIndex;
  }

  tasks[0]();

  // waitting the subtasks to be done
  while (true) {
    bool running = false;
    for (int i = 1; i < neededThreadsNum; ++i) {
      running = running || taskReadyPtrs[i]->load();
    }
    if (!running) {
      break;
    }
  }

  // merge the result of sub threads
  long long firstBestScore = std::numeric_limits<long long>::min();
  for (ThreadTmpResult &tmpRes : threadsTmpResult) {
    if (firstBestScore < tmpRes.firstBestScore) {
      firstBestScore = tmpRes.firstBestScore;
      firstBestRows.clear();
      firstBestVars.clear();
      for (int j = 0; j < tmpRes.firstBestRows.size(); ++j) {
        firstBestRows.push_back(tmpRes.firstBestRows[j]);
        firstBestVars.push_back(tmpRes.firstBestVars[j]);
      }
    } else if (firstBestScore == tmpRes.firstBestScore) {
      for (int j = 0; j < tmpRes.firstBestRows.size(); ++j) {
        firstBestRows.push_back(tmpRes.firstBestRows[j]);
        firstBestVars.push_back(tmpRes.firstBestVars[j]);
      }
    }

    if (bestScore < tmpRes.bestScore) {
      bestScore = tmpRes.bestScore;
      bestRows.clear();
      bestVars.clear();
      for (int j = 0; j < tmpRes.bestRows.size(); ++j) {
        bestRows.push_back(tmpRes.bestRows[j]);
        bestVars.push_back(tmpRes.bestVars[j]);
      }
    } else if (bestScore == tmpRes.bestScore) {
      for (int j = 0; j < tmpRes.bestRows.size(); ++j) {
        bestRows.push_back(tmpRes.bestRows[j]);
        bestVars.push_back(tmpRes.bestVars[j]);
      }
    }
  }

  if (bestScore > 0) {
    unsigned ran = mersenne.next(bestRows.size());
    replaceParallel(bestVars[ran], bestRows[ran]);
    return;
  }

  if (mersenne.nextClosed() < 0.001) {
    replaceRow(mersenne.next(array.size()), uncoveredTuples.encode(base));
    return;
  }
  if (firstBestRows.size() != 0) {
    unsigned ran = mersenne.next(firstBestRows.size());
    replaceParallel(firstBestVars[ran], firstBestRows[ran]);
    return;
  }

  const unsigned tupleEncode =
      uncoveredTuples.encode(base % uncoveredTuples.size());
  const std::vector<unsigned> &tuple = coverage.getTuple(tupleEncode);
  const std::vector<unsigned> &columns = coverage.getColumns(tupleEncode);

  bestRows.clear();

  std::vector<unsigned> changedVars;
  std::vector<unsigned> org_vars;
  for (unsigned lineIndex = 0; lineIndex < array.size(); ++lineIndex) {
    org_vars.clear();
    changedVars.clear();
    std::vector<unsigned> &line = array[lineIndex];
    for (unsigned i = 0; i < tuple.size(); ++i) {
      if (line[columns[i]] != tuple[i]) {
        org_vars.push_back(line[columns[i]]);
        changedVars.push_back(tuple[i]);
      }
    }
    if (changedVars.size() == 0) {
      continue;
    }
    // Tabu
    bool isTabu = false;
    for (auto v : changedVars) {
      unsigned diffOption = specificationFile.getOptions().option(v);
      if (entryTabu.isTabu(Entry(lineIndex, diffOption))) {
        isTabu = true;
        break;
      }
    }
    if (isTabu) {
      continue;
    }
    // check constraint, before tmpScore or after it?
    bool need_to_check = false;
    for (auto v : changedVars) {
      const Options &options = specificationFile.getOptions();
      if (option_constrained_indicator[options.option(v)]) {
        need_to_check = true;
        break;
      }
    }

    // my check
    if (need_to_check && !validater.valida_row(array[lineIndex], changedVars)) {
      continue;
    }
    // greedy
    long long tmpScore;
    if (changedVars.size() > 1) {
      tmpScore = multiVarScoreOfRow(changedVars, lineIndex);
    } else {
      tmpScore = varScoreOfRow(changedVars[0], lineIndex);
    }
    if (bestScore < tmpScore) {
      bestScore = tmpScore;
      bestRows.clear();
      bestRows.push_back(lineIndex);
    } else if (bestScore == tmpScore) {
      bestRows.push_back(lineIndex);
    }
  }
  // need to handle when bestRows.size() == 0
  if (bestRows.size() != 0) {
    // std::cout << "change_mutivar" << std::endl;
    unsigned lineIndex = bestRows[mersenne.next(bestRows.size())];
    changedVars.clear();
    for (unsigned i = 0; i < tuple.size(); ++i) {
      if (array[lineIndex][columns[i]] != tuple[i]) {
        changedVars.push_back(tuple[i]);
      }
    }
    if (changedVars.size() > 1) {
      // multiVarReplace(changedVars, lineIndex);
      for (auto v : changedVars) {
        replaceParallel(v, lineIndex);
      }
    } else {
      replaceParallel(changedVars[0], lineIndex);
    }
    return;
  }
  replaceRow(mersenne.next(array.size()), tupleEncode);
}

void CoveringArray::tabugwSubTask(const size_t start_index,
                                  const size_t end_index, const unsigned &base,
                                  ThreadTmpResult &threadTmpResult) {
  size_t index = start_index;
  size_t tupleIndex = index / array.size();
  size_t lineIndex;
  unsigned tupleEncode =
      uncoveredTuples.encode((base + tupleIndex) % uncoveredTuples.size());
  const std::vector<unsigned> *tuple = &coverage.getTuple(tupleEncode);
  const std::vector<unsigned> *columns = &coverage.getColumns(tupleEncode);

  for (; index < end_index; ++index) {
    lineIndex = index % array.size();
    if (lineIndex == 0) {
      tupleIndex = index / array.size();
      tupleEncode =
          uncoveredTuples.encode((base + tupleIndex) % uncoveredTuples.size());
      tuple = &coverage.getTuple(tupleEncode);
      columns = &coverage.getColumns(tupleEncode);
    }
    std::vector<unsigned> &line = array[lineIndex];
    unsigned diffCount = 0;
    unsigned diffVar;
    for (unsigned i = 0; i < tuple->size(); ++i) {
      if (line[columns->at(i)] != tuple->at(i)) {
        diffCount++;
        diffVar = tuple->at(i);
      }
    }
    if (diffCount > 1) {
      continue;
    }
    unsigned diffOption = specificationFile.getOptions().option(diffVar);
    // Tabu
    if (entryTabu.isTabu(Entry(lineIndex, diffOption))) {
      continue;
    }
    if (testSet.isExistedOption(lineIndex, diffOption)) {
      continue;
    }
    // my check
    if (!validater.valida_change(lineIndex, diffOption, line[diffOption],
                                 diffVar)) {
      continue;
    }
    long long tmpScore = varScoreOfRow(diffVar, lineIndex);
    if (threadTmpResult.bestScore < tmpScore) {
      threadTmpResult.bestScore = tmpScore;
      if (tupleIndex == 0) {
        threadTmpResult.firstBestScore = tmpScore;
        threadTmpResult.firstBestRows.clear();
        threadTmpResult.firstBestRows.push_back(lineIndex);
        threadTmpResult.firstBestVars.clear();
        threadTmpResult.firstBestVars.push_back(diffVar);
      }
      threadTmpResult.bestRows.clear();
      threadTmpResult.bestRows.push_back(lineIndex);
      threadTmpResult.bestVars.clear();
      threadTmpResult.bestVars.push_back(diffVar);
    } else if (threadTmpResult.bestScore == tmpScore) {
      threadTmpResult.bestRows.push_back(lineIndex);
      threadTmpResult.bestVars.push_back(diffVar);
      if (tupleIndex == 0) {
        threadTmpResult.firstBestRows.push_back(lineIndex);
        threadTmpResult.firstBestVars.push_back(diffVar);
      }
    }
  }
}

long long
CoveringArray::multiVarScoreOfRow(const std::vector<unsigned> &sortedMultiVars,
                                  const unsigned lineIndex) {
  const Options &options = specificationFile.getOptions();
  std::vector<unsigned> &line = array[lineIndex];
  std::vector<unsigned> changedColumns;
  for (auto var : sortedMultiVars) {
    int c = options.option(var);
    changedColumns.push_back(c);
  }

  long long coverChangeCount = 0;
  for (auto tupleEncode : uncoveredTuples) {
    const std::vector<unsigned> &tuple = coverage.getTuple(tupleEncode);
    const std::vector<unsigned> &columns = coverage.getColumns(tupleEncode);

    bool needChange = true;
    size_t i = 0, j = 0;
    while (i != columns.size() && j != changedColumns.size()) {
      if (columns[i] == changedColumns[j]) {
        if (tuple[i] != sortedMultiVars[j]) {
          needChange = false;
          break;
        }
        i++;
        j++;
      } else {
        if (columns[i] < changedColumns[j]) {
          if (tuple[i] != line[columns[i]]) {
            needChange = false;
            break;
          }
          i++;
        } else {
          j++;
        }
      }
    }
    for (; i != columns.size(); ++i) {
      if (tuple[i] != line[columns[i]]) {
        needChange = false;
        break;
      }
    }
    if (needChange) {
      coverChangeCount++;
    }
  }
  std::vector<long long> overlapCount(changedColumns.size() + 1, 0);
  for (auto cc : changedColumns) {
    for (auto &ecEntry : oneCoveredTuples.getECbyLineVar(lineIndex, line[cc])) {
      unsigned tupleEncode = ecEntry.encode;
      const std::vector<unsigned> &tuple = coverage.getTuple(tupleEncode);
      const std::vector<unsigned> &columns = coverage.getColumns(tupleEncode);
      bool needChange = true;
      for (size_t i = 0; i < columns.size(); ++i) {
        if (line[columns[i]] != tuple[i]) {
          needChange = false;
          break;
        }
      }
      if (needChange) {
        long long olc = 0;
        size_t j = 0, k = 0;
        while (j != columns.size() && k != changedColumns.size()) {
          if (columns[j] == changedColumns[k]) {
            olc++;
            j++;
            k++;
          } else {
            columns[j] < changedColumns[k] ? j++ : k++;
          }
        }
        overlapCount[olc]++;
      }
    }
  }
  for (size_t olc = 1; olc < overlapCount.size(); ++olc) {
    coverChangeCount -= overlapCount[olc] / olc;
  }
  return coverChangeCount;
}

long long
CoveringArray::multiVarRow(const std::vector<unsigned> &sortedMultiVars,
                           const unsigned lineIndex, const bool change) {
  const Options &options = specificationFile.getOptions();
  const unsigned strength = specificationFile.getStrenth();
  long long score = 0;

  std::vector<unsigned> varColumns;
  for (auto var : sortedMultiVars) {
    varColumns.push_back(options.option(var));
  }
  if (change) {
    for (auto column : varColumns) {
      entryTabu.insert(Entry(lineIndex, column));
    }
  }
  std::vector<unsigned> &line = array[lineIndex];

  // must from the end to the begining
  for (unsigned i = sortedMultiVars.size(); i--;) {
    std::swap(line[line.size() - sortedMultiVars.size() + i],
              line[varColumns[i]]);
  }

  std::vector<unsigned> tmpSortedColumns(strength);
  std::vector<unsigned> tmpSortedTupleToCover(strength);
  std::vector<unsigned> tmpSortedTupleToUncover(strength);
  unsigned tmpToCoverEncode;
  unsigned tmpToUncoverEncode;

  if (sortedMultiVars.size() >= strength) {
    for (std::vector<unsigned> changedColums = combinadic.begin(strength);
         changedColums[strength - 1] < sortedMultiVars.size();
         combinadic.next(changedColums)) {
      for (unsigned i = 0; i < strength; ++i) {
        tmpSortedTupleToCover[i] = sortedMultiVars[changedColums[i]];
        tmpSortedTupleToUncover[i] =
            line[line.size() - sortedMultiVars.size() + changedColums[i]];
      }
      std::sort(tmpSortedTupleToCover.begin(), tmpSortedTupleToCover.end());
      std::sort(tmpSortedTupleToUncover.begin(), tmpSortedTupleToUncover.end());
      for (unsigned i = 0; i < strength; ++i) {
        tmpSortedColumns[i] = options.option(tmpSortedTupleToCover[i]);
      }
      tmpToCoverEncode =
          coverage.encode(tmpSortedColumns, tmpSortedTupleToCover);
      tmpToUncoverEncode =
          coverage.encode(tmpSortedColumns, tmpSortedTupleToUncover);
      if (change) {
        cover(tmpToCoverEncode, lineIndex);
        uncover(tmpToUncoverEncode, lineIndex);
      } else {
        if (coverage.coverCount(tmpToCoverEncode) == 0) {
          ++score;
        }
        if (coverage.coverCount(tmpToUncoverEncode) == 1) {
          --score;
        }
      }
    }
  }

  for (unsigned curRelevantCount = 1,
                maxRelevantCount = std::min(
                    strength - 1,
                    (const unsigned)(line.size() - sortedMultiVars.size()));
       curRelevantCount <= maxRelevantCount; ++curRelevantCount) {
    for (std::vector<unsigned> relevantColumns =
             combinadic.begin(curRelevantCount);
         relevantColumns[curRelevantCount - 1] <
         line.size() - sortedMultiVars.size();
         combinadic.next(relevantColumns)) {
      for (std::vector<unsigned> changedColums =
               combinadic.begin(strength - curRelevantCount);
           changedColums[strength - curRelevantCount - 1] <
           sortedMultiVars.size();
           combinadic.next(changedColums)) {

        for (unsigned i = 0; i < curRelevantCount; ++i) {
          tmpSortedTupleToCover[i] = tmpSortedTupleToUncover[i] =
              line[relevantColumns[i]];
        }
        for (unsigned i = 0; i < strength - curRelevantCount; ++i) {
          tmpSortedTupleToCover[curRelevantCount + i] =
              sortedMultiVars[changedColums[i]];
          tmpSortedTupleToUncover[curRelevantCount + i] =
              line[line.size() - sortedMultiVars.size() + changedColums[i]];
        }
        std::sort(tmpSortedTupleToCover.begin(), tmpSortedTupleToCover.end());
        std::sort(tmpSortedTupleToUncover.begin(),
                  tmpSortedTupleToUncover.end());
        for (unsigned i = 0; i < strength; ++i) {
          tmpSortedColumns[i] = options.option(tmpSortedTupleToCover[i]);
        }
        tmpToCoverEncode =
            coverage.encode(tmpSortedColumns, tmpSortedTupleToCover);
        tmpToUncoverEncode =
            coverage.encode(tmpSortedColumns, tmpSortedTupleToUncover);
        if (change) {
          cover(tmpToCoverEncode, lineIndex);
          uncover(tmpToUncoverEncode, lineIndex);
        } else {
          if (coverage.coverCount(tmpToCoverEncode) == 0) {
            ++score;
          }
          if (coverage.coverCount(tmpToUncoverEncode) == 1) {
            --score;
          }
        }
      }
    }
  }
  // must from the begining to the end
  for (unsigned i = 0; i < sortedMultiVars.size(); ++i) {
    std::swap(line[line.size() - sortedMultiVars.size() + i],
              line[varColumns[i]]);
  }

  if (change) {
    for (unsigned i = 0; i < sortedMultiVars.size(); ++i) {
      line[varColumns[i]] = sortedMultiVars[i];
    }
  }
  return score;
}

void CoveringArray::multiVarReplace(
    const std::vector<unsigned> &sortedMultiVars, const unsigned lineIndex) {
  const Options &options = specificationFile.getOptions();
  std::vector<unsigned> org_vars;
  for (auto var : sortedMultiVars) {
    auto option = options.option(var);
    org_vars.push_back(array[lineIndex][option]);
  }
  validater.change_mutivar(lineIndex, org_vars, sortedMultiVars);
  multiVarRow(sortedMultiVars, lineIndex, true);
}

long long CoveringArray::varScoreOfRow(const unsigned var,
                                       const unsigned lineIndex) {
  const Options &options = specificationFile.getOptions();
  std::vector<unsigned> &line = array[lineIndex];
  const unsigned varOption = options.option(var);
  if (line[varOption] == var) {
    return 0;
  }
  long long coverChangeCount = 0;
  for (auto tupleEncode : uncoveredTuples) {
    const std::vector<unsigned> &tuple = coverage.getTuple(tupleEncode);
    const std::vector<unsigned> &columns = coverage.getColumns(tupleEncode);
    bool match = false;
    bool needChange = true;
    for (size_t i = 0; i < columns.size(); ++i) {
      if (columns[i] == varOption) {
        match = true;
        if (tuple[i] != var) {
          needChange = false;
          break;
        }
      } else if (line[columns[i]] != tuple[i]) {
        needChange = false;
        break;
      }
    }
    if (match && needChange) {
      coverChangeCount++;
    }
  }
  return coverChangeCount -
         oneCoveredTuples.getECbyLineVar(lineIndex, line[varOption]).size();
  for (auto &ecEntry :
       oneCoveredTuples.getECbyLineVar(lineIndex, line[varOption])) {
    unsigned tupleEncode = ecEntry.encode;
    const std::vector<unsigned> &tuple = coverage.getTuple(tupleEncode);
    const std::vector<unsigned> &columns = coverage.getColumns(tupleEncode);
    bool needChange = true;
    for (size_t i = 0; i < columns.size(); ++i) {
      if (line[columns[i]] != tuple[i]) {
        needChange = false;
        break;
      }
    }
    if (needChange) {
      coverChangeCount--;
    }
  }
  return coverChangeCount;
}
bool CoveringArray::checkCovered(unsigned encode) {
  for (auto ec : uncoveredTuples) {
    if (ec == encode) {
      return true;
    }
  }
  return false;
}

void CoveringArray::replace(const unsigned var, const unsigned lineIndex) {
  const Options &options = specificationFile.getOptions();
  const unsigned strength = specificationFile.getStrenth();
  std::vector<unsigned> &line = array[lineIndex];
  const unsigned varOption = options.option(var);

  entryTabu.insert(Entry(lineIndex, varOption));

  if (line[varOption] == var) {
    return;
  }
  validater.change_var(lineIndex, varOption, line[varOption], var);
  std::swap(line[line.size() - 1], line[varOption]);

  std::vector<unsigned> tmpSortedColumns(strength);
  std::vector<unsigned> tmpSortedTupleToCover(strength);
  std::vector<unsigned> tmpSortedTupleToUncover(strength);
  for (std::vector<unsigned> columns = combinadic.begin(strength - 1);
       columns[strength - 2] < line.size() - 1; combinadic.next(columns)) {
    for (unsigned i = 0; i < columns.size(); ++i) {
      tmpSortedTupleToUncover[i] = tmpSortedTupleToCover[i] = line[columns[i]];
    }
    tmpSortedTupleToCover[strength - 1] = var;
    tmpSortedTupleToUncover[strength - 1] = line[line.size() - 1];
    std::sort(tmpSortedTupleToCover.begin(), tmpSortedTupleToCover.end());
    std::sort(tmpSortedTupleToUncover.begin(), tmpSortedTupleToUncover.end());
    for (unsigned i = 0; i < tmpSortedTupleToCover.size(); ++i) {
      tmpSortedColumns[i] = options.option(tmpSortedTupleToCover[i]);
    }
    unsigned tmpTupleToCoverEncode =
        coverage.encode(tmpSortedColumns, tmpSortedTupleToCover);
    unsigned tmpTupleToUncoverEncode =
        coverage.encode(tmpSortedColumns, tmpSortedTupleToUncover);
    cover(tmpTupleToCoverEncode, lineIndex);
    uncover(tmpTupleToUncoverEncode, lineIndex);
  }
  std::swap(line[line.size() - 1], line[varOption]);
  line[varOption] = var;
}

void CoveringArray::replaceParallel(const unsigned var,
                                    const unsigned lineIndex) {
  std::vector<unsigned> &line = array[lineIndex];
  const Options &options = specificationFile.getOptions();
  const unsigned strength = specificationFile.getStrenth();
  const unsigned varOption = options.option(var);

  size_t taskNum = 0;
  taskNum = pt.nCr(line.size() - 1, strength - 1);

  //  std::cout << taskNum << std::endl;
  int neededThreadsNum = std::min(
      threadsNum, (int)std::ceil((double)taskNum / minReplaceTaskSize));
  if (threadsNum == 1 || neededThreadsNum < 2) {
    replace(var, lineIndex);
    return;
  }

  entryTabu.insert(Entry(lineIndex, varOption));

  if (line[varOption] == var) {
    return;
  }

  validater.change_var(lineIndex, varOption, line[varOption], var);
  std::swap(line[line.size() - 1], line[varOption]);

  size_t taskSize = taskNum / neededThreadsNum;
  int left = taskNum % neededThreadsNum;
  size_t count;
  std::vector<unsigned> columns = combinadic.begin(strength - 1);
  columns = combinadic.begin(strength - 1);
  for (int t = 0; t < neededThreadsNum; ++t) {
    if (left > 0) {
      count = taskSize + 1;
      --left;
    } else {
      count = taskSize;
    }

    tasks[t] = std::packaged_task<void()>(
        [&var, &lineIndex, &options, &strength, &line, columns, count, this] {
          this->tabugwReplaceSubTask(var, lineIndex, options, strength, line,
                                     columns, count);
        });
    taskReadyPtrs[t]->store(true);

    for (int i = 0; i < count; ++i) {
      combinadic.next(columns);
    }
  }

  tasks[0]();

  // waitting the subtasks to be done
  while (true) {
    bool running = false;
    for (int i = 1; i < neededThreadsNum; ++i) {
      running = running || taskReadyPtrs[i]->load();
    }
    if (!running) {
      break;
    }
  }

  std::swap(line[line.size() - 1], line[varOption]);
  line[varOption] = var;
}

void CoveringArray::tabugwReplaceSubTask(
    const unsigned &var, const unsigned &lineIndex, const Options &options,
    const unsigned &strength, std::vector<unsigned> &line,
    std::vector<unsigned> columns, size_t count) {
  {
    std::vector<unsigned> tmpSortedColumns(strength);
    std::vector<unsigned> tmpSortedTupleToCover(strength);
    std::vector<unsigned> tmpSortedTupleToUncover(strength);
    for (size_t done = 0; done < count; ++done, combinadic.next(columns)) {
      for (unsigned i = 0; i < columns.size(); ++i) {
        tmpSortedTupleToUncover[i] = tmpSortedTupleToCover[i] =
            line[columns[i]];
      }

      tmpSortedTupleToCover[strength - 1] = var;
      tmpSortedTupleToUncover[strength - 1] = line[line.size() - 1];
      std::sort(tmpSortedTupleToCover.begin(), tmpSortedTupleToCover.end());
      std::sort(tmpSortedTupleToUncover.begin(), tmpSortedTupleToUncover.end());
      for (unsigned i = 0; i < tmpSortedTupleToCover.size(); ++i) {
        tmpSortedColumns[i] = options.option(tmpSortedTupleToCover[i]);
      }
      unsigned tmpTupleToCoverEncode =
          coverage.encode(tmpSortedColumns, tmpSortedTupleToCover);
      unsigned tmpTupleToUncoverEncode =
          coverage.encode(tmpSortedColumns, tmpSortedTupleToUncover);
      // need not check coverCount, cover(encode) will do this
      cover_with_lock(tmpTupleToCoverEncode, lineIndex);
      uncover_with_lock(tmpTupleToUncoverEncode, lineIndex);
    }
  }
}

void CoveringArray::cover(const unsigned encode, const unsigned oldLineIndex) {
  coverage.cover(encode);
  unsigned coverCount = coverage.coverCount(encode);
  if (coverCount == 1) {
    uncoveredTuples.pop(encode);
    oneCoveredTuples.push(encode, oldLineIndex, coverage.getTuple(encode));
  }
  if (coverCount == 2) {
    const std::vector<unsigned> &tuple = coverage.getTuple(encode);
    const std::vector<unsigned> &columns = coverage.getColumns(encode);
    for (size_t lineIndex = 0; lineIndex < array.size(); ++lineIndex) {
      if (lineIndex == oldLineIndex) {
        continue;
      }
      auto &line = array[lineIndex];
      bool match = true;
      for (size_t i = 0; i < columns.size(); ++i) {
        if (tuple[i] != line[columns[i]]) {
          match = false;
          break;
        }
      }
      if (match) {
        oneCoveredTuples.pop(encode, lineIndex, tuple);
        break;
      }
    }
  }
}

void CoveringArray::cover_with_lock(const unsigned encode,
                                    const unsigned oldLineIndex) {
  coverage.cover(encode);
  unsigned coverCount = coverage.coverCount(encode);
  if (coverCount == 1) {
    {
      std::lock_guard<std::mutex> lg(uncoveredTuplesMutex);
      uncoveredTuples.pop(encode);
    }
    {
      std::lock_guard<std::mutex> lg(oneCoveredTuplesMutex);
      oneCoveredTuples.push(encode, oldLineIndex, coverage.getTuple(encode));
    }
  }
  if (coverCount == 2) {
    const std::vector<unsigned> &tuple = coverage.getTuple(encode);
    const std::vector<unsigned> &columns = coverage.getColumns(encode);
    for (size_t lineIndex = 0; lineIndex < array.size(); ++lineIndex) {
      if (lineIndex == oldLineIndex) {
        continue;
      }
      auto &line = array[lineIndex];
      bool match = true;
      for (size_t i = 0; i < columns.size(); ++i) {
        if (tuple[i] != line[columns[i]]) {
          match = false;
          break;
        }
      }
      if (match) {
        {
          std::lock_guard<std::mutex> lg(oneCoveredTuplesMutex);
          oneCoveredTuples.pop(encode, lineIndex, tuple);
        }
        break;
      }
    }
  }
}

void CoveringArray::uncover(const unsigned encode,
                            const unsigned oldLineIndex) {
  coverage.uncover(encode);
  unsigned coverCount = coverage.coverCount(encode);
  if (coverCount == 0) {
    uncoveredTuples.push(encode);
    oneCoveredTuples.pop(encode, oldLineIndex, coverage.getTuple(encode));
  }
  if (coverCount == 1) {
    const std::vector<unsigned> &tuple = coverage.getTuple(encode);
    const std::vector<unsigned> &columns = coverage.getColumns(encode);
    for (size_t lineIndex = 0; lineIndex < array.size(); ++lineIndex) {
      if (lineIndex == oldLineIndex) {
        continue;
      }
      auto &line = array[lineIndex];
      bool match = true;
      for (size_t i = 0; i < columns.size(); ++i) {
        if (tuple[i] != line[columns[i]]) {
          match = false;
          break;
        }
      }
      if (match) {
        oneCoveredTuples.push(encode, lineIndex, tuple);
        break;
      }
    }
  }
}

void CoveringArray::uncover_with_lock(const unsigned encode,
                                      const unsigned oldLineIndex) {
  coverage.uncover(encode);
  unsigned coverCount = coverage.coverCount(encode);
  if (coverCount == 0) {
    {
      std::lock_guard<std::mutex> lg(uncoveredTuplesMutex);
      uncoveredTuples.push(encode);
    }
    {
      std::lock_guard<std::mutex> lg(oneCoveredTuplesMutex);
      oneCoveredTuples.pop(encode, oldLineIndex, coverage.getTuple(encode));
    }
  }
  if (coverCount == 1) {
    const std::vector<unsigned> &tuple = coverage.getTuple(encode);
    const std::vector<unsigned> &columns = coverage.getColumns(encode);
    for (size_t lineIndex = 0; lineIndex < array.size(); ++lineIndex) {
      if (lineIndex == oldLineIndex) {
        continue;
      }
      auto &line = array[lineIndex];
      bool match = true;
      for (size_t i = 0; i < columns.size(); ++i) {
        if (tuple[i] != line[columns[i]]) {
          match = false;
          break;
        }
      }
      if (match) {
        {
          std::lock_guard<std::mutex> lg(oneCoveredTuplesMutex);
          oneCoveredTuples.push(encode, lineIndex, tuple);
        }
        break;
      }
    }
  }
}

void CoveringArray::tmpPrint() {
  struct timeval end_time;
  gettimeofday(&end_time, NULL);
  double start_time_sec = start_time.tv_sec + start_time.tv_usec / 1000000.0;
  double end_time_sec = end_time.tv_sec + end_time.tv_usec / 1000000.0;
  std::cout << end_time_sec - start_time_sec << '\t' << array.size() << '\t'
            << step << std::endl;
}

bool CoveringArray::verify(
    const std::vector<std::vector<unsigned>> &resultArray) {
  const unsigned strength = specificationFile.getStrenth();
  const Options &options = specificationFile.getOptions();
  Coverage tmpCoverage(specificationFile);
  tmpCoverage.initialize(satSolver, option_constrained_indicator);

  std::vector<unsigned> tuple(strength);
  unsigned lineIndex = 0;
  for (auto &line : resultArray) {
    for (unsigned column = 0; column < line.size(); ++column) {
      if (line[column] < options.firstSymbol(column) ||
          line[column] > options.lastSymbol(column)) {
        std::cerr << "error line: " << lineIndex;
        std::cerr << " option: " << column << std::endl;
        std::cerr << "should be " << options.firstSymbol(column)
                  << " <= var <= " << options.lastSymbol(column) << std::endl;
        for (auto var : line) {
          std::cerr << var << ' ';
        }
        std::cerr << std::endl;
        return false;
      }
    }
    for (std::vector<unsigned> columns = combinadic.begin(strength);
         columns[strength - 1] < line.size(); combinadic.next(columns)) {
      for (unsigned i = 0; i < strength; ++i) {
        tuple[i] = line[columns[i]];
      }
      unsigned encode = tmpCoverage.encode(columns, tuple);
      if (tmpCoverage.coverCount(encode) < 0) {
        std::cerr << "violate constraints" << std::endl;
        return false;
      }
      tmpCoverage.cover(encode);
    }
    ++lineIndex;
  }
  return tmpCoverage.allIsCovered();
}
