//
// Created by hebing on 19-11-26.
//

#ifndef FASTCA_THREADTMPRESULT_H
#define FASTCA_THREADTMPRESULT_H

#include <limits>
#include <vector>
class ThreadTmpResult {
public:
  long long firstBestScore = std::numeric_limits<long long>::min();
  long long bestScore = std::numeric_limits<long long>::min();
  std::vector<unsigned> firstBestRows;
  std::vector<unsigned> firstBestVars;
  std::vector<unsigned> bestRows;
  std::vector<unsigned> bestVars;
};

class MediaLineVarTuple {
public:
  const unsigned encode;
  const unsigned oldLineIndex;
  const std::vector<unsigned> tuple;
  MediaLineVarTuple(const unsigned encode, const unsigned oldLineIndex,
                    const std::vector<unsigned> tuple)
      : encode(encode), oldLineIndex(oldLineIndex), tuple(tuple){};
};
class ReplaceTmpResult {
public:
  size_t popedUncoveredTuplesNum = 0;
  size_t popedOneCoveredTuplesNum = 0;
  std::vector<MediaLineVarTuple> pushedUncoveredTuples;
  std::vector<MediaLineVarTuple> pushedOneCoveredTuples;
};
#endif // FASTCA_THREADTMPRESULT_H
