#include "Options.h"

void Options::initialize(const std::vector<unsigned> &value_counts) {
  cumulativeValueCounts.resize(value_counts.size());
  unsigned cumulativeValueCount = 0;
  for (std::vector<unsigned>::size_type i = 0; i < value_counts.size(); ++i) {
    cumulativeValueCount += value_counts[i];
    cumulativeValueCounts[i] = cumulativeValueCount;
  }
  owingOptions.resize(cumulativeValueCount);
  for (unsigned i = 0, j = 0; i < value_counts.size(); ++i) {
    for (unsigned k = 0; k < value_counts[i]; ++k) {
      owingOptions[j++] = i;
    }
  }
}

#ifndef NDEBUG
void Options::print() {
  std::cout << "********Debuging Options********" << std::endl;
  for (auto i : cumulativeValueCounts) {
    std::cout << i << ' ';
  }
  std::cout << std::endl;
  for (auto i : owingOptions) {
    std::cout << i << ' ';
  }
  std::cout << std::endl;
  std::cout << "********end of Debuging Options********" << std::endl;
}
#endif
