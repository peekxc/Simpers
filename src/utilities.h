#include <Rcpp.h>
using namespace Rcpp;

#include <functional>
#include <vector>
#include <algorithm>
using std::size_t;
using std::vector;

template <typename Lambda>
void apply_combinations(const size_t n, const size_t k, Lambda f){
  vector< bool > v(n);
  std::fill(v.begin(), v.begin() + k, true);
  vector< size_t > idx;
  size_t i, cc;
  do {
    i = 0, cc = 0;
    std::for_each(v.begin(), v.end(), [&i, &idx, &cc](bool val){
      if (val){ idx[cc++] = i++; }
    });
    f(idx); // gives 0-based indices to lambda
  } while (std::prev_permutation(v.begin(), v.end()));
}

List apply_combinations(const size_t n, const size_t k, Function f) {
  apply_combinations(n, k, [&f](const vector< size_t > idx){ f(idx); });
}

