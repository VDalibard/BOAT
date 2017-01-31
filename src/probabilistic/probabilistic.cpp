#include "probabilistic.hpp"
#include <algorithm>
// Static data

namespace boat {
std::mt19937 generator;
Averager* Averager::current_averager = nullptr;
std::string dag_model_fake_name = "dontusethisasname1234";

double log_sum_exp(const std::vector<double>& log_values) {
  double max_log_value = log_values[0];
  for (const auto& log_value : log_values) {
    if (max_log_value < log_value) {
      max_log_value = log_value;
    }
  }
  double sum = 0.0;
  for (const auto& log_value : log_values) {
    sum += exp(log_value - max_log_value);
  }
  return max_log_value + log(sum);
}

double get_normalized_likelihoods(std::vector<double>& normalized_likelihoods,
                                  std::vector<double>& log_likelihoods) {
  double normalization = log_sum_exp(log_likelihoods);
  normalized_likelihoods.clear();
  for (auto& log_likelihood : log_likelihoods) {
    log_likelihood -= normalization;
    normalized_likelihoods.push_back(exp(log_likelihood));
  }
  return normalization;
}

void draw_bins_rd(std::vector<size_t>& bins,
               const std::vector<double>& normalized_likelihoods) {
  bins.resize(normalized_likelihoods.size());
  for (auto& bin : bins) {
    bin = 0;
  }
  std::vector<double> normalized_likelihoods_cumul(normalized_likelihoods);
  for (size_t i = 1; i < normalized_likelihoods.size(); i++) {
    normalized_likelihoods_cumul[i] += normalized_likelihoods_cumul[i - 1];
  }

  // sum may not exactly be 1.0
  std::uniform_real_distribution<double> distribution(
      0.0, normalized_likelihoods_cumul.back());

  for (size_t b = 0; b < bins.size(); b++) {
    double draw = distribution(generator);
    auto it = std::upper_bound(normalized_likelihoods_cumul.begin(),
                               normalized_likelihoods_cumul.end(), draw);
    size_t p = it - normalized_likelihoods_cumul.begin();
    bins[p]++;
  }
}

std::vector<int> draw_bins(
    const std::vector<double>& normalized_likelihoods,
    int num_draws) {
  std::vector<int> bins(normalized_likelihoods.size(), 0);
  double num_draws_inv = 1.0 / double(num_draws);
  //double u_i =
  //    std::uniform_real_distribution<double>(0.0, num_draws_inv)(generator);
  // This may break some of the inference properties, but it guarantees that
  // when we use this function for an average, we always use the same particles
  // Most likely fine though as particles are shuffled
  double u_i = num_draws_inv * 0.5;
  int index = 0;
  for (int i = 0; i < num_draws; i++) {
    while (u_i > normalized_likelihoods[index]) {
      u_i -= normalized_likelihoods[index];
      index++;
      assert(index < SIZE(normalized_likelihoods));
    }
    bins[index]++;
    u_i += num_draws_inv;
  }
  return bins;
}

double calc_ess(const std::vector<double>& normalized_likelihoods) {
  double sum = 0.0;
  for (const auto& normalized_likelihood : normalized_likelihoods) {
    sum += pow(normalized_likelihood, 2.0);
  }
  return 1.0 / sum;
}

double normal_lnp(double x, double mean, double variance) {
  assert(!std::isnan(variance));
  assert(variance>0.0);
  constexpr double log_two_pi = log(2.0 * M_PI);
  double xmms = pow(x - mean, 2.0);
  double res = -0.5 * (xmms / variance + log_two_pi + log(variance));
  return res;
}

double average_normalized(const std::vector<double>& values,
                          const std::vector<double>& nl){
  double res = 0.0;
  for(int i=0; i<SIZE(nl); i++) {
    if(nl[i] == 0.0){continue;}
    res += values[i] * nl[i];
  }
  return res;
}

std::vector<double> average_normalized(
    const std::vector<std::vector<double>>& values,
    const std::vector<double>& nl) {
  assert(nl.size() == values.size());
  std::vector<double> res;
  for(int i=0; i<SIZE(nl); i++) {
    if(nl[i] != 0.0){
      if(res.empty()) {
        res = std::vector<double>(values[i].size(), 0.0);
      }
      assert(values[i].size() == res.size());
      for(int j=0; j<SIZE(res); j++){
        res[j] += values[i][j] * nl[i];
      }
    }
  }
  assert(!res.empty());
  return res;
}
}
