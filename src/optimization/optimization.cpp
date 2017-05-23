#include <limits>
#include <numeric>
#include "optimization.hpp"
#include "numerical_optimization.hpp"
namespace boat{
const std::vector<double> EmptyHolder::dummy_;

IterationResult<double>::IterationResult(ParameterInstance* instance,
                                         double utility_result)
    : instance_(instance), utility_result_(utility_result) {}

double IterationResult<double>::utility() const { return utility_result_; }

IterationResult<std::unordered_map<std::string, double>>::IterationResult(
  ParameterInstance* instance,
  std::unordered_map<std::string, double> utility_result)
    : instance_(instance), utility_result_(std::move(utility_result)) {}

double IterationResult<std::unordered_map<std::string, double>>::utility() const {
   auto it = utility_result_.find("objective");
   assert(it != utility_result_.end());
   return it->second;
}


OptimizationBase::OptimizationBase()
    : minimizing_(true),
      max_num_iterations_(20),
      stoprel_num_iters_(999),
      stoprel_ratio_(0.0),
      best_so_far_(0),
      best_objective_so_far_(std::numeric_limits<double>::max()),
      best_objective_delayed_(std::numeric_limits<double>::max()) {}

OptimizationBase& OptimizationBase::set_minimizing() {
  assert(scalar_results_.empty());
  minimizing_ = true;
  best_objective_so_far_ = std::numeric_limits<double>::max();
  best_objective_delayed_ = std::numeric_limits<double>::max();
  return *this;
}

OptimizationBase& OptimizationBase::set_maximizing() {
  assert(scalar_results_.empty());
  minimizing_ = false;
  best_objective_so_far_ = std::numeric_limits<double>::lowest();
  best_objective_delayed_ = std::numeric_limits<double>::lowest();
  return *this;
}

OptimizationBase& OptimizationBase::set_max_num_iterations(
    size_t num_iterations) {
  max_num_iterations_ = num_iterations;
  return *this;
}

OptimizationBase& OptimizationBase::set_stopping_relative_rate(
    size_t num_iterations, double ratio) {
  stoprel_num_iters_ = num_iterations;
  stoprel_ratio_ = ratio;
  return *this;
}


OptimizationBase& OptimizationBase::set_stopping_criterion_function(
    std::function<bool()> stopping_criterion_function) {
  stopping_criterion_function_ = std::move(stopping_criterion_function);
  return *this;
}

bool OptimizationBase::minimizing() const { return minimizing_; }

size_t OptimizationBase::max_num_iterations() const {
  return max_num_iterations_;
}

double OptimizationBase::stopping_relative_rate() const {
  return stoprel_ratio_;
}

void OptimizationBase::add_scalar_result(double new_utility) {
  assert(!std::isnan(new_utility));
  scalar_results_.push_back(new_utility);
  bool new_best = minimizing_ ? new_utility < best_objective_so_far_
                              : new_utility > best_objective_so_far_;
  if (new_best) {
    best_so_far_ = scalar_results_.size() - 1;
    best_objective_so_far_ = new_utility;
  }
  if (scalar_results_.size() > stoprel_num_iters_) {
    // We will compare later the current best to the best utility we had up to
    // that point
    size_t to_consider = scalar_results_.size() - 1 - stoprel_num_iters_;
    double to_consider_utility = scalar_results_[to_consider];
    bool new_best_delayed = minimizing_
                                ? to_consider_utility < best_objective_delayed_
                                : to_consider_utility > best_objective_delayed_;
    if (new_best_delayed) {
      best_objective_delayed_ = to_consider_utility;
    }
  }
}

bool OptimizationBase::has_more_iterations() const {
  size_t ic = iteration_count();
  if (ic == 0) {
    return true;
  }
  if(stopping_criterion_function_ && stopping_criterion_function_()) {
    return false;
  }
  if (ic >= max_num_iterations_) {
    return false;
  }

  if (stoprel_ratio_ == 0.0) {
    return true;
  }
  if(scalar_results_.size() <= stoprel_num_iters_){
    return true;
  }
  bool stopping = std::abs(best_objective_so_far_ - best_objective_delayed_) <
                  stoprel_ratio_ * (std::abs(best_objective_so_far_) +
                                    std::abs(best_objective_delayed_)) *
                      0.5;
  return !stopping;
}

size_t OptimizationBase::iteration_count() const {
  return scalar_results_.size();
}

double OptimizationBase::best_objective() const {
  return best_objective_so_far_;
}

size_t OptimizationBase::best_iteration() const { return best_so_far_; }

double OptimizationBase::last_objective() const {
  return scalar_results_.back();
}


double ugaussian_pdf(double x){
  static constexpr double inv_sqrt_2pi = 1.0 / sqrt(2.0 * M_PI);
  return inv_sqrt_2pi * std::exp(-0.5 * x * x);
}

double ugaussian_cdf(double x) {
  static constexpr double inv_sqrt_2 = 1.0 / sqrt(2.0);
  return 0.5 * (1.0 + std::erf(x * inv_sqrt_2));
}

double expected_improvement(double mean, double var, double best,
                            bool minimizing) {
  assert(!std::isnan(mean));
  assert(!std::isnan(var));
  assert(!std::isnan(best));
  double gamma;
  if(minimizing){
    gamma = best - mean;
  } else {
    gamma = mean - best;
  }
  gamma /= sqrt(var);
  double sum = gamma * ugaussian_cdf(gamma) + ugaussian_pdf(gamma);
  double res = sum * sqrt(var);
  if(std::isnan(res)){
    // var is close to zero
    return 0.0;
  } else {
    return res;
  }
  
}

double exp_imp(double mean, double var, double best,
                            bool minimizing) {
  return expected_improvement(mean, var, best, minimizing);
}
}
