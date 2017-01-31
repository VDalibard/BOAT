#ifndef OPTIMIZATION_HPP_INCLUDED
#define OPTIMIZATION_HPP_INCLUDED

#include <vector>
#include <memory>
#include <assert.h>
#include <random>
#include "../parameters/parameters.hpp"
namespace boat{
template <class ObjectiveResult = double, class... SuboptResult>
class BayesOpt;

template <class ItResult = double>
class SimpleOpt;

template <class ObjectiveResult = double, class... SuboptResult>
struct IterationResult;

// C++ forces us to write all classes specialization
template <>
struct IterationResult<double> {
  ParameterInstance* instance_;
  double utility_result_;

  IterationResult(ParameterInstance* instance, double utility_result);

  double utility() const;
};

template <>
struct IterationResult<std::unordered_map<std::string, double>> {
  ParameterInstance* instance_;
  std::unordered_map<std::string, double> utility_result_;

  IterationResult(ParameterInstance* instance,
                  std::unordered_map<std::string, double> utility_result);

  double utility() const;
};

template <class ObjectiveResult>
struct IterationResult<ObjectiveResult> {
  ParameterInstance* instance_;
  ObjectiveResult utility_result_;
  double utility_value_;

  IterationResult(ParameterInstance* instance, ObjectiveResult utility_result);

  double utility() const;
};

template <class SuboptResult>
struct IterationResult<double, SuboptResult> {
  ParameterInstance* instance_;
  SuboptResult instantiation_result_;
  double utility_result_;

  IterationResult(ParameterInstance* instance,
                  SuboptResult instantitation_result,
                  double utility_result);

  double utility() const;
};

template <class ObjectiveResult, class SuboptResult>
struct IterationResult<ObjectiveResult, SuboptResult> {
  ParameterInstance* instance_;
  SuboptResult instantiation_result_;
  ObjectiveResult utility_result_;
  double utility_value_;

  IterationResult(ParameterInstance* instance,
                  SuboptResult instantitation_result,
                  ObjectiveResult utility_result);

  double utility() const;
};

class OptimizationBase {
  template <class T, class... U>
  friend class Optimization;

 public:
  OptimizationBase();
  OptimizationBase& set_minimizing();
  OptimizationBase& set_maximizing();
  OptimizationBase& set_max_num_iterations(size_t num_iterations);
  OptimizationBase& set_stopping_relative_rate(size_t num_iterations,
                                               double ratio);
  OptimizationBase& set_stopping_criterion_function(
      std::function<bool()> stopping_criterion_function);

  bool minimizing() const;
  double stopping_relative_rate() const;
  size_t max_num_iterations() const;
  bool has_more_iterations() const;
  size_t iteration_count() const;
  double best_objective() const;
  size_t best_iteration() const;
  double last_objective() const;

  void add_scalar_result(double res);

 private:
  // Stopping conditions
  bool minimizing_;
  size_t max_num_iterations_;
  size_t stoprel_num_iters_;
  double stoprel_ratio_;

  // Results so far
  size_t best_so_far_;
  double best_objective_so_far_;
  double best_objective_delayed_;
  std::vector<double> scalar_results_;
  std::function<bool()> stopping_criterion_function_;
};

template <class ObjectiveResult, class... SuboptResult>
class Optimization : public OptimizationBase {
 public:
  typedef IterationResult<ObjectiveResult, SuboptResult...>
      iteration_result_type;

  ObjectiveResult run_optimization();
  void add_result(iteration_result_type ir);

  const std::vector<iteration_result_type>& get_all_results();
  virtual void iteration_internal(ParameterInstance* instance) = 0;

 private:
  std::vector<iteration_result_type> all_results_;
};

template <class ItResult>
class SimpleOpt : public Optimization<ItResult> {
 public:
  void set_objective_function(std::function<ItResult()> func);
  void set_objective_function(std::function<ItResult(int)> func);

 private:
  void iteration_internal(ParameterInstance* instance) override;
  std::function<ItResult(int)> iteration_function_;
};

template <class ObjectiveResult, class... SuboptResult>
class BayesOpt;

template <class ObjectiveResult>
class BayesOpt<ObjectiveResult> : public Optimization<ObjectiveResult> {
 public:
  BayesOpt();

  void set_subopt_function(std::function<void()> func);

  void set_objective_function(std::function<ObjectiveResult()> func);

  void set_learning_function(std::function<
      void(const std::vector<IterationResult<ObjectiveResult> >&)> func);
  void set_learning_function(
      std::function<void(const IterationResult<ObjectiveResult>&)> func);
  void set_learning_function(std::function<void(const ObjectiveResult&)> func);

  void set_utility_permanent(bool is_permanent);
  void set_learn_from_last(bool learn);

 private:
  std::function<void()> instantiation_function_;
  std::function<ObjectiveResult()> utility_function_;
  std::function<void(const std::vector<IterationResult<ObjectiveResult> >&)>
      learning_function_;

  void iteration_internal(ParameterInstance* instance) override;
  bool utility_permanent_;
  bool learn_from_last_;
};

template <class ObjectiveResult, class SuboptResult>
class BayesOpt<ObjectiveResult, SuboptResult>
    : public Optimization<ObjectiveResult, SuboptResult> {
 public:
  BayesOpt();

  void set_subopt_function(std::function<SuboptResult()> func);

  void set_objective_function(
      std::function<ObjectiveResult(const SuboptResult&)> func);
  void set_objective_function(std::function<ObjectiveResult()> func);

  void set_learning_function(std::function<void(const std::vector<
      IterationResult<ObjectiveResult, SuboptResult> >&)> func);
  void set_learning_function(std::function<
      void(const IterationResult<ObjectiveResult, SuboptResult>&)> func);
  void set_learning_function(std::function<
      void(const SuboptResult&, const ObjectiveResult&)> func);

  void set_utility_permanent(bool is_permanent);
  void set_learn_from_last(bool learn);

 private:
  std::function<SuboptResult()> instantiation_function_;
  std::function<ObjectiveResult(const SuboptResult&)> utility_function_;
  std::function<void(
      const std::vector<IterationResult<ObjectiveResult, SuboptResult> >&)>
      learning_function_;

  void iteration_internal(ParameterInstance* instance) override;
  bool utility_permanent_;
  bool learn_from_last_;
};


template <class ItResult>
class SimulatedAnnealingOpt : public Optimization<ItResult> {
 public:
  SimulatedAnnealingOpt();
  void set_generation_function(std::function<void()>);
  void set_mutation_function(std::function<void(double, ParameterInstance*)>);
  void set_objective_function(std::function<ItResult()>);

  void set_cooling_rate(double cooling_rate);
  void set_initial_temperature(double initial_temperature);

 private:
  void iteration_internal(ParameterInstance* instance) override;
  void test_new_result(double utility);

  std::function<void(double, ParameterInstance*)> mutation_function_;
  std::function<void()> generation_function_;
  std::function<ItResult()> utility_function_;

  double previous_utility_;
  double initial_temperature_;
  double temperature_;
  double cooling_rate_;
  ParameterInstance* previous_instance_;

  // TODO: Add utility permanent
};


// Useful for optimizations
double expected_improvement(double mean, double var, double best, bool minimizing = true);

template <class T>
void uniform_random_instantiate(RangeParameter<T>& param);
}
#include "optimization_t.hpp"

#endif  // OPTIMIZATION_HPP_INCLUDED
