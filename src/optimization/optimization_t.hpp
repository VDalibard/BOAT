#ifndef OPTIMIZATION_T_HPP_INCLUDED
#define OPTIMIZATION_T_HPP_INCLUDED
#include <algorithm>
#include "../utilities/useful_macros.hpp"
// Iteration Result
// void, double
// In optimization.cpp
namespace boat{
// void, ObjectiveResult
template <class ObjectiveResult>
IterationResult<ObjectiveResult>::IterationResult(ParameterInstance* instance,
                                                ObjectiveResult utility_result)
    : instance_(instance),
      utility_result_(std::move(utility_result)),
      utility_value_(utility_result_.utility()) {}

template <class ObjectiveResult>
double IterationResult<ObjectiveResult>::utility() const {
  return utility_value_;
}

// InstantiationResult, double
template <class InstantiationResult>
IterationResult<double, InstantiationResult>::IterationResult(
    ParameterInstance* instance, InstantiationResult instantitation_result,
    double utility_result)
    : instance_(instance),
      instantiation_result_(std::move(instantitation_result)),
      utility_result_(utility_result) {}

template <class InstantiationResult>
double IterationResult<double, InstantiationResult>::utility() const {
  return utility_result_;
}

// InstantiationResult, ObjectiveResult
template <class ObjectiveResult, class InstantiationResult>
IterationResult<ObjectiveResult, InstantiationResult>::IterationResult(
    ParameterInstance* instance, InstantiationResult instantitation_result,
    ObjectiveResult utility_result)
    : instance_(instance),
      instantiation_result_(std::move(instantitation_result)),
      utility_result_(std::move(utility_result)),
      utility_value_(utility_result_.utility()) {}

template <class ObjectiveResult, class InstantiationResult>
double IterationResult<ObjectiveResult, InstantiationResult>::utility() const {
  return utility_value_;
}

// Optimization
template <class ObjectiveResult, class... InstantiationResult>
ObjectiveResult
Optimization<ObjectiveResult, InstantiationResult...>::run_optimization() {
  assert(ParameterInstance::current_instance != nullptr);
  ParameterInstance* optimization_instance =
      ParameterInstance::current_instance;
  while (has_more_iterations()) {
    ParameterInstance* instance = optimization_instance->create_child();
    iteration_internal(instance);
  }
  optimization_instance->promote(all_results_[best_so_far_].instance_);
  return std::move(all_results_[best_so_far_].utility_result_);
}

template <class ObjectiveResult, class... InstantiationResult>
void Optimization<ObjectiveResult, InstantiationResult...>::add_result(
    IterationResult<ObjectiveResult, InstantiationResult...> ir) {
  all_results_.push_back(std::move(ir));
  add_scalar_result(all_results_.back().utility());
  // Check how we are doing with best so far and best delayed
}

template <class ObjectiveResult, class... InstantiationResult>
const std::vector<IterationResult<ObjectiveResult, InstantiationResult...> >&
Optimization<ObjectiveResult, InstantiationResult...>::get_all_results() {
  return all_results_;
}

// SimpleOpt

template <class ItResult>
void SimpleOpt<ItResult>::set_objective_function(
    std::function<ItResult(int)> func) {
  iteration_function_ = std::move(func);
}

template <class ItResult>
void SimpleOpt<ItResult>::set_objective_function(
    std::function<ItResult()> func) {
  iteration_function_ = [func{std::move(func)}](int dummy) { return func(); };
}

template <class ItResult>
void SimpleOpt<ItResult>::iteration_internal(ParameterInstance* instance) {
  this->add_result(IterationResult<ItResult>(
      instance, iteration_function_(this->iteration_count())));
}

// BayesOpt
// No InstantiationResult specialization

template <class ObjectiveResult>
BayesOpt<ObjectiveResult>::BayesOpt()
    : utility_permanent_(true), learn_from_last_(true) {}

template <class ObjectiveResult>
void BayesOpt<ObjectiveResult>::set_subopt_function(
    std::function<void()> func) {
  instantiation_function_ = std::move(func);
}

template <class ObjectiveResult>
void BayesOpt<ObjectiveResult>::set_objective_function(
    std::function<ObjectiveResult()> func) {
  utility_function_ = std::move(func);
}

template <class ObjectiveResult>
void BayesOpt<ObjectiveResult>::set_learning_function(std::function<
    void(const std::vector<IterationResult<ObjectiveResult> >&)> func) {
  learning_function_ = std::move(func);
}

template <class ObjectiveResult>
void BayesOpt<ObjectiveResult>::set_learning_function(
    std::function<void(const IterationResult<ObjectiveResult>&)> func) {
  learning_function_ = [func{std::move(func)}](
      const std::vector<IterationResult<ObjectiveResult> >& all_results) {
    func(all_results.back());
  };
}

template <class ObjectiveResult>
void BayesOpt<ObjectiveResult>::set_learning_function(
    std::function<void(const ObjectiveResult&)> func) {
  learning_function_ = [func{std::move(func)}](
      const std::vector<IterationResult<ObjectiveResult> >& all_results) {
    func(all_results.back().utility_result_);
  };
}

template <class ObjectiveResult>
void BayesOpt<ObjectiveResult>::set_utility_permanent(bool utility_permanent) {
  utility_permanent_ = utility_permanent;
}

template <class ObjectiveResult>
void BayesOpt<ObjectiveResult>::iteration_internal(ParameterInstance* instance) {
  // Perform the iteration
  instantiation_function_();

  // Sometimes, in order to evaluate the utility, we may wish to assign
  // some other parameters
  // We create a new parameter instance so that these instantiations
  // are not visible once this optimization returns
  if (!utility_permanent_) {
    instance->create_child();
  }  // Tempory instance
  ObjectiveResult utility_result = utility_function_();
  if (!utility_permanent_) {
    instance->remove_children();
  }

  this->add_result(
      IterationResult<ObjectiveResult>(instance, std::move(utility_result)));
  // If this isn't the last iteration, we learn

  if (learn_from_last_ || this->has_more_iterations()) {
    learning_function_(this->get_all_results());
  }
}

//--------------

template <class ObjectiveResult, class InstantiationResult>
BayesOpt<ObjectiveResult, InstantiationResult>::BayesOpt()
    : utility_permanent_(true) {}

template <class ObjectiveResult, class InstantiationResult>
void BayesOpt<ObjectiveResult, InstantiationResult>::set_subopt_function(
    std::function<InstantiationResult()> func) {
  instantiation_function_ = std::move(func);
}

template <class ObjectiveResult, class InstantiationResult>
void BayesOpt<ObjectiveResult, InstantiationResult>::set_objective_function(
    std::function<ObjectiveResult(const InstantiationResult&)> func) {
  utility_function_ = std::move(func);
}

template <class ObjectiveResult, class InstantiationResult>
void BayesOpt<ObjectiveResult, InstantiationResult>::set_objective_function(
    std::function<ObjectiveResult()> func) {
  utility_function_ = [func{std::move(func)}](
      const InstantiationResult& inst_res) {
    return func();
  };
}

template <class ObjectiveResult, class InstantiationResult>
void BayesOpt<ObjectiveResult, InstantiationResult>::set_learning_function(
    std::function<void(const std::vector<
        IterationResult<ObjectiveResult, InstantiationResult> >&)> func) {
  learning_function_ = std::move(func);
}

template <class ObjectiveResult, class InstantiationResult>
void BayesOpt<ObjectiveResult, InstantiationResult>::set_learning_function(
    std::function<void(
        const IterationResult<ObjectiveResult, InstantiationResult>&)> func) {
  learning_function_ = [func{std::move(func)}](
      const std::vector<IterationResult<ObjectiveResult> >& all_results) {
    func(all_results.back());
  };
}

template <class ObjectiveResult, class InstantiationResult>
void BayesOpt<ObjectiveResult, InstantiationResult>::set_learning_function(
    std::function<void(const InstantiationResult&, const ObjectiveResult&)>
        func) {
  learning_function_ = [func{std::move(func)}](const std::vector<
      IterationResult<ObjectiveResult, InstantiationResult> >& all_results) {
    func(all_results.back().instantiation_result_,
         all_results.back().utility_result_);
  };
}

template <class ObjectiveResult, class InstantiationResult>
void BayesOpt<ObjectiveResult, InstantiationResult>::set_utility_permanent(
    bool utility_permanent) {
  utility_permanent_ = utility_permanent;
}

template <class ObjectiveResult, class InstantiationResult>
void BayesOpt<ObjectiveResult, InstantiationResult>::set_learn_from_last(
    bool learn) {
  learn_from_last_ = learn;
}

template <class ObjectiveResult, class InstantiationResult>
void BayesOpt<ObjectiveResult, InstantiationResult>::iteration_internal(
    ParameterInstance* instance) {
  // Perform the iteration
  InstantiationResult instantiation_result = instantiation_function_();

  // Sometimes, in order to evaluate the utility, we may wish to assign
  // some other parameters
  // We create a new parameter instance so that these instantiations
  // are not visible once this optimization returns
  if (!utility_permanent_) {
    instance->create_child();
  }  // Tempory instance
  ObjectiveResult utility_result = utility_function_(instantiation_result);
  if (!utility_permanent_) {
    instance->remove_children();
  }

  this->add_result(IterationResult<ObjectiveResult, InstantiationResult>(
      instance, std::move(instantiation_result), std::move(utility_result)));
  // If this isn't the last iteration, we learn

  if (learn_from_last_ || this->has_more_iterations()) {
    learning_function_(this->get_all_results());
  }
}

//-----

template <class ItResult>
SimulatedAnnealingOpt<ItResult>::SimulatedAnnealingOpt()
    : initial_temperature_(1.0),
      temperature_(initial_temperature_),
      cooling_rate_(0.95) {}

template <class ItResult>
void SimulatedAnnealingOpt<ItResult>::set_generation_function(
    std::function<void()> generation_func) {
  generation_function_ = std::move(generation_func);
}

template <class ItResult>
void SimulatedAnnealingOpt<ItResult>::set_mutation_function(
    std::function<void(double, ParameterInstance*)> mutation_func) {
  mutation_function_ = std::move(mutation_func);
}

template <class ItResult>
void SimulatedAnnealingOpt<ItResult>::set_objective_function(
    std::function<ItResult()> utility_func) {
  utility_function_ = std::move(utility_func);
}

template <class ItResult>
void SimulatedAnnealingOpt<ItResult>::set_cooling_rate(double cooling_rate) {
  assert(0.0 < cooling_rate && cooling_rate <= 1.0);
  cooling_rate_ = cooling_rate;
}

template <class ItResult>
void SimulatedAnnealingOpt<ItResult>::set_initial_temperature(
    double initial_temperature) {
  assert(initial_temperature > 0.0);
  initial_temperature_ = initial_temperature;
}

template <class ItResult>
void SimulatedAnnealingOpt<ItResult>::iteration_internal(
    ParameterInstance* instance) {
  if (this->iteration_count() == 0) {
    temperature_ = initial_temperature_;
    generation_function_();
    this->add_result(IterationResult<ItResult>(instance, utility_function_()));
    previous_utility_ = this->last_objective();
    previous_instance_ = ParameterInstance::current_instance;
  } else {
    mutation_function_(temperature_, previous_instance_);
    this->add_result(IterationResult<ItResult>(instance, utility_function_()));
    test_new_result(this->last_objective());
  }
  temperature_ = temperature_ * cooling_rate_;
}
extern std::mt19937 generator;  // in probabilistic/probabilistic.cpp
template <class ItResult>
void SimulatedAnnealingOpt<ItResult>::test_new_result(double utility) {
  bool accepted = (utility > previous_utility_) ^ this->minimizing();
  if (!accepted) {
    double prob =
        1.0 /
        (1.0 + std::exp(std::abs(utility - previous_utility_) / temperature_));
    accepted = prob > std::uniform_real_distribution<>(0.0, 1.0)(generator);
  }
  if (accepted) {
    previous_utility_ = utility;
    previous_instance_ = ParameterInstance::current_instance;
  }
}


template <class T>
void uniform_random_instantiate(RangeParameter<T>& param){
  param.assign(
    std::uniform_real_distribution<double>(param.get_lower(),param.get_upper())
                                           (generator));
}
}

#endif  // OPTIMIZATION_T_HPP_INCLUDED
