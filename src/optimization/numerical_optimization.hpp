#ifndef NUMERICAL_OPTIMIZATION_HPP_INCLUDED
#define NUMERICAL_OPTIMIZATION_HPP_INCLUDED

#include "../optimization/optimization.hpp"
#include "../parameters/parameters.hpp"
#include "../utilities/utilities.hpp"
namespace boat{
class InstantiableHolder {
 public:
  virtual void assign(const std::vector<double>& input) const = 0;
  virtual std::vector<double> value() const = 0;
  virtual size_t size() const = 0;
  virtual const std::vector<double>& get_lower_bound() const = 0;
  virtual const std::vector<double>& get_upper_bound() const = 0;
  virtual std::unique_ptr<InstantiableHolder> clone() const = 0;
  virtual ~InstantiableHolder(){};
};

template <class... Args>
class ParameterHolder : public InstantiableHolder {
 public:
  ParameterHolder(RangeParameter<Args>&... args)
      : instantiables_(args...),
        lower_bound_(lower_bound_helper(gen_seq<sizeof...(Args)>())),
        upper_bound_(upper_bound_helper(gen_seq<sizeof...(Args)>())) {}

  void assign(const std::vector<double>& input) const override {
    auto s = gen_seq<sizeof...(Args)>();
    instantitate_helper(input, s);
  }

  std::vector<double> value() const override {
    auto s = gen_seq<sizeof...(Args)>();
    return value_helper(s);
  }

  size_t size() const override { return std::tuple_size<tuple_type>::value; }

  const std::vector<double>& get_lower_bound() const override {
    return lower_bound_;
  }

  const std::vector<double>& get_upper_bound() const override {
    return upper_bound_;
  }

  std::unique_ptr<InstantiableHolder> clone() const override {
    return std::make_unique<ParameterHolder>(*this);
  }

  ~ParameterHolder() override{};

 private:
  template <int... S>
  std::vector<double> lower_bound_helper(seq<S...>) const {
    std::vector<double> lb(std::tuple_size<tuple_type>::value);
    auto l = {(lb[S] = std::get<S>(instantiables_).get_lower(), 0)...};
    (void)l;  // Avoids unused var warning;
    return lb;
  }

  template <int... S>
  std::vector<double> upper_bound_helper(seq<S...>) const {
    std::vector<double> ub(std::tuple_size<tuple_type>::value);
    auto l = {(ub[S] = std::get<S>(instantiables_).get_upper(), 0)...};
    (void)l;  // Avoids unused var warning;
    return ub;
  }

  template <int... S>
  void instantitate_helper(const std::vector<double>& input, seq<S...>) const {
    assert(input.size() == std::tuple_size<tuple_type>::value);
    auto l = {(std::get<S>(instantiables_).instantiate_in_range(input[S]), 0)...};
    (void)l;  // Avoids unused var warning;
  }

  template <int... S>
  std::vector<double> value_helper(seq<S...>) const {
    std::vector<double> val(std::tuple_size<tuple_type>::value);
    auto l = {(val[S] = std::get<S>(instantiables_).value(), 0)...};
    (void)l;  // Avoids unused var warning;
    return val;
  }


  typedef std::tuple<RangeParameter<Args>&...> tuple_type;
  tuple_type instantiables_;
  std::vector<double> lower_bound_;
  std::vector<double> upper_bound_;
};

class MergedHolder : public InstantiableHolder {
 public:
  typedef std::vector<std::unique_ptr<InstantiableHolder>> Input;
  MergedHolder(Input holders)
      : holders_(std::move(holders)) {
  //Could check if some of them are MergedHolders and flatten
    for (auto& holder : holders_) {
      lower_bound_.insert(lower_bound_.end(), holder->get_lower_bound().begin(),
                          holder->get_lower_bound().end());
      upper_bound_.insert(upper_bound_.end(), holder->get_upper_bound().begin(),
                          holder->get_upper_bound().end());
    }
  }
  void assign(const std::vector<double>& input) const override {
    size_t ind = 0;
    for (auto& holder : holders_) {
      std::vector<double> tmp(input.begin() + ind,
                              input.begin() + ind + holder->size());
      holder->assign(tmp);
      ind += holder->size();
    }
  }

  std::vector<double> value() const override {
    std::vector<double> res;
    for (auto& holder : holders_) {
      auto tmp = holder->value();
      res.insert(res.end(), tmp.begin(), tmp.end());
    }
    return res;
  }

  size_t size() const override { return lower_bound_.size(); }

  const std::vector<double>& get_lower_bound() const override {
    return lower_bound_;
  }

  const std::vector<double>& get_upper_bound() const override {
    return upper_bound_;
  }

  std::unique_ptr<InstantiableHolder> clone() const override {
    auto res = std::unique_ptr<MergedHolder>(new MergedHolder());
    res->lower_bound_ = lower_bound_;
    res->upper_bound_ = upper_bound_;
    for(auto& ptr : holders_){
      res->holders_.push_back(ptr->clone());
    }
    return res;
  }

 private:

  MergedHolder() = default;

  std::vector<std::unique_ptr<InstantiableHolder>> holders_;
  std::vector<double> lower_bound_;
  std::vector<double> upper_bound_;
};

class EmptyHolder : public InstantiableHolder {
  void assign(const std::vector<double>& input) const override {
    assert(input.empty());
  }
  size_t size() const override {
    return 0;
  }
  const std::vector<double>& get_lower_bound() const override {
    return dummy_;
  }
  const std::vector<double>& get_upper_bound() const override {
    return dummy_;
  }

  std::unique_ptr<InstantiableHolder> clone() const override {
    return std::make_unique<EmptyHolder>();
  }

  std::vector<double> value() const override {
    return std::vector<double>();
  }

  static const std::vector<double> dummy_;
};

template <class InstantiationFunc>
std::unique_ptr<InstantiableHolder> make_holder(
    InstantiationFunc func, std::vector<double> lower_bound,
    std::vector<double> upper_bound) {
  struct FunctionalHolder : public InstantiableHolder {
    FunctionalHolder(InstantiationFunc f, std::vector<double> lb,
                     std::vector<double> ub)
        : func_(std::move(f)),
          lower_bound_(std::move(lb)),
          upper_bound_(std::move(ub)) {
      assert(lower_bound_.size() == upper_bound_.size());
    }

    void assign(const std::vector<double>& input) const override {
      func_(input);
    }


    std::vector<double> value() const override {
      //TODO:: provide a class that also has a value() function
      assert(false && "Can not use value with a FunctionalHolder");
      return std::vector<double>();
    }

    size_t size() const override { return lower_bound_.size(); }

    const std::vector<double>& get_lower_bound() const override {
      return lower_bound_;
    }
    const std::vector<double>& get_upper_bound() const override {
      return upper_bound_;
    }
    std::unique_ptr<InstantiableHolder> clone() const override {
      return std::make_unique<FunctionalHolder>(*this);
    }
    InstantiationFunc func_;
    std::vector<double> lower_bound_;
    std::vector<double> upper_bound_;
  };

  return std::make_unique<FunctionalHolder>(
      std::move(func), std::move(lower_bound), std::move(upper_bound));
}

template <class... Args>
std::unique_ptr<InstantiableHolder> make_holder(RangeParameter<Args>&... args) {
  return std::make_unique<ParameterHolder<Args...>>(args...);
}

template <class ObjectiveResult>
class NumericalOpt : public OptimizationBase {
 public:
  NumericalOpt(std::shared_ptr<InstantiableHolder> instantiables);
  NumericalOpt& set_initial_guess(std::vector<double> initial_guess);
  NumericalOpt& set_utility_permanent(bool is_permanent);

  void set_objective_function(std::function<ObjectiveResult()> func);
  void set_objective_function(
      std::function<ObjectiveResult(const std::vector<double>&)> func);

  ObjectiveResult run_optimization();

 protected:
  double iteration_internal(const std::vector<double>& input);
  virtual void run_optimization_internal() = 0;
  std::shared_ptr<InstantiableHolder> instantiables_;
  std::vector<double> initial_guess_;

 private:
  std::function<ObjectiveResult(const std::vector<double>&)> utility_function_;
  std::vector<IterationResult<ObjectiveResult>> all_results_;
  bool utility_permanent_;
  ParameterInstance* optimization_instance_;
};

template <class ObjectiveResult>
NumericalOpt<ObjectiveResult>::NumericalOpt(
    std::shared_ptr<InstantiableHolder> instantiables)
    : instantiables_(std::move(instantiables)),
      utility_permanent_(true),
      optimization_instance_(nullptr) {}

template <class ObjectiveResult>
NumericalOpt<ObjectiveResult>& NumericalOpt<ObjectiveResult>::set_initial_guess(
    std::vector<double> initial_guess) {
  initial_guess_ = std::move(initial_guess);
  return *this;
}

template <class ObjectiveResult>
NumericalOpt<ObjectiveResult>& NumericalOpt<ObjectiveResult>::set_utility_permanent(
    bool utility_permanent) {
  utility_permanent_ = utility_permanent;
  return *this;
}

template <class ObjectiveResult>
void NumericalOpt<ObjectiveResult>::set_objective_function(
    std::function<ObjectiveResult()> func) {
  utility_function_ = [func{std::move(func)}](const std::vector<double>&) {
    return func();
  };
}

template <class ObjectiveResult>
void NumericalOpt<ObjectiveResult>::set_objective_function(
    std::function<ObjectiveResult(const std::vector<double>&)> func) {
  utility_function_ = std::move(func);
}

template <class ObjectiveResult>
ObjectiveResult NumericalOpt<ObjectiveResult>::run_optimization() {
  assert(optimization_instance_ == nullptr);
  optimization_instance_ = ParameterInstance::current_instance;
  // Set the initial guess if it hasn't been done

  if (initial_guess_.empty()) {
    for (size_t i = 0; i < instantiables_->size(); i++) {
      initial_guess_.push_back((instantiables_->get_upper_bound()[i] +
                                instantiables_->get_lower_bound()[i]) /
                               2.0);
    }
  }
  assert(initial_guess_.size() == instantiables_->size());

  run_optimization_internal();
  optimization_instance_->promote(all_results_[best_iteration()].instance_);
  return std::move(all_results_[best_iteration()].utility_result_);
}

template <class ObjectiveResult>
double NumericalOpt<ObjectiveResult>::iteration_internal(
    const std::vector<double>& input) {
  //Sanity check
  for(double d : input) {
    assert(!std::isnan(d));
  }
  assert(utility_function_);
  ParameterInstance* instance = optimization_instance_->create_child();
  instantiables_->assign(input);
  if (!utility_permanent_) {
    instance->create_child();
  }  // Tempory instance
  all_results_.emplace_back(instance, utility_function_(input));
  if (!utility_permanent_) {
    instance->remove_children();
  }
  add_scalar_result(all_results_.back().utility());
  assert(!std::isnan(all_results_.back().utility()));
  return all_results_.back().utility();
}
}
#endif  // NUMERICAL_OPTIMIZATION_HPP_INCLUDED
