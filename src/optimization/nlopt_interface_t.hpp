#ifndef NLOPT_INTERFACE_T_HPP_INCLUDED
#define NLOPT_INTERFACE_T_HPP_INCLUDED
namespace boat{
template <class ObjectiveResult>
template <class... Args>
NLOpt<ObjectiveResult>::NLOpt(RangeParameter<Args> &... args)
    : NLOpt(make_holder(args...)) {}

template <class ObjectiveResult>
NLOpt<ObjectiveResult>::NLOpt(std::shared_ptr<InstantiableHolder> instantiables)
    : NumericalOpt<ObjectiveResult>(std::move(instantiables)),
      algorithm_(nlopt::GN_DIRECT_L) {}

template <class ObjectiveResult>
void NLOpt<ObjectiveResult>::run_optimization_internal() {
  auto lam = [&](const std::vector<double> &input) {
    double u;
    // NLOpt catches exceptions, so it is easier to handle them here
    try {
      u = this->iteration_internal(input);
    } catch (const std::exception &exc) {
      std::cerr << exc.what();
      assert(false);
    }
    if (!this->has_more_iterations()) {
      throw nlopt::forced_stop();
    }
    return u;
  };

  struct FuncPtrWrapper {
    static double wrapper(const std::vector<double> &x,
                          std::vector<double> &grad, void *my_func_data) {
      return (*((decltype(lam) *)my_func_data))(x);
    }
  };

  std::vector<double> lb(this->instantiables_->get_lower_bound());
  std::vector<double> ub(this->instantiables_->get_upper_bound());
  // NLopt sometimes goes outside bounds
  for (auto &e : lb) {
    e += 1e-8;
  }
  for (auto &e : ub) {
    e -= 1e-8;
  }

  nlopt::opt opt(algorithm_, this->instantiables_->size());
  opt.set_lower_bounds(lb);
  opt.set_upper_bounds(ub);

  if (this->minimizing()) {
    opt.set_min_objective(FuncPtrWrapper::wrapper, &lam);
  } else {
    opt.set_max_objective(FuncPtrWrapper::wrapper, &lam);
  }
  if (this->max_num_iterations() > 0) {
    opt.set_maxeval(this->max_num_iterations());
  } else {
    assert(this->stopping_relative_rate() > 0.0);
    // The optimization will terminate with a forced termination
    opt.set_maxeval(999999);  // Big number
  }
  try {
    double minf;
    opt.optimize(this->initial_guess_, minf);
  } catch (std::runtime_error &e) {
    // std::cerr << e.what();
  }
}

template <class ObjectiveResult>
NLOpt<ObjectiveResult> &NLOpt<ObjectiveResult>::set_algorithm(
    nlopt::algorithm algorithm) {
  algorithm_ = algorithm;
  return *this;
}
}
#endif  // NLOPT_INTERFACE_T_HPP_INCLUDED
