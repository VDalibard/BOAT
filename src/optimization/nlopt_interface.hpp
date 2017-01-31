#ifndef NLOPT_INTERFACE_HPP_INCLUDED
#define NLOPT_INTERFACE_HPP_INCLUDED

#include <nlopt.hpp>
#include "../optimization/numerical_optimization.hpp"
namespace boat{
template <class ObjectiveResult = double>
class NLOpt : public NumericalOpt<ObjectiveResult> {
 public:
  template <class... Args>
  NLOpt(RangeParameter<Args>&... args);
  NLOpt(std::shared_ptr<InstantiableHolder> instantiables);
  NLOpt& set_algorithm(nlopt::algorithm algorithm);

 private:
  void run_optimization_internal() override;
  nlopt::algorithm algorithm_;
};
}
#include "nlopt_interface_t.hpp"

#endif  // NLOPT_INTERFACE_HPP_INCLUDED
