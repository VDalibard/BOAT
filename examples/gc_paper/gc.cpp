
#include "boat.hpp"

using namespace std;

using boat::ProbEngine;
using boat::SemiParametricModel;
using boat::DAGModel;
using boat::GPParams;
using boat::NLOpt;
using boat::BayesOpt;
using boat::SimpleOpt;
using boat::generator;
using boat::RangeParameter;

struct GCRateModel : public SemiParametricModel<GCRateModel> {

  GCRateModel() {
    allocated_mbs_per_sec = std::uniform_real_distribution<>(0.0, 5000.0)(generator);

    p_.default_noise(0.0);
    p_.mean(uniform_real_distribution<>(0.0, 10.0)(generator));
    p_.stdev(uniform_real_distribution<>(0.0, 200.0)(generator));
    p_.linear_scales({uniform_real_distribution<>(0.0, 15.0)(generator)});
    set_params(p_);
  }

  double parametric(double eden_size) const {
    return allocated_mbs_per_sec / eden_size;
  }

  double allocated_mbs_per_sec;
  GPParams p_;
};


int main() {
  ProbEngine<GCRateModel> eng;  

  eng.observe(0.40, 1024);
  eng.observe(0.25, 2048);

  std::cout << eng.predict(1536) << std::endl;
}


