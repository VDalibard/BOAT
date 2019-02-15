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

struct BHParams{
  BHParams() : x1_(-5, 10), x2_(0, 15){}
  RangeParameter<double> x1_;
  RangeParameter<double> x2_;
};

/// The Branin-Hoo objective function
double branin_hoo(const BHParams& p){
  static constexpr double a = 1;
  static constexpr double b = 5.1 / (4 * pow(M_PI, 2));
  static constexpr double c = 5.0 / M_PI;
  static constexpr double r = 6.0;
  static constexpr double s = 10;
  static constexpr double t = 1.0 / (8 * M_PI);

  double x1 = p.x1_.value();
  double x2 = p.x2_.value();
  double fx = a * pow(x2 - b * pow(x1, 2) + c * x1 - r, 2) +
              s * (1.0 - t) * std::cos(x1) + s;
  return fx;
}

double branin_hoo_first_term(const BHParams& p){
  static constexpr double s = 10;
  static constexpr double t = 1.0 / (8 * M_PI);

  double x1 = p.x1_.value();
  double fx = s * (1.0 - t) * std::cos(x1) + s;
  return fx;
}


/// Naive way of optimizing it, model with simple GP prior
struct Param : public SemiParametricModel<Param> {
  Param() {
    p_.default_noise(0.0);
    p_.mean(uniform_real_distribution<>(0.0, 10.0)(generator));
    p_.stdev(uniform_real_distribution<>(0.0, 200.0)(generator));
    p_.linear_scales({uniform_real_distribution<>(0.0, 15.0)(generator),
                     uniform_real_distribution<>(0.0, 15.0)(generator)});
    set_params(p_);
  }

  GPParams p_;
};

struct FullModel : public DAGModel<FullModel> {
  FullModel(){
    eng_.set_num_particles(100);
  }
  void model(const BHParams& p) {
    output("objective", eng_, p.x1_.value(), p.x2_.value());
  }

  void print() {
    PR(AVG_PROP(eng_, p_.mean()));
    PR(AVG_PROP(eng_, p_.stdev()));
    PR(AVG_PROP(eng_, p_.linear_scales()[0]));
    PR(AVG_PROP(eng_, p_.linear_scales()[1]));
  }
  ProbEngine<Param> eng_;
};

void maximize_ei(FullModel& m, BHParams& p, double incumbent) {
  NLOpt<> opt(p.x1_, p.x2_);
  auto obj = [&]() {
    double r = m.expected_improvement("objective", incumbent, p);
    return r;
  };
  opt.set_objective_function(obj);
  opt.set_max_num_iterations(10000);
  opt.set_maximizing();
  opt.run_optimization();
}

void bo_naive_optim() {
  FullModel m;
  m.set_num_particles(100);
  BHParams p;
  BayesOpt<unordered_map<string, double> > opt;
  auto subopt = [&]() {
    maximize_ei(m, p, opt.best_objective());
  };
  auto util = [&](){
    unordered_map<string, double> res;
    res["objective"] = branin_hoo(p);
    PR(p.x1_.value(), p.x2_.value(), res["objective"]);
    return res;
  };
  auto learn = [&](const unordered_map<string, double>& r){
    m.observe(r, p);
  };
  opt.set_subopt_function(subopt);
  opt.set_objective_function(util);
  opt.set_learning_function(learn);
  opt.set_minimizing();
  opt.set_max_num_iterations(25);
  opt.run_optimization();
}

/// Structured way of optimizing it, more accurate model
struct FirstTermModel : public SemiParametricModel<FirstTermModel> {
  FirstTermModel() {
    alpha_ = uniform_real_distribution<>(0.0, 20.0)(generator);

    p_.default_noise(0.0);
    p_.mean(uniform_real_distribution<>(0.0, 10.0)(generator));
    p_.stdev(uniform_real_distribution<>(0.0, 200.0)(generator));
    p_.linear_scales({uniform_real_distribution<>(0.0, 15.0)(generator)});
    set_params(p_);
  }

  double parametric(double x1) const {
    return alpha_ * cos(x1);
  }

  double alpha_;
  GPParams p_;
};

struct SecondTermModel : public SemiParametricModel<SecondTermModel> {
  SecondTermModel() {
    beta_ = uniform_real_distribution<>(0.0, 20.0)(generator);

    p_.default_noise(0.0);
    p_.mean(uniform_real_distribution<>(0.0, 10.0)(generator));
    p_.stdev(uniform_real_distribution<>(0.0, 200.0)(generator));
    p_.linear_scales({uniform_real_distribution<>(0.0, 15.0)(generator),
                      uniform_real_distribution<>(0.0, 15.0)(generator)});
    set_params(p_);
  }
  double parametric(double first_term, double x1, double x2) const {
    return first_term + pow(x2, 2.0) + beta_ * pow(x1, 4.0);
  }
  vector<double> to_vec(double first_term, double x1, double x2) const {
    return {x1, x2};
  }
  double beta_;
  GPParams p_;
};

struct BHModel : public DAGModel<BHModel> {
  BHModel() {
    first_.set_num_particles(100000);
    second_.set_num_particles(100000);
  }
  void model(const BHParams& p) {
    double ft = output("first", first_, p.x1_.value());
    output("objective", second_, ft, p.x1_.value(), p.x2_.value());
  }
  ProbEngine<FirstTermModel> first_;
  ProbEngine<SecondTermModel> second_;
};

void maximize_ei(BHModel& m, BHParams& p, double incumbent) {
  NLOpt<> opt(p.x1_, p.x2_);
  auto obj = [&]() {
    double r = m.expected_improvement("objective", incumbent, p);
    return r;
  };
  opt.set_objective_function(obj);
  opt.set_max_num_iterations(10000);
  opt.set_maximizing();
  opt.run_optimization();
}

void bo_struct_optim() {
  BHModel m;
  m.set_num_particles(100);
  BHParams p;
  BayesOpt<unordered_map<string, double> > opt;
  auto subopt = [&]() {
    maximize_ei(m, p, opt.best_objective());
  };
  auto util = [&](){
    unordered_map<string, double> res;
    res["first"] = branin_hoo_first_term(p);
    res["objective"] = branin_hoo(p);
    PR(p.x1_.value(), p.x2_.value(),  res["first"], res["objective"]);
    return res;
  };
  auto learn = [&](const unordered_map<string, double>& r){
    m.observe(r, p);
  };
  opt.set_subopt_function(subopt);
  opt.set_objective_function(util);
  opt.set_learning_function(learn);
  opt.set_minimizing();
  opt.set_max_num_iterations(25);
  opt.run_optimization();
}

/// Better optimization procedure: sample from the model and optimize its
/// prediction
void maximize_sample_prediction(BHModel& m, BHParams& p) {
  NLOpt<> opt(p.x1_, p.x2_);
  auto util = [&]() {
    return m.sample_predict("objective", p);
  };
  opt.set_objective_function(util);
  opt.set_minimizing();
  opt.set_max_num_iterations(10000);
  opt.run_optimization();
}

void maximize_with_samples(BHModel& m, BHParams& p, double incumbent) {
  SimpleOpt<> opt;
  auto util = [&]() {
    auto s = m.sample();
    maximize_sample_prediction(*s, p);
    return m.expected_improvement("objective", incumbent, p);
  };
  opt.set_objective_function(util);
  opt.set_maximizing();
  opt.set_max_num_iterations(20);
  opt.run_optimization();
}

void bo_struct_optim_efficient() {
  BHModel m;
  m.set_num_particles(100);
  BHParams p;
  BayesOpt<unordered_map<string, double> > opt;
  auto subopt = [&]() {
    maximize_with_samples(m, p, opt.best_objective());
  };
  auto util = [&](){
    unordered_map<string, double> res;
    res["first"] = branin_hoo_first_term(p);
    res["objective"] = branin_hoo(p);
    PR(p.x1_.value(), p.x2_.value(),  res["first"], res["objective"]);
    return res;
  };
  auto learn = [&](const unordered_map<string, double>& r){
    m.observe(r, p);
  };
  opt.set_subopt_function(subopt);
  opt.set_objective_function(util);
  opt.set_learning_function(learn);
  opt.set_minimizing();
  opt.set_max_num_iterations(25);
  opt.run_optimization();
}

int main() {
  //bo_naive_optim();
  //bo_struct_optim();
  bo_struct_optim_efficient();
}
