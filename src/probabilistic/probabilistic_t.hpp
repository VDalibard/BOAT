#ifndef PROBABILISTIC_T_HPP_INCLUDED
#define PROBABILISTIC_T_HPP_INCLUDED
#include "../utilities/useful_macros.hpp"
#include <iostream>
#include <fstream>
#include <algorithm>
#include <functional>

namespace boat {
extern std::mt19937 generator;

double log_sum_exp(const std::vector<double>& log_values);
double get_normalized_likelihoods(std::vector<double>& normalized_likelihoods,
                                  std::vector<double>& log_likelihoods);
void draw_bins_rd(std::vector<size_t>& bins,
               const std::vector<double>& normalized_likelihoods);

std::vector<int> draw_bins(
    const std::vector<double>& normalized_likelihoods,
    int num_draws);

double calc_ess(const std::vector<double>& normalized_likelihoods);


double average_normalized(const std::vector<double>& values,
                          const std::vector<double>& nl);

std::vector<double> average_normalized(
                          const std::vector<std::vector<double>>& values,
                          const std::vector<double>& nl);

double exp_imp(double mean, double var, double best, bool minimizing);

struct Averager {
  static Averager* current_averager;
  Averager(int num_particles)
      : num_particles_(num_particles),
        index_(0) {
    assert(current_averager == nullptr);
    current_averager = this;
  }

  ~Averager() {
    assert(current_averager == this);
    current_averager = nullptr;
  }

  int num_particles_;
  int index_;
};


template <class Particle>
ProbEngine<Particle>::ProbEngine(int num_particles)
    : initialized_(false),
      num_particles_(num_particles),
      num_particles_inv_(1.0 / (double)num_particles_),
      mll_(0.0) {}

template <class Particle>
ProbEngine<Particle>::ProbEngine(std::vector<Particle> particles)
    : ProbEngine(particles.size()) {
  initialize();
  particles_ = std::move(particles);
}

template <class Particle>
void ProbEngine<Particle>::initialize() const {
  assert(!initialized_);
  assert(particles_.empty());
  initialized_ = true;
  normalized_likelihoods_.reserve(num_particles_);
  log_likelihoods_.reserve(num_particles_);
  particles_.reserve(num_particles_);
  double log_likelihood = log(num_particles_inv_);
  for (size_t i = 0; i < num_particles_; i++) {
    normalized_likelihoods_.push_back(num_particles_inv_);
    log_likelihoods_.push_back(log_likelihood);
  }
}

template <class Particle>
ProbEngine<Particle>& ProbEngine<Particle>::set_num_particles(
    int num_particles) {
  if (!initialized_) {
    num_particles_ = num_particles;
    num_particles_inv_ = 1.0 / (double)num_particles_;
  } else {
    ProbEngine<Particle> tmp(num_particles);
    tmp.create_from_engines(*this);
    *this = std::move(tmp);
  }
  return *this;
}

template <class Particle>
int ProbEngine<Particle>::get_num_particles() const {
  return num_particles_;
}

template <class Particle>
ProbEngine<Particle> ProbEngine<Particle>::single_particle_engine() const {
  ProbEngine<Particle> res(1);
  const Particle& s = sample();
  res.create_using_constructor(s);
  return res;
}

template <class Particle>
template <class... Args>
void ProbEngine<Particle>::create_using_constructor(const Args&... args) {
  initialize();
  for (size_t i = 0; i < num_particles_; i++) {
    particles_.emplace_back(args...);
  }
}

template <class Particle>
template <class ForwardIt>
void ProbEngine<Particle>::create_from_range(ForwardIt begin, ForwardIt end) {
  initialize();
  ForwardIt it = begin;
  for (size_t i = 0; i < num_particles_; i++) {
    particles_.emplace_back(*it);
    it++;
    if (it == end) {
      it = begin;
    }
  }
}

template <class Particle>
template <class... Args>
void ProbEngine<Particle>::create_from_engines(
    const ProbEngine<Args>&... args) {
  initialize();
  auto inputs = std::make_tuple(args.get_random_selection(num_particles_)...);
  auto s = gen_seq<sizeof...(Args)>();
  cfe_helper(inputs, s);
}

template <class Particle>
void ProbEngine<Particle>::check_ess() {
  if (calc_ess(normalized_likelihoods_) < num_particles_ / 2.0) {
    resample();
  }
}


template <class Particle>
void ProbEngine<Particle>::check_constructed() const {
  if (particles_.empty()) {
    construction_helper();
  }
}

template <class Particle>
std::vector<const Particle*>
ProbEngine<Particle>::get_random_selection(int new_np) const {
  check_constructed();
  std::vector<const Particle*> res;
  res.reserve(new_np);
  std::vector<int> bins = draw_bins(normalized_likelihoods_, new_np);
  for(int i=0; i<SIZE(bins); i++) {
    res.insert(res.end(), bins[i], &particles_[i]);
  }
  std::shuffle(res.begin(), res.end(), generator);
  return res;
}

template <class Particle>
template <class... Args>
void ProbEngine<Particle>::observe(const Args&... args) {
  observe_i(false, args...);
}

template <class Particle>
template <class... Args>
void ProbEngine<Particle>::observe_deterministic(const Args&... args) {
  observe_i(true, args...);
}

template <class Particle>
template <class... Args>
void ProbEngine<Particle>::observe_i(bool deterministic, const Args&... args) {
  check_constructed();
  particles_for_average_.clear();
  if (!deterministic) {
    check_ess();
  }
  for (size_t p = 0; p < num_particles_; p++) {
    if (normalized_likelihoods_[p] == 0.0) {
      continue;
    }
    log_likelihoods_[p] += particles_[p].observe(args...);
    assert(!std::isnan(log_likelihoods_[p]));
    log_likelihoods_[p] = std::max(log_likelihoods_[p], -900.0);
  }

  mll_ += get_normalized_likelihoods(normalized_likelihoods_, log_likelihoods_);
}

template <class Particle>
void ProbEngine<Particle>::resample() {
  check_constructed();
  PRS(Resampling-------------------------------------------);
  std::vector<int> bins = draw_bins(normalized_likelihoods_,
                                    num_particles_);
  size_t to_copy = 0;
  size_t to_replace = 0;
  auto find_next = [&]() {
    while (to_replace < num_particles_ && bins[to_replace] > 0) {
      to_replace++;
    }
    while (to_copy < num_particles_ && bins[to_copy] < 2) {
      to_copy++;
    }
  };
  find_next();
  while (to_replace < num_particles_) {
    particles_[to_replace] = particles_[to_copy];
    bins[to_replace]++;
    bins[to_copy]--;
    find_next();
  }
  assert(to_copy == num_particles_);
  for (auto& l : normalized_likelihoods_) {
    l = num_particles_inv_;
  }
  double ll = log(num_particles_inv_);
  for (auto& log_likelihood : log_likelihoods_) {
    log_likelihood = ll;
  }
}

template <class Particle>
const Particle& ProbEngine<Particle>::sample() const {
  check_constructed();
  double rd = std::uniform_real_distribution<>(0.0, 1.0)(generator);
  double sum = 0.0;
  size_t i;
  //-1: we don't test the last one as that must be it.
  // Also prevents floating point errors is the sum doesn't add up to 1.0
  for (i = 0; i < particles_.size() - 1; i++) {
    sum += normalized_likelihoods_[i];
    if (rd < sum) {
      break;
    }
  }
  return particles_[i];
}

template <class Particle>
const Particle& ProbEngine<Particle>::value() const {
  check_constructed();
  assert(Averager::current_averager != nullptr);
  init_particles_for_average(Averager::current_averager->num_particles_);
  return *particles_for_average_[Averager::current_averager->index_];
}

template <class Particle>
void ProbEngine<Particle>::init_particles_for_average(int pfa) const {
  if(SIZE(particles_for_average_) != pfa) {
    particles_for_average_ = get_random_selection(pfa);
  }
}

template <class Particle>
template <class Func>
typename std::result_of<Func(Particle)>::type ProbEngine<Particle>::average(Func func) const {
  typedef typename std::result_of<Func(Particle)>::type RetType;
  check_constructed();
  std::vector<RetType> res;
  res.reserve(particles_.size());
  for (size_t i = 0; i < num_particles_; i++) {
    if (normalized_likelihoods_[i] == 0.0) {
      res.emplace_back();
      continue;
    }
    res.emplace_back(func(particles_[i]));
  }
  return average_normalized(res, normalized_likelihoods_);
}

template <class Particle>
template <class Func>
void ProbEngine<Particle>::execute(Func func) {
  check_constructed();
  for (size_t i = 0; i < num_particles_; i++) {
    func(particles_[i]);
  }
}

template <class Particle>
template <class... Args>
double ProbEngine<Particle>::predict(const Args&... args) const {
  if(Averager::current_averager == nullptr){
    return average([&](const Particle& p) { return p.predict(args...); });
  } else {
    return value().predict(args...);
  }
}

template <class Particle>
template <class T>
void ProbEngine<Particle>::construction_helper(typename std::enable_if<
    std::is_default_constructible<T>::value>::type* dummy) const {
  initialize();
  particles_.resize(num_particles_);
}

template <class Particle>
template <class T>
void ProbEngine<Particle>::construction_helper(typename std::enable_if<
    !std::is_default_constructible<T>::value>::type* dummy) const {
  std::cout << "Error: no constructor provided and no default constructor"
            << std::endl;
  assert(false);
}

template <class Particle>
template <class T, int... S>
void ProbEngine<Particle>::cfe_helper(const T& inputs, seq<S...>) {
  for (size_t i = 0; i < num_particles_; i++) {
    particles_.emplace_back((*std::get<S>(inputs)[i])...);
  }
}

template <class Particle>
double ProbEngine<Particle>::get_marginal_log_likelihood() const {
  return mll_;
}

template <class Particle>
template <class Archive>
void ProbEngine<Particle>::serialize(Archive& ar, const unsigned int version) {
  ar& initialized_;
  ar& num_particles_;
  ar& num_particles_inv_;
  ar& log_likelihoods_;
  ar& normalized_likelihoods_;
  ar& particles_;
  ar& mll_;
}

template <class Particle>
const Particle& ProbEngine<Particle>::first_particle() const {
  check_constructed();
  return particles_[0];
}

template <class Func>
auto average(Func func, int num_particles) {
  typedef typename std::result_of<Func()>::type RetType;
  Averager av(num_particles);
  std::vector<RetType> results;
  results.reserve(num_particles);
  for(av.index_ = 0; av.index_<num_particles; av.index_++){
    results.emplace_back(func());
  }
  RetType res = average_normalized(
      results, std::vector<double>(num_particles, 1.0 / double(num_particles)));
  return res;
}

//SemiParametricParticle

template <class Child>
const Child* SemiParametricParticle<Child>::crtp_this() const {
  return static_cast<const Child*>(this);
}

template <class Child>
Child* SemiParametricParticle<Child>::crtp_this() {
  return static_cast<Child*>(this);
}

template <class Child>
void SemiParametricParticle<Child>::set_params(GPParams params) {
  default_noise_ = params.default_noise();
  gp_.set_params(std::move(params));
}

template <class Child>
template <class... Args>
double SemiParametricParticle<Child>::predict_mean(const Args&... args) const {
  double p = crtp_this()->parametric(args...);
  std::vector<double> v = crtp_this()->to_vec(args...);
  return p + gp_.predict_mean(v);
}


template <class Child>
template <class... Args>
double SemiParametricParticle<Child>::parametric(const Args&... args) const {
 return 0.0;
}

template <class Child>
template <class... Args>
double SemiParametricParticle<Child>::parametric_noise(const Args&... args) const {
  return default_noise_;
}

template <class Child>
template <class... Args>
std::vector<double> SemiParametricParticle<Child>::to_vec(const Args&... args) const {
  return {double(args)...};
}


template <class Child>
template <class... Args>
GaussianDistrib SemiParametricParticle<Child>::predict_distrib(
    const Args&... args) const {
  double p = crtp_this()->parametric(args...);
  std::vector<double> v = crtp_this()->to_vec(args...);
  GaussianDistrib d = gp_.predict_distrib(v);
  d.mu_ += p;
  return d;
}

template <class Child>
template <class... Args>
double SemiParametricParticle<Child>::observe(double result,
                                             const Args&... args) const {
  double p = crtp_this()->parametric(args...);
  double n = crtp_this()->parametric_noise(args...);
  std::vector<double> v = crtp_this()->to_vec(args...);
  return gp_.observe(v, result - p, n);
}

// DAGModel

template <class Child>
DAGModel<Child>::DAGModel()
    : is_sample_(false),
      num_particles_(10000),
      observe_(false),
      predict_(false),
      target_(nullptr),
      past_target_(false),
      ei_(false),
      minimizing_(true),
      sample_predict_(false),
      test_pass_(false)
      {}

template <class Child>
void DAGModel<Child>::set_num_particles(int np) {
  assert(!is_sample_);
  num_particles_ = np;
}

template <class Child>
template <class... Args>
void DAGModel<Child>::observe(
    const std::unordered_map<std::string, double>& measurements,
    const Args&... args) {
  assert(!is_sample_);
  check_unset();
  observe_ = true;
  measurements_ = &measurements;
  crtp_this()->model(args...);
  assert(measurements.size() == observed_measurements_.size());
  observed_measurements_.clear();
  observe_ = false;
  check_unset();
}

template <class Child>
template <class... Args>
std::vector<Distrib> DAGModel<Child>::predict(const Args&... args){
  assert(!is_sample_);
  check_unset();
  predict_ = true;
  target_ = &dag_model_fake_name;
  run_test_pass(args...);
  reset_generator();
  auto func = [&]() {
    output_list_index_ = 0;
    crtp_this()->model(args...);
    assert(output_list_index_ == SIZE(output_list_));
    std::vector<double> v;
    std::swap(v, predict_mus_and_vars_);
    assert(v.size() == 3* output_list_.size());
    return v;
  };
  std::vector<double> avgs = average(func, num_particles_);
  assert(avgs.size() == 3* output_list_.size());
  std::vector<Distrib> res;
  for(int i = 0; i<SIZE(output_list_); i++){
    double mu = avgs[3*i];
    double var = avgs[3*i + 1] + avgs[3*i + 2] - mu * mu;
    res.emplace_back(Distrib{output_list_[i].first, mu, var});
  }
  clear_test_pass();
  predict_ = false;
  target_ = nullptr;
  check_unset();
  return res;
}

template <class Child>
template <class... Args>
double DAGModel<Child>::expected_improvement(const std::string& target,
                                             double incumbent,
                                             const Args&... args) {
  assert(!is_sample_);
  check_unset();
  ei_ = true;
  target_ = &target;
  incumbent_ = incumbent;
  run_test_pass(args...);
  reset_generator();
  auto func = [&]() {
    ei_res_ = NAN;
    output_list_index_ = 0;
    crtp_this()->model(args...);
    assert(output_list_index_ == SIZE(output_list_));
    assert(past_target_);
    past_target_ = false;
    assert(!std::isnan(ei_res_));
    return ei_res_;
  };
  double av_ei = average(func, num_particles_);
  clear_test_pass();
  ei_ = false;
  target_ = nullptr;
  check_unset();
  return av_ei;
}

template <class Child>
std::unique_ptr<Child> DAGModel<Child>::sample() {
  assert(!is_sample_);
  check_unset();
  std::unique_ptr<Child> p = std::make_unique<Child>(*crtp_this());
  p->is_sample_ = true;
  return p;
}

template <class Child>
template <class... Args>
double DAGModel<Child>::sample_predict(const std::string& target,
                                       const Args&... args) {
  assert(is_sample_);
  check_unset();
  sample_predict_ = true;
  target_ = &target;
  generator = std::mt19937(time(NULL));
  sample_predict_value_ = NAN;
  crtp_this()->model(args...);
  sample_predict_ = false;
  target_ = nullptr;
  past_target_ = false;
  check_unset();
  return sample_predict_value_;
}

template <class Child>
template <class... Args>
void DAGModel<Child>::run_test_pass(const Args&... args) {
  assert(!test_pass_ && !is_sample_ &&
         need_observe_.empty() &&
         seen_outputs_.empty() &&
         duplicated_engines_.empty() &&
         output_list_.empty());
  test_pass_ = true;
  crtp_this()->model(args...);
  test_pass_ = false;
  past_target_ = false;
  for(const auto& e: seen_outputs_){
    // All but the last output calls will need to be observed
    for(int i=0; i<SIZE(e.second) - 1; i++) {
      need_observe_.insert(e.second[i]);
    }
  }
}

template <class Child>
void DAGModel<Child>::clear_test_pass() {
  need_observe_.clear();
  seen_outputs_.clear();
  duplicated_engines_.clear();
  output_list_.clear();
  output_list_index_ = 0;
}

template <class Child>
template <class Particle, class... Args>
double DAGModel<Child>::output(const std::string& name,
                               ProbEngine<Particle>& eng,
                               const Args&... args){
  if(test_pass_) {
    return output_test_pass(name, eng, args...);
  } else if(past_target_){
    return output_past_target(name, eng, args...);
  } else if (observe_) {
    return output_observe(name, eng, args...);
  } else if(ei_) {
    return output_ei(name, eng, args...);
  } else if(predict_) {
    return output_predict(name, eng, args...);
  } else if (is_sample_) {
    return output_sample_predict(name, eng, args...);
  }
  assert(false);
}

template <class Child>
double DAGModel<Child>::output(const std::string& name,
                               double val){
  if(test_pass_) {
    return output_test_pass(name, val);
  } else if(past_target_){
    return output_past_target(name, val);
  } else if (observe_) {
    return output_observe(name, val);
  } else if(ei_) {
    return output_ei(name, val);
  } else if(predict_) {
    return output_predict(name, val);
  } else if (is_sample_) {
    return output_sample_predict(name, val);
  }
  assert(false);
}

template <class Child>
template <class Particle, class... Args>
double DAGModel<Child>::output_test_pass(const std::string& name,
                                         ProbEngine<Particle>& eng,
                                         const Args&... args){
  void* eng_addr = static_cast<void*>(&eng);
  output_list_.push_back(make_pair(name, eng_addr));
  if(!past_target_) {
    auto it = seen_outputs_.find(eng_addr);
    if(it == seen_outputs_.end()) {
      seen_outputs_.emplace(eng_addr, std::vector<std::string>{name});
    } else {
      // If it's the first time we see a repeat, we make a duplicate of the
      // engine
      if(it->second.size() == 1) {
        auto ptr = std::make_shared<ProbEngine<Particle>>(eng);
        ptr->set_num_particles(num_particles_);
        ptr->resample();
        duplicated_engines_.emplace(eng_addr, std::move(ptr));
      }
      it->second.push_back(name);
    }
    if(name == *target_) {
      past_target_ = true;
    }
  }
  return eng.first_particle().predict_mean(args...);
}

template <class Child>
double DAGModel<Child>::output_test_pass(const std::string& name,
                                         double val){
  output_list_.push_back(make_pair(name, nullptr));
  if(name == *target_) {
    past_target_ = true;
  }
  return val;
}

template <class Child>
template <class Particle, class... Args>
double DAGModel<Child>::output_predict(const std::string& name,
                                       ProbEngine<Particle>& eng,
                                       const Args&... args){
  void* eng_addr = static_cast<void*>(&eng);
  check_output_list(name, eng_addr);

  const Particle* p;
  auto it = duplicated_engines_.find(eng_addr);
  if(it != duplicated_engines_.end()) {
    ProbEngine<Particle>* de = static_cast<ProbEngine<Particle> *>(it->second.get());
    p = &de->value();
  } else {
    assert(need_observe_.count(name) == 0);
    p = &eng.value();
  }

  GaussianDistrib d = p->predict_distrib(args...);
  double dr = std::normal_distribution<>(d.mu_, sqrt(d.var_))(generator);
  if(need_observe_.count(name) > 1) {
    p->observe(dr, args...);
  }
  predict_mus_and_vars_.push_back(d.mu_);
  predict_mus_and_vars_.push_back(d.mu_ * d.mu_);
  predict_mus_and_vars_.push_back(d.var_);
  return dr;
}

template <class Child>
double DAGModel<Child>::output_predict(const std::string& name,
                                       double val){
  check_output_list(name, nullptr);
  predict_mus_and_vars_.push_back(val);
  predict_mus_and_vars_.push_back(val * val);
  predict_mus_and_vars_.push_back(0.0);
  return val;
}

template <class Child>
template <class Particle, class... Args>
double DAGModel<Child>::output_ei(const std::string& name,
                                  ProbEngine<Particle>& eng,
                                  const Args&... args){
  void* eng_addr = static_cast<void*>(&eng);
  check_output_list(name, eng_addr);

  const Particle* p;
  auto it = duplicated_engines_.find(eng_addr);
  if(it != duplicated_engines_.end()) {
    ProbEngine<Particle>* de = static_cast<ProbEngine<Particle> *>(
      it->second.get());
    p = &de->value();
  } else {
    assert(need_observe_.count(name) == 0);
    p = &eng.value();
  }

  GaussianDistrib d = p->predict_distrib(args...);
  if(name == *target_) {
    // We have found the target, now measure the expected improvement
    assert(need_observe_.count(name) == 0);
    ei_res_ = exp_imp(d.mu_, d.var_, incumbent_, minimizing_);
    past_target_ = true;
    return d.mu_;
  }else {
    double dr = std::normal_distribution<>(d.mu_, sqrt(d.var_))(generator);
    if(need_observe_.count(name) > 1) {
      p->observe(dr, args...);
    }
    return dr;
  }
}

template <class Child>
double DAGModel<Child>::output_ei(const std::string& name,
                                  double val){
  check_output_list(name, nullptr);
  if(name == *target_) {
    if(minimizing_) {
      ei_res_ = std::max(incumbent_ - val, 0.0);
    } else {
      ei_res_ = std::max(val - incumbent_, 0.0);
    }
    past_target_ = true;
  }
  return val;
}

template <class Child>
template <class Particle, class... Args>
double DAGModel<Child>::output_sample_predict(const std::string& name,
                                              ProbEngine<Particle>& eng,
                                              const Args&... args){
  void* eng_addr = static_cast<void*>(&eng);
  if(eng.get_num_particles() > 1) {
    assert(fixed_draws_.count(eng_addr) == 0);
    eng = eng.single_particle_engine();
    double d = std::normal_distribution<>(0.0, 1.0)(generator);
    fixed_draws_.emplace(eng_addr, d);
  }
  auto it = fixed_draws_.find(eng_addr);
  assert(it != fixed_draws_.end());
  GaussianDistrib d = eng.first_particle().predict_distrib(args...);
  double res = d.mu_ + it->second * sqrt(d.var_);
  if(name == *target_){
    sample_predict_value_ = res;
    past_target_ = true;
  }
  return res;
}

template <class Child>
double DAGModel<Child>::output_sample_predict(const std::string& name,
                                              double val){
  if(name == *target_){
    sample_predict_value_ = val;
    past_target_ = true;
  }
  return val;
}

template <class Child>
template <class Particle, class... Args>
double DAGModel<Child>::output_past_target(const std::string& name,
                                           ProbEngine<Particle>& eng,
                                           const Args&... args){
  if(ei_) {
    void* eng_addr = static_cast<void*>(&eng);
    check_output_list(name, eng_addr);
  } else {
    assert(sample_predict_);
  }
  return eng.value().predict_mean(args...);
}

template <class Child>
double DAGModel<Child>::output_past_target(const std::string& name,
                                           double val){
  if(ei_) {
    check_output_list(name, nullptr);
  } else {
    assert(sample_predict_);
  }
  return val;
}

template <class Child>
template <class Particle, class... Args>
double DAGModel<Child>::output_observe(const std::string& name,
                                       ProbEngine<Particle>& eng,
                                       const Args&... args){
  auto it = measurements_->find(name);
  assert(it != measurements_->end());
  assert(observed_measurements_.count(name) == 0);
  observed_measurements_.insert(name);
  eng.observe_deterministic(it->second, args...);
  return it->second;
}

template <class Child>
double DAGModel<Child>::output_observe(const std::string& name,
                                       double val){
  return val;
}

template <class Child>
void DAGModel<Child>::check_output_list(const std::string& name,
                                          void* eng_addr){
  // Checks the names and engines of the output is consistent within an average
  assert(output_list_index_ >= 0 && output_list_index_< SIZE(output_list_));
  assert(name == output_list_[output_list_index_].first);
  assert(eng_addr == output_list_[output_list_index_].second);
  output_list_index_++;
}

template <class Child>
void DAGModel<Child>::reset_generator(){
  generator = std::mt19937();
}

template <class Child>
void DAGModel<Child>::check_unset() {
  assert(!observe_);
  assert(!predict_);
  assert(!past_target_);
  assert(!ei_);
  assert(!sample_predict_);
  assert(!test_pass_);
  assert(target_ == nullptr);
  assert(need_observe_.empty() &&
         seen_outputs_.empty() &&
         duplicated_engines_.empty() &&
         output_list_.empty());
}

template <class Child>
const Child* DAGModel<Child>::crtp_this() const{
  return static_cast<const Child*>(this);
}

template <class Child>
Child* DAGModel<Child>::crtp_this(){
  return static_cast<Child*>(this);
}

}
#endif  // PROBABILISTIC_T_HPP_INCLUDED
