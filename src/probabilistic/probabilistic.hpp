#ifndef PROBABILISTIC_HPP_INCLUDED
#define PROBABILISTIC_HPP_INCLUDED
#include <vector>
#include <memory>
#include <random>
#include <math.h>
#include <cassert>
#include <unordered_map>
#include <unordered_set>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include "../probabilistic/gp.hpp"
#include "../utilities/utilities.hpp"

namespace boat {
double normal_lnp(double x, double mean, double variance);

template <class Particle> class SemiParamericEngine;
extern std::string dag_model_fake_name;

template <class Particle>
class ProbEngine {
  friend class boost::serialization::access;
  friend class SemiParamericEngine<Particle>;
 public:
  ProbEngine(int num_particles = 50000);

  ProbEngine(std::vector<Particle> particles);

  ProbEngine& set_num_particles(int num_particles);  // Resamples if
                                                        // initialized
  int get_num_particles() const;

  ProbEngine single_particle_engine() const ;


  template <class... Args>
  void create_using_constructor(const Args&... args);

  template <class ForwardIt>
  void create_from_range(ForwardIt begin, ForwardIt end);

  template <class... Args>
  void create_from_engines(const ProbEngine<Args>&... args);

  const Particle& value() const;

  template <class... Args>
  void observe(const Args&... args);

  template <class... Args>
  void observe_deterministic(const Args&... args);

  template <class... Args>
  void observe_i(bool deterministic, const Args&... args);

  template <class Func>
  typename std::result_of<Func(Particle)>::type average(Func func) const;

  template <class Func>
  void execute(Func func);

  template <class... Args>
  double predict(const Args&... args) const;

  double get_marginal_log_likelihood() const;

  template <class Archive>
  void serialize(Archive& ar, const unsigned int version);

  void resample();

  const Particle& sample() const;
  const Particle& first_particle() const;

private:
  void initialize() const;
  void check_ess();
  void check_constructed() const;
  std::vector<const Particle*> get_random_selection(int new_ns) const;
  void init_particles_for_average(int pfa) const;

  template <class T = Particle>
  void construction_helper(
      typename std::enable_if<std::is_default_constructible<T>::value>::type*
          dummy = nullptr) const;

  template <class T = Particle>
  void construction_helper(
      typename std::enable_if<!std::is_default_constructible<T>::value>::type*
          dummy = nullptr) const;

  template <class T, int... S>
  void cfe_helper(const T& inputs, seq<S...>);


  mutable bool initialized_;
  size_t num_particles_;
  double num_particles_inv_;
  mutable std::vector<double> log_likelihoods_;
  mutable std::vector<double> normalized_likelihoods_;
  mutable std::vector<Particle> particles_;
  mutable std::vector<const Particle*> particles_for_average_;
  mutable double mll_;
};

// Semi-parametric utilities

template <class Particle>
class SemiParamericEngine;

template <class Child>
class DAGModel;
// Todo, turn this into a class with multiple fields.
// Add an expected improvement possibility. When the name of the output is
// "objective", return the expected improvement rather than an sample.

struct Distrib {
  std::string name_;
  double mu_;
  double var_;
  friend std::ostream &operator<<(std::ostream& os, const Distrib& d){
    os << d.name_ << ":  mu " << d.mu_ << " std " << sqrt(d.var_);
    return os;
  }
};

template <class Child>
class SemiParametricModel {
  template <class Particle> friend class SemiParamericEngine;
public:
  template<class... Args>
  GaussianDistrib predict_distrib(const Args&... args) const;

  template<class... Args>
  double predict_mean(const Args&... args) const;

  template<class... Args>
  double observe(double result, const Args&... args) const;
protected:

  template <class... Args>
  double parametric(const Args&... args) const;

  template <class... Args>
  double parametric_noise(const Args&... args) const;

  template <class... Args>
  std::vector<double> to_vec(const Args&... args) const;

  void set_params(GPParams params);

//private:
  Child* crtp_this();
  const Child* crtp_this() const;

  double default_noise_;
  mutable TreedGPS gp_;
};

template <class Child>
class DAGModel {
public:
  DAGModel();

  void set_num_particles(int np);
  template <class... Args>
  void observe(const std::unordered_map<std::string, double>& measurements,
               const Args&... args);

  template <class... Args>
  std::vector<Distrib> predict(const Args&... args);

  template <class... Args>
  double expected_improvement(const std::string& target, double incumbent,
                              const Args&... args);

  std::unique_ptr<Child> sample();
  template <class... Args>
  double sample_predict(const std::string& target, const Args&... args);


protected:
  template <class Particle, class... Args>
  double output(const std::string& name, ProbEngine<Particle>& eng,
                const Args&... args);

  double output(const std::string& name, double value);
private:

  template <class... Args>
  void run_test_pass(const Args&... args);

  void clear_test_pass();
  //Output functions for each situation
  template <class Particle, class... Args>
  double output_test_pass(const std::string& name, ProbEngine<Particle>& eng,
                          const Args&... args);

  double output_test_pass(const std::string& name, double val);

  template <class Particle, class... Args>
  double output_predict(const std::string& name, ProbEngine<Particle>& eng,
                          const Args&... args);

  double output_predict(const std::string& name, double val);

  template <class Particle, class... Args>
  double output_ei(const std::string& name, ProbEngine<Particle>& eng,
                          const Args&... args);

  double output_ei(const std::string& name, double val);

  template <class Particle, class... Args>
  double output_sample_predict(const std::string& name,
                               ProbEngine<Particle>& eng,
                               const Args&... args);

  double output_sample_predict(const std::string& name, double val);

  template <class Particle, class... Args>
  double output_past_target(const std::string& name, ProbEngine<Particle>& eng,
                          const Args&... args);

  double output_past_target(const std::string& name, double val);

  template <class Particle, class... Args>
  double output_observe(const std::string& name, ProbEngine<Particle>& eng,
                          const Args&... args);

  double output_observe(const std::string& name, double val);


  void check_unset();
  void check_output_list(const std::string& name, void* eng_addr);
  // This following function resets the random number generator so
  // results are the same if we give the same input. This makes debugging easier
  void reset_generator();
  Child* crtp_this();
  const Child* crtp_this() const;


  bool is_sample_;
  int num_particles_;


  //Observe
  bool observe_;
  const std::unordered_map<std::string, double>* measurements_;
  std::unordered_set<std::string> observed_measurements_;


  // Predict
  bool predict_;
  std::vector<double> predict_mus_and_vars_;

  // Useful for the next two that have targets
  const std::string* target_;
  bool past_target_;

  //EI
  bool ei_;
  double incumbent_;
  double ei_res_;
  bool minimizing_;

  //Sample predict
  bool sample_predict_;
  std::unordered_map<void *, double> fixed_draws_;
  double sample_predict_value_;


  // Copying necessary engines
  bool test_pass_;
  std::unordered_set<std::string> need_observe_;
  std::unordered_map<void *, std::vector<std::string> > seen_outputs_;
  std::unordered_map<void*, std::shared_ptr<void> > duplicated_engines_;
  std::vector<std::pair<std::string, void*> > output_list_;
  int output_list_index_;
};


}
#include "probabilistic_t.hpp"

#endif  // PROBABILISTIC_HPP_INCLUDED
