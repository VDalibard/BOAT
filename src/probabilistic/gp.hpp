#ifndef GP_HPP_INCLUDED
#define GP_HPP_INCLUDED
#include <vector>
#include <memory>
#include <iostream>
#include <Eigen/Dense>
#include <Eigen/Cholesky>

namespace boat {
double normal_lnp(double x, double mean, double variance);
typedef enum { SQUARED_EXP, MATERN52 } kernel_t;

class GPParams{
  public:
  GPParams():mean_(0.0), amplitude_(1.0), default_noise_(0.1),
  kernel_(MATERN52) {}

  void amplitude(double amplitude) {
    amplitude_ = amplitude;
  }

  void stdev(double stdev) {
    amplitude_ = stdev * stdev;
  }

  void default_noise(double default_noise) {
    assert(default_noise >= 0.0);
    //default_noise_ = std::max(default_noise, 1e-4);
    default_noise_ = default_noise;
  }

  void mean(double mean) {
    mean_ = mean;
  }

  void linear_scales(const std::vector<double>& linear_scales){
    linear_scales_ = linear_scales;
    inv_linear_scales_ = Eigen::VectorXd::Map(linear_scales.data(), linear_scales.size())
                           .array()
                           .inverse();
  }

  void kernel(kernel_t kernel){
    kernel_ = kernel;
  }


  double stdev() const {
    return sqrt(amplitude_);
  }

  double amplitude() const {
    return amplitude_;
  }

  double mean() const {
    return mean_;
  }

  double default_noise() const {
    return default_noise_;
  }

  const Eigen::VectorXd& inv_linear_scales() const {
    return inv_linear_scales_;
  }

  const std::vector<double>& linear_scales() const {
    return linear_scales_;
  }

  kernel_t kernel() const {
    return kernel_;
  }

  private:
  double mean_;
  double amplitude_;
  double default_noise_;
  std::vector<double> linear_scales_;
  Eigen::VectorXd inv_linear_scales_;
  kernel_t kernel_;
};

struct GaussianDistrib {
  double mu_;
  double var_;
  double mu() {
    return mu_;
  }
  double var() {
    return var_;
  }
  double stdev() {
    return sqrt(var_);
  }
};

class GP {
  friend class TreedGPS;
  public:
  GP();
  GP(GPParams params);
  void set_params(GPParams params);
  int num_dims() const;
  double observe(const std::vector<double>& x_new, double y_new);
  double observe(const std::vector<double>& x_new, double y_new, double noise);

  double predict_mean(const std::vector<double>& x_new) const;
  GaussianDistrib predict_distrib(const std::vector<double>& x_new) const;
  void print() const;

   //private:
   double observe_i(const std::vector<double>& x_new,
                    double y_new, double noise_var);
   void compute_cholesky() const;
   int at_index(const Eigen::RowVectorXd& v) const;
   void remove_observation(int i);
   Eigen::MatrixXd covariance(const Eigen::MatrixXd& x1,
                              const Eigen::MatrixXd& x2) const;

  void compute_cholesky_from_scratch() const;

  // Parameters
  GPParams params_;

  // Processed data
  Eigen::MatrixXd observed_x_;
  Eigen::VectorXd observed_y_;
  Eigen::VectorXd noise_var_;

  mutable Eigen::MatrixXd covariance_cholesky_;
  mutable Eigen::VectorXd alpha_;
};

class TreedGPS {
  public:
  TreedGPS();
  TreedGPS(GPParams params, int observation_thresh = 64, int overlap = 3);
  TreedGPS(const TreedGPS& other);
  TreedGPS& operator=(const TreedGPS& other);
  void set_params(GPParams params);

  double observe(const std::vector<double>& x_new, double y_new);
  double observe(const std::vector<double>& x_new, double y_new, double noise);

  double predict_mean(const std::vector<double>& x_new) const;
  GaussianDistrib predict_distrib(const std::vector<double>& x_new) const;

  private:
  bool side(const std::vector<double>& x_new) const;

  //Lots of helpers for the splitting
  std::pair<std::vector<int>, std::vector<int>> compute_thresh_and_divide(double max_diff);
  std::pair<std::vector<int>, std::vector<int>> compute_thresh_props();
  void check_thresh();
  void copy_observation(bool s, int i);
  void compute_children_cholesky();

  int observation_thresh_;
  int overlap_;

  int thresh_dim_;
  double thresh_;
  double inv_ls_along_thresh_dim_;
  std::unique_ptr<GP> gp_;
  std::unique_ptr<TreedGPS> left_;
  std::unique_ptr<TreedGPS> right_;
};
}
#endif  // GP_HPP_INCLUDED
