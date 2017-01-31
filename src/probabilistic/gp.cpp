#include <iostream>
#include "../utilities/useful_macros.hpp"
#include "gp.hpp"
using Eigen::MatrixXd;
using Eigen::VectorXd;
using Eigen::RowVectorXd;
using Eigen::LLT;
using Eigen::Lower;
using Eigen::Upper;
using std::make_unique;
using std::pair;
using std::vector;
using std::min;
using std::max;
namespace boat {
extern std::mt19937 generator;
double normal_lnp(double x, double mean, double variance);

MatrixXd squared_distance(const MatrixXd& x1, const MatrixXd& x2) {
  // x1 is NxD, x2 is MxD
  // returns and NxM matrix of the distances between each row
  MatrixXd p = -2.0 * x1 * x2.transpose();
  p.colwise() += x1.rowwise().squaredNorm();
  p.rowwise() += x2.rowwise().squaredNorm().transpose();
  return p;
}

MatrixXd squared_exp(const MatrixXd& x1, const MatrixXd& x2) {
  MatrixXd dist = squared_distance(x1, x2);
  dist *= -0.5;
  dist = dist.array().exp();
  return dist;
}

MatrixXd matern52(const MatrixXd& x1, const MatrixXd& x2) {
  thread_local MatrixXd r2;
  thread_local MatrixXd r;
  thread_local MatrixXd sq5r;
  r2 = squared_distance(x1, x2).cwiseAbs();
  r = r2.array().sqrt();
  sq5r = sqrt(5.0) * r;
  MatrixXd res = MatrixXd::Ones(r.rows(), r.cols()); //Gets RVO
  res += sq5r;
  res += (5.0 / 3.0) * r2;
  res.array() *= (-sq5r).array().exp();
  return res;
}

int GPScalars2::num_dims() const {
  assert(!params_.linear_scales().empty());
  return params_.linear_scales().size();
}

MatrixXd GPScalars2::covariance(const MatrixXd& x1, const MatrixXd& x2) const{
  // Inputs must already be scaled
  static thread_local MatrixXd cov;
  if (params_.kernel() == SQUARED_EXP) {
    cov = squared_exp(x1, x2);
  } else if (params_.kernel() == MATERN52) {
    cov = matern52(x1, x2);
  } else {
    assert(false);
  }
  cov *= params_.amplitude();
  return cov;
}

void GPScalars2::compute_cholesky() const {
  int new_size = observed_x_.rows();
  int old_size = covariance_cholesky_.rows();
  if (old_size == new_size) {
    return;
  }
  assert(old_size == new_size - 1);
  Eigen::VectorXd new_cov =
      covariance(observed_x_, observed_x_.row(new_size - 1));



  Eigen::VectorXd new_dec =
      covariance_cholesky_.triangularView<Eigen::Lower>().solve(
          new_cov.head(new_size - 1));

  double self_cov = new_cov(new_size - 1) +
                    noise_var_(new_size - 1) +
                    1e-6 * params_.amplitude();
  double corner = sqrt(self_cov - new_dec.squaredNorm());
  covariance_cholesky_.conservativeResize(new_size, new_size);
  covariance_cholesky_.bottomLeftCorner(1, new_size - 1) = new_dec.transpose();
  covariance_cholesky_.bottomRightCorner<1, 1>()(0, 0) = corner;

  alpha_ = covariance_cholesky_.triangularView<Lower>().solve(observed_y_);
  alpha_ = covariance_cholesky_.transpose().triangularView<Upper>()
      .solve(alpha_);
}

void GPScalars2::compute_cholesky_from_scratch() const{
  int s = observed_x_.rows();
  MatrixXd covariance_mat = covariance(observed_x_, observed_x_);
  for(int i=0; i<s; i++) {
    covariance_mat(i, i) += noise_var_(i) +  1e-6 * params_.amplitude();
  }
  LLT<MatrixXd> llt_of_cov(covariance_mat);
  covariance_cholesky_ = llt_of_cov.matrixLLT();
  alpha_ = covariance_cholesky_.triangularView<Lower>().solve(observed_y_);
  alpha_ = covariance_cholesky_.transpose().triangularView<Upper>()
      .solve(alpha_);
}

double GPScalars2::predict_mean(const std::vector<double>& x_new) const {
  assert(SIZE(x_new) == num_dims());
  if (observed_x_.rows() == 0) {
    return params_.mean();
  }
  compute_cholesky();
  static thread_local RowVectorXd eigen_x_new;
  eigen_x_new = RowVectorXd::Map(x_new.data(), num_dims());
  eigen_x_new *= params_.inv_linear_scales().asDiagonal();
  VectorXd cov_new_old = covariance(observed_x_, eigen_x_new);
  double mu = alpha_.dot(cov_new_old);
  return params_.mean() + mu;
}

GaussianDistrib GPScalars2::predict_distrib(const std::vector<double>& x_new) const {
  if (observed_x_.rows() == 0) {
    return {params_.mean(), params_.amplitude()};
  }

  compute_cholesky();
  static thread_local RowVectorXd eigen_x_new;
  eigen_x_new = RowVectorXd::Map(x_new.data(), num_dims());
  eigen_x_new *= params_.inv_linear_scales().asDiagonal();
  VectorXd cov_new_old = covariance(observed_x_, eigen_x_new);
  double cov_new_new = covariance(eigen_x_new, eigen_x_new)(0, 0);

  double mu = alpha_.dot(cov_new_old);

  VectorXd v =
      covariance_cholesky_.triangularView<Eigen::Lower>().solve(cov_new_old);
  double var = cov_new_new - v.squaredNorm();
  return {params_.mean() + mu, var};
}

int GPScalars2::at_index(const Eigen::RowVectorXd& v) const {
  for(int i = 0; i < observed_x_.rows(); i++) {
    if(observed_x_.row(i) == v){
      return i;
    }
  }
  return -1;
}

void GPScalars2::remove_observation(int i) {

  // Now we need to update the Cholesky decomposition
  // The bottom right corner is the tricky one
  int k = observed_x_.rows() - i -1;
  if(k > 0) {
    VectorXd removed_col = covariance_cholesky_.block(i+1, i, k, 1);

    // Internal function to do the rank update, otherwise not exposed in
    // Eigen's API
    MatrixXd br = covariance_cholesky_.bottomRightCorner(k, k);

    Eigen::internal::llt_inplace<MatrixXd::Scalar, Eigen::Lower>::
        rankUpdate(br, removed_col, 1);
    //Debug
    /*
    MatrixXd tmp = br * br.transpose();
    tmp += removed_col * removed_col.transpose();
    LLT<MatrixXd> lltoftmp(tmp);
    MatrixXd right_result = lltoftmp.matrixLLT();
    PR(br);
    PR(right_result);
    */
    covariance_cholesky_.bottomRightCorner(k, k) = br;
  }

  for(int j = i+1; j < observed_x_.rows(); j++) {
    observed_x_.row(j-1) = observed_x_.row(j);
    observed_y_(j-1) = observed_y_(j);
    noise_var_(j-1) = noise_var_(j);
    covariance_cholesky_.row(j-1) = covariance_cholesky_.row(j);
    covariance_cholesky_.col(j-1) = covariance_cholesky_.col(j);
  }

  int ns = observed_x_.rows() -1;
  observed_x_.conservativeResize(ns, num_dims());
  observed_y_.conservativeResize(ns);
  noise_var_.conservativeResize(ns);
  covariance_cholesky_.conservativeResize(ns, ns);
}

double GPScalars2::observe_i(const std::vector<double>& x_new, double y_new,
                             double noise_var_new){
  GaussianDistrib distrib = predict_distrib(x_new);
  double var = distrib.var_ + noise_var_new;

  // Debug
  if(var < 0.0){
    //PR(observed_x_.size(), y_new, distrib.first, distrib.var_, noise);
    PRS(GOT negative variance - ASSIGNING P=0 TO THIS PARTICLE);
    return -900.0;
  }

  double lnp = normal_lnp(y_new, distrib.mu_, var);

  // Debug
  if(std::isnan(lnp)){
    //PR(observed_x_.size(), y_new, distrib.first, distrib.var_, noise);
    PRS(GOT NAN LIKELIHOOD - ASSIGNING P=0 TO THIS PARTICLE);
    return -900.0;
  }

  // Scale the new input
  static thread_local RowVectorXd eigen_x_new;
  eigen_x_new = RowVectorXd::Map(x_new.data(), num_dims());
  eigen_x_new *= params_.inv_linear_scales().asDiagonal();

  int ai = at_index(eigen_x_new);
  if (ai != -1) {
    double old_y = observed_y_(ai);
    double old_noise_var = noise_var_(ai);
    remove_observation(ai);
    double noise_var_new2 = 1.0 / (1.0 / old_noise_var + 1.0 / noise_var_new);
    y_new = noise_var_new2 * (old_y / old_noise_var + y_new / noise_var_new);
    noise_var_new = noise_var_new2;
  }
  // Add the new input to the data
  int new_size = observed_x_.rows() + 1;
  observed_x_.conservativeResize(new_size, num_dims());
  observed_x_.row(new_size - 1) = eigen_x_new;
  observed_y_.conservativeResize(new_size);
  observed_y_(new_size - 1) = y_new - params_.mean();
  noise_var_.conservativeResize(new_size);
  noise_var_(new_size - 1) = noise_var_new;
  return lnp;
}

double GPScalars2::observe(const std::vector<double>& x_new, double y_new){
  return observe(x_new, y_new, params_.default_noise());
}

double GPScalars2::observe(const std::vector<double>& x_new, double y_new,
                           double noise){
  return observe_i(x_new, y_new, noise * noise);
}

GPScalars2::GPScalars2() : GPScalars2(GPParams()) {}

GPScalars2::GPScalars2(GPParams params) : params_(std::move(params)) {}

void GPScalars2::set_params(GPParams params) {
  assert(observed_x_.rows() == 0);
  params_ = std::move(params);
}

//TreedGPS

TreedGPS::TreedGPS():TreedGPS(GPParams()){}

TreedGPS::TreedGPS(GPParams params, int observation_thresh, int overlap)
    : observation_thresh_(observation_thresh),
      overlap_(overlap),
      gp_(make_unique<GPScalars2>(std::move(params))){
  assert(overlap_ < observation_thresh_ / 2);
}

TreedGPS::TreedGPS(const TreedGPS& other) {
  *this = other;
}

TreedGPS& TreedGPS::operator=(const TreedGPS& other){
  observation_thresh_ = other.observation_thresh_;
  overlap_ = other.overlap_;
  thresh_dim_ = other.thresh_dim_;
  thresh_ = other.thresh_;
  inv_ls_along_thresh_dim_ = other.inv_ls_along_thresh_dim_;
  if(other.gp_){
    assert(!other.right_ && !other.left_);
    gp_ = make_unique<GPScalars2>(*other.gp_);
  } else {
    assert(other.right_ && other.left_);
    left_ = make_unique<TreedGPS>(*other.left_);
    right_ = make_unique<TreedGPS>(*other.right_);
  }
  return *this;
}

void TreedGPS::set_params(GPParams params) {
  assert(gp_);
  gp_->set_params(std::move(params));
}

double TreedGPS::predict_mean(const std::vector<double>& x_new) const {
  if(gp_) {
    return gp_->predict_mean(x_new);
  } else {
    assert(right_ && left_);
    if(side(x_new)) {
      return right_->predict_mean(x_new);
    } else {
      return left_->predict_mean(x_new);
    }
  }
}

GaussianDistrib TreedGPS::predict_distrib(const std::vector<double>& x_new) const {
  if(gp_) {
    return gp_->predict_distrib(x_new);
  } else {
    assert(right_ && left_);
    if(side(x_new)) {
      return right_->predict_distrib(x_new);
    } else {
      return left_->predict_distrib(x_new);
    }
  }
}

double TreedGPS::observe(const std::vector<double>& x_new, double y_new) {
  check_thresh();
  if(gp_) {
    return gp_->observe(x_new, y_new);
  } else {
    assert(right_ && left_);
    if(side(x_new)) {
      return right_->observe(x_new, y_new);
    } else {
      return left_->observe(x_new, y_new);
    }
  }
}

double TreedGPS::observe(const std::vector<double>& x_new,
                         double y_new, double noise) {
  check_thresh();
  if(gp_) {
    return gp_->observe(x_new, y_new, noise);
  } else {
    assert(right_ && left_);
    if(side(x_new)) {
      return right_->observe(x_new, y_new, noise);
    } else {
      return left_->observe(x_new, y_new, noise);
    }
  }
}

pair<vector<int>, vector<int>> TreedGPS::compute_thresh_and_divide(double max_diff) {
  vector<pair<double, int> > v;
  assert(observation_thresh_ == gp_->observed_x_.rows());
  for(int i = 0; i<observation_thresh_; i++) {
    v.emplace_back(gp_->observed_x_(i, thresh_dim_), i);
  }

  auto comp = [](const auto& a, const auto& b){
    return a.first < b.first;
  };

  std::sort(v.begin(), v.end(), comp);
  assert(v[observation_thresh_ - 1].first - v[0].first == max_diff);

  //Splitting twice, we have a box centered on the median
  constexpr double r = 0.5 * (3.0 - sqrt(5.0));
  int middle = observation_thresh_ * r + 1;

  // If all the values were different, we would be almost done
  // But we must be more careful in case some values are identical
  auto p = std::equal_range(v.begin(), v.end(), v[middle], comp);
  int low = std::distance(v.begin(), p.first);
  int high = std::distance(v.begin(), p.second);
  assert(low <= middle && high > middle);
  int index;
  if(low == 0 || high - 1 - middle < middle - low) {
    index = high;
  }else {
    index = low;
  }
  assert(index > 0 && index < observation_thresh_ &&
         v[index].first > v[index - 1].first);
  thresh_ = (v[index].first + v[index -1].first)/2.0;

  if(thresh_ <= v[index -1].first) {
    // Can happen due to limited precision of doubles
    // thresh_ shoudl be stricly greater than v(index -1))
    thresh_ = std::nexttoward(v[index -1].first, v[index].first);
  }

  // Compute the set of elements that go on each side
  int end_left = min(observation_thresh_-2, index + overlap_ -1);
  int start_right = max(1, index - overlap_);
  pair<vector<int>, vector<int>> res;
  for(int i = 0; i<observation_thresh_; i++) {
    if(i < end_left) {
      res.first.push_back(i);
    }
    if(i > start_right) {
      res.second.push_back(i);
    }
  }
  return res;
}

pair<vector<int>, vector<int>> TreedGPS::compute_thresh_props() {
  // Find the dimension
  thresh_dim_ = 0;
  double max_diff = -1.0;
  for(int i = 0; i < gp_->num_dims(); i++) {
    double diff = gp_->observed_x_.col(i).maxCoeff() -
                    gp_->observed_x_.col(i).minCoeff();
    if(diff > max_diff){
      thresh_dim_ = i;
      max_diff = diff;
    }
  }
  assert(max_diff > 0.0);
  auto res = compute_thresh_and_divide(max_diff);
  inv_ls_along_thresh_dim_ = gp_->params_.inv_linear_scales()(thresh_dim_);
  return res;
}

void TreedGPS::copy_observation(bool s, int i) {
  GPScalars2& child = s ?  *(right_->gp_) : *(left_->gp_);
  int ns = child.observed_x_.rows() + 1;
  child.observed_x_.conservativeResize(ns, gp_->num_dims());
  child.observed_y_.conservativeResize(ns);
  child.noise_var_.conservativeResize(ns);
  child.observed_x_.row(ns - 1) = gp_->observed_x_.row(i);
  child.observed_y_(ns-1) = gp_->observed_y_(i);
  child.noise_var_(ns-1) = gp_->noise_var_(i);
}

void TreedGPS::check_thresh() {
  if(!gp_) {return;}
  if(gp_->observed_x_.rows() < observation_thresh_){ return;}
  assert(gp_->observed_x_.rows() == observation_thresh_);

  auto split = compute_thresh_props();

  // Now we need to construct the two children and give them their datapoints
  left_ = make_unique<TreedGPS>(gp_->params_, observation_thresh_, overlap_);
  right_ = make_unique<TreedGPS>(gp_->params_, observation_thresh_, overlap_);
  for(int i : split.first) {
    copy_observation(false, i);
  }
  for(int i : split.second) {
    copy_observation(true, i);
  }

  // Setup the cholesky
  left_->gp_->compute_cholesky_from_scratch();
  right_->gp_->compute_cholesky_from_scratch();

  //Finally delete the local GP
  gp_ = nullptr;
}

bool TreedGPS::side(const vector<double>& x_new) const {
  assert(thresh_dim_ >= 0 && thresh_dim_ < SIZE(x_new));
  return x_new[thresh_dim_] * inv_ls_along_thresh_dim_ >= thresh_;
}
}
