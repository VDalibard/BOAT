#ifndef GP_T_HPP_INCLUDED
#define GP_T_HPP_INCLUDED

namespace boat {
double normal_lnp(double x, double mean, double variance);
extern std::mt19937 generator;

// GP
template <class Input>
double GP<Input>::observe(Input x_new, double y_new) {
  return observe(std::move(x_new), y_new, default_noise_);
}

/*
template <class Input>
void GP<Input>::set_covariance_function(
    std::function<Eigen::VectorXd(const std::vector<Input const*>& x1,
                                  const Input& x2)> covariance_function) {
  assert(observed_x_.empty());
  covariance_function_ = std::move(covariance_function);
}

template <class Input>
void GP<Input>::set_covariance_function(
    std::function<double(const Input&, const Input& x2)> covariance_function) {
  set_covariance_function([cf = std::move(covariance_function)](
      const std::vector<Input const*>& x1, const Input& x2) {
    thread_local Eigen::VectorXd res;
    res.resize(SIZE(x1));
    for (int i = 0; i < SIZE(x1); i++) {
      res(i) = cf(*x1[i], x2);
    }
    return res;
  });
}
*/
template <class Input>
double GP<Input>::observe(Input x_new, double y_new, double noise) {
  assert(!global_drawing_ && !local_drawing_);
  auto distrib = predict_distrib(x_new);
  if(distrib.second + noise * noise + 1e-6 * y_new < 0.0){
    PR(observed_x_.size(), y_new, distrib.first, distrib.second, noise);
    PRS(GOT negative variance - ASSIGNING P=0 TO THIS PARTICLE);
    return -900.0;
  }

  double lnp = normal_lnp(y_new, distrib.first,
                          distrib.second + noise * noise + 1e-6 * y_new);
  if(std::isnan(lnp)){
    PR(observed_x_.size(), y_new, distrib.first, distrib.second, noise);
    PRS(GOT NAN LIKELIHOOD - ASSIGNING P=0 TO THIS PARTICLE);
    return -900.0;
  }
  observed_x_.emplace_back(std::move(x_new));
  observe_i(y_new, noise);
  return lnp;
}

template <class Input>
std::pair<double, double> GP<Input>::predict_distrib(const Input& x_new) const {
  assert(!global_drawing_ && !local_drawing_);
  if (observed_x_.empty()) {
    return std::make_pair(mean_, amplitude_);
  }
  compute_cholesky();
  Eigen::VectorXd cov_new_old = covariance(get_vec_ptrs(observed_x_), x_new);
  double cov_new_new = covariance({&x_new}, x_new)(0, 0);
  return predict_distrib_i(cov_new_old, cov_new_new);
}

template <class Input>
double GP<Input>::predict(const Input& x_new) const {
  if(global_drawing_){
    assert(!local_drawing_);
    return predict_with_draw(x_new, global_draw_);
  } else if (local_drawing_) {
    return predict_with_draw(x_new, local_draw_);
  }
  if (observed_x_.empty()) {
    return mean_;
  }
  compute_cholesky();
  Eigen::VectorXd cov_new_old = covariance(get_vec_ptrs(observed_x_), x_new);
  return predict_i(cov_new_old);
}

template <class Input>
double GP<Input>::predict_with_draw(const Input& x_new, double draw) const {
  if (observed_x_.empty()) {
    return mean_ + draw * sqrt(amplitude_);
  }
  compute_cholesky();
    //Slightly duplicate code wih predict_distrib
  Eigen::VectorXd cov_new_old = covariance(get_vec_ptrs(observed_x_), x_new);
  double cov_new_new = covariance({&x_new}, x_new)(0, 0);
  auto p = predict_distrib_i(cov_new_old, cov_new_new);
  return p.first + draw * std::sqrt(p.second);
}

template <class Input>
double GP<Input>::sample(Input x_new) {
  return sample(std::move(x_new), default_noise_);
}

template <class Input>
double GP<Input>::sample(Input x_new, double noise) {
  auto p = predict_distrib(x_new);
  p.second += noise * noise;
  double s = std::normal_distribution<>(p.first, sqrt(p.second))(generator);
  //TODO: two versions of sample, one that observes and the other not
  //observe(std::move(x_new), s, noise);
  return s;
}

template <class Input>
void GP<Input>::compute_cholesky() const {
  if (covariance_cholesky_.rows() == SIZE(observed_x_)) {
    return;
  }
  Eigen::VectorXd new_cov =
      covariance(get_vec_ptrs(observed_x_), observed_x_.back());
  compute_cholesky_i(new_cov);
}

template <class Input>
template <class Other>
void GP<Input>::copy_observes(const std::vector<int>& indices, Other& other) {
  for (const auto i : indices) {
    other.observe(observed_x_[i], observed_y_(i) + mean_, noise_(i));
  }
}

template <class Input>
Eigen::VectorXd GP<Input>::covariance(const std::vector<Input const*>& x1,
                                      const Input& x2) const {
  Eigen::VectorXd cov = covariance_function(x1, x2);
  cov *= amplitude_;
  return cov;
}

template <class Input>
const std::vector<Input>& GP<Input>::get_inputs() const{
  return observed_x_;
}

// Tree
template <class Model>
Tree<Model>::Tree()
    : Tree(64) {}

template <class Model>
Tree<Model>::Tree(int observations_thresh)
    : observations_thresh_(observations_thresh) {}

template <class Model>
Tree<Model>::Tree(const Tree& other){
  *this = other;
}

template <class Model>
Tree<Model>& Tree<Model>::operator=(const Tree<Model>& other) {
  //Does a deep copy
  //other.check_state(); // Removed so that we can use unset trees
  observations_thresh_ = other.observations_thresh_;
  split_function_ = other.split_function_;
  side_function_ = other.side_function_;
  if(side_function_){
    //It is split
    left_ = other.left_->clone_full();
    right_ = other.right_->clone_full();
    model_ = nullptr;
  } else {
    model_ = std::make_unique<Model>(*other.model_);
    left_ = nullptr;
    right_ = nullptr;
  }
  return *this;
}

template <class Model>
std::unique_ptr<Tree<Model> > Tree<Model>::clone_full() {
  return std::make_unique<Tree<Model> > (*this);
}

template <class Model>
void Tree<Model>::set_split_function(SplitFunction split_function) {
  split_function_ = std::make_shared<SplitFunction>(std::move(split_function));
}

template <class Model>
void Tree<Model>::set_split_function(std::shared_ptr<SplitFunction> split_function) {
  split_function_ = std::move(split_function);
}

template <class Model>
void Tree<Model>::set_model(std::unique_ptr<Model> model) {
  assert(!model_);
  assert(model->get_inputs().empty());
  model_ = std::move(model);
}

template <class Model>
template <class... Args>
double Tree<Model>::observe(Input input, Args... args) {
  check_state();
  if (model_) {
    if (split_condition(input)) {
      split(input);
      // We recurs on the same object but we will go down the tree this time
      return observe(std::move(input), std::forward<Args>(args)...);
    } else {
      return model_->observe(std::move(input), std::forward<Args>(args)...);
    }
  } else {
    assert(left_ && right_);
    bool side = side_function_(input);
    if (side) {
      return right_->observe(std::move(input), std::forward<Args>(args)...);
    } else {
      return left_->observe(std::move(input), std::forward<Args>(args)...);
    }
  }
}

template <class Model>
template <class... Args>
double Tree<Model>::sample(Input input, Args... args) {
  check_state();
  if (model_) {
    if (split_condition(input)) {
      split(input);
      // We recurs on the same object but we will go down the tree this time
      return sample(std::move(input), std::forward<Args>(args)...);
    } else {
      return model_->sample(std::move(input), std::forward<Args>(args)...);
    }
  } else {
    assert(left_ && right_);
    bool side = side_function_(input);
    if (side) {
      return right_->sample(std::move(input), std::forward<Args>(args)...);
    } else {
      return left_->sample(std::move(input), std::forward<Args>(args)...);
    }
  }
}

template <class Model>
void Tree<Model>::split(const Input& new_input) {
  auto split = (*split_function_)(model_->get_inputs(), new_input);
  assert(!side_function_);
  side_function_ = std::move(std::get<2>(split));
  left_ = clone_with_params_only();
  model_->copy_observes(std::get<0>(split), *(left_->model_));

  right_ = clone_with_params_only();
  model_->copy_observes(std::get<1>(split), *(right_->model_));

  model_ = nullptr;
}

template <class Model>
template <class Func>
void Tree<Model>::execute(const Input& input, Func func) const {
  check_state();
  if (model_) {
    func(*model_);
  } else {
    bool side = side_function_(input);
    if (side) {
      return right_->execute(input, std::move(func));
    } else {
      return left_->execute(input, std::move(func));
    }
  }
}

template <class Model>
template <class Func>
void Tree<Model>::execute_all(Func func) {
  check_state();
  if(model_){
    func(*model_);
  } else {
    left_->execute_all(func);
    right_->execute_all(func);
  }
}

template <class Model>
double Tree<Model>::predict(const Input& input) const{
  double result;
  execute(input, [&](const Model& m) { result = m.predict(input); });
  return result;
}

template <class Model>
void Tree<Model>::check_state() const {
  assert(split_function_);
  if (side_function_) {
    assert(left_);
    assert(right_);
    assert(!model_);
  } else {
    assert(!left_);
    assert(!right_);
    assert(model_);
  }
}

template <class Model>
void Tree<Model>::set_observation_thresh(int observations_thresh){
  assert(!side_function_ && (!model_ || model_->get_inputs().empty()));
  observations_thresh_ = observations_thresh;
}

template <class Model>
Model& Tree<Model>::model(){
  if(!model_){
    //Construct it ourselves
    assert(!side_function_);  //Check we are just initializing
    model_ = std::make_unique<Model>();
  }
  return *model_;
}

template <class Model>
bool Tree<Model>::split_condition(const Input& input){
  return SIZE(model_->get_inputs()) == observations_thresh_;
}

template <class Model>
std::unique_ptr<Tree<Model> > Tree<Model>::clone_with_params_only() {
  auto ptr = std::make_unique<Tree<Model> >();
  ptr->copy_tree_parameters(*this);
  return std::move(ptr);
}

template <class Model>
void Tree<Model>::copy_tree_parameters(const Tree<Model>& other){
  set_observation_thresh(other.observations_thresh_);
  set_split_function(other.split_function_);
  set_model(std::make_unique<Model>());
  model().copy_model_parameters(*other.model_);
}
}
#endif  // GP_T_HPP_INCLUDED
