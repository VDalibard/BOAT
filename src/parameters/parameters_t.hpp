#ifndef PARAMETERS_T_HPP_INCLUDED
#define PARAMETERS_T_HPP_INCLUDED

#include <memory>
namespace boat{
template <class T>
ParameterInstance* TopParameterPtr<T>::init_previous() {
  // Computed before we initialize the parameter
  assert(computation_type == MAIN);
  ParameterInstance* previous_instance = ParameterInstance::current_instance;
  if (previous_instance != nullptr) {
    previous_instance->set_fixed();
  }
  ParameterInstance::current_instance = root_instance_.get();

  return previous_instance;
}

template <class T>
template <class... Args>
TopParameterPtr<T>::TopParameterPtr(Args&&... args)
    : root_instance_(new ParameterInstance()),
      previous_instance_(init_previous()),
      param_(std::make_unique<T>(std::forward<Args>(args)...)) {
  // Watch out for init_previous above
}

template <class T>
void TopParameterPtr<T>::set_previous_instance_current() {
  // This is true in most parctical cases so leaving it here for debug for now
  // But should probably be removed
  assert(ParameterInstance::current_instance == root_instance_.get());

  // Should probably use the latest child,
  // but again same in most practical cases
  root_instance_->set_fixed();
  ParameterInstance::current_instance = previous_instance_;
}

template <class T>
TopParameterPtr<T>::~TopParameterPtr() {
  // This is true in most parctical cases so leaving it here for debug for now
  // But should probably be removed
  assert(! root_instance_->has_children());

  if(ParameterInstance::current_instance == root_instance_.get()){
    ParameterInstance::current_instance = nullptr;
  }
}

template <class T>
const std::unordered_map<ParameterInstance*, T> Parameter<T>::dummy_;

template <class T>
Parameter<T>::Parameter()
    : fixed_instance_(nullptr), cached_iterator_(dummy_.end()) {
  // Some debug
  if (computation_type == PARAMETER) {
    assert(parameter_function_handler->is_calling_new());
  }
}

template <class T>
Parameter<T>::~Parameter() {
  while (!instances_.empty()) {
    remove_instance_internal(instances_.begin()->first, true);
  }
}

template <class T>
void Parameter<T>::remove_instance(ParameterInstance* instance) {
  remove_instance_internal(instance, false);
}

template <class T>
void Parameter<T>::remove_instance_internal(ParameterInstance* instance,
                                            bool deleting_obj) {
  auto it = instances_.find(instance);
  assert(it != instances_.end());
  instances_.erase(it);
  cached_iterator_ = dummy_.end();
  instance->remove_parameter(this);
}

template <class T>
void Parameter<T>::promote_instance(ParameterInstance* parent,
                                    ParameterInstance* child) {
  // We are making the child instance replace the parent one

  assert(child->parent() == parent);
  auto it = instances_.find(child);
  assert(it != instances_.end());
  assert(instances_.find(parent) == instances_.end());

  // Move let us do it with unique ptrs
  instances_.emplace(parent, std::move(it->second));
  instances_.erase(it);

  child->remove_parameter(this);
  parent->add_parameter(this);

  cached_iterator_ = dummy_.end();
}

template <class T>
ParameterInstance* Parameter<T>::is_instantiated(
    ParameterInstance* inst) const {
  auto it = find_instance(inst);
  if (it == instances_.end()) {
    return nullptr;
  } else {
    return it->first;
  }
}

template <class T>
void Parameter<T>::set_fixed() {
  fixed_instance_ = ParameterInstance::current_instance;
}

template <class T>
void Parameter<T>::unset_fixed() {
  fixed_instance_ = nullptr;
}

template <class T>
void Parameter<T>::assign(T val) {
  assert(check_value(val) && "Assigning an unallowed value");

  // TODO: Debug this all works, Paramfuncs only recently allowed to assign
  // TODO: Delete outdate comment just below
  // The only instantiable parameter threads can assign are ParameterPtrs
  // And they bypass this function by calling emplace straight away

  ParameterInstance* instance;
  if (computation_type == MAIN) {
    instance = ParameterInstance::current_instance;
  } else {
    instance = parameter_function_handler->instance_;
  }
  instance->add_parameter(this);
  instantiate_internal(std::move(val), instance);
}

template <class T>
void Parameter<T>::instantiate_internal(T val, ParameterInstance* instance) {
  assert(find_instance(instance) == instances_.end() &&
         "Instantiable already instantiated");
  auto it = instances_.emplace(instance, std::move(val));
  assert(it.second);  // Just means there was no item already there

  cached_iterator_ = it.first;
}

template <class T>
bool Parameter<T>::is_set() const {
  if (fixed_instance_ == nullptr) {
    return is_set(ParameterInstance::current_instance);
  } else {
    return is_set(fixed_instance_);
  }
}

template <class T>
bool Parameter<T>::is_set(ParameterInstance* instance) const {
  auto it = find_instance(instance);
  return it != instances_.end();
}

template <class T>
const T& Parameter<T>::value() const {
  if (fixed_instance_ == nullptr) {
    return value(ParameterInstance::current_instance);
  } else {
    return value(fixed_instance_);
  }
}

template <class T>
const T& Parameter<T>::value(ParameterInstance* instance) const {
  if (computation_type == PARAMETER) {
    assert(instance == ParameterInstance::current_instance);
    return value_param();
  } else if (computation_type == MAIN) {
    return value_main(instance);
  } else {
    assert(false);
  }
}

template <class T>
const T& Parameter<T>::value_main(ParameterInstance* instance) const {
  auto it = find_instance(instance);
  assert(it != instances_.end() && "Parameter accessed but not instantiated");
  return it->second;
}

template <class T>
const T& Parameter<T>::value_param() const {
  while (true) {  // This is so we can spawn mult
    auto it = find_instance(ParameterInstance::current_instance);
    // If it's a param, should have been instantiated
    // If it's a paramptr, the API only exposes deref which makes sure it's
    // instantiated
    assert(it != instances_.end() && "Parameter accessed but not instantiated");
    if (it->first->is_child_of(parameter_function_handler->instance_)) {
      // We must spawn a lower level function
      // We stay in this loop until the coroutine dies through promotion
      parameter_function_handler->spawn_lower_instance_function(it->first);
    } else {
      return it->second;
    }
  }
}

template <class T>
typename Parameter<T>::instance_map::const_iterator Parameter<T>::find_instance(
    ParameterInstance* instance) const {
  // Check the cached iterator first
  if (cached_iterator_ != dummy_.end() &&
      (instance == cached_iterator_->first ||
       instance->is_child_of(cached_iterator_->first))) {
    return cached_iterator_;
  }

  // Traverse up the ParemeterInstance tree
  while (instance != nullptr) {
    typename Parameter<T>::instance_map::const_iterator it =
        instances_.find(instance);
    if (it != instances_.end()) {
      cached_iterator_ = it;
      return it;
    }
    instance = instance->parent();
  }

  // Could not find it
  return instances_.end();
}

template <class T>
ParameterPtr<T>::ParameterPtr()
    : Parameter<std::unique_ptr<T> >() {
  //TODO: This is wrong if the parent ParameterSpace also constructs
  // ParameterSpaces
  // PRS(Deprecated use of ParameterPtr());
  owner_ = ParameterSpace::last_constructed_parameter_space_;
  owner_->register_paramptr(this);
}

template <class T>
ParameterPtr<T>::ParameterPtr(ParameterSpace* owner)
    : Parameter<std::unique_ptr<T> >(),
      owner_(owner) {
  owner_->register_paramptr(this);
}

template <class T>
bool ParameterPtr<T>::check_value(const std::unique_ptr<T>& val) {
  assert(false);  // This should never happen
  return true;
}

template <class T>
T& ParameterPtr<T>::deref() const {
  if (this->fixed_instance_ == nullptr) {
    return deref(ParameterInstance::current_instance);
  } else {
    return deref(this->fixed_instance_);
  }
}

template <class T>
void ParameterPtr<T>::request_instantiated(ParameterInstance* instance) const {
  if (!is_set(instance)) {
    if (computation_type == PARAMETER) {
      // Our dependence must be to a different parameter function
      assert(owner_ != parameter_function_handler->owner_);
    }
    owner_->request(this);
  }
}

template <class T>
T& ParameterPtr<T>::deref(ParameterInstance* instance) const {
  request_instantiated(instance);
  return *value(instance);
}

template <class T>
bool ParameterPtr<T>::is_null(ParameterInstance* instance) const {
  request_instantiated(instance);
  return value(instance) == nullptr;
}

template <class T>
template <class U, class... Args>
void ParameterPtr<T>::new_parameter(Args&&... args) {
  // Check the parameter thread being ran is one associated with this parameter
  assert(computation_type != MAIN &&
         "Cannot assign a ParameterPtr from the main function");

  assert(!parameter_function_handler->is_calling_new());
  assert(parameter_function_handler->owner_->owns(this));
  // Are we retracing ?
  auto it = find_instance(parameter_function_handler->instance_);
  if (it == this->instances_.end()) {
    // First time instantiating

    parameter_function_handler->calling_new_ = true;
    T* new_param = new U(std::forward<Args>(args)...);
    parameter_function_handler->calling_new_ = false;

    instantiate_internal(std::unique_ptr<T>(new_param),
                         parameter_function_handler->instance_);
    parameter_function_handler->instance_->add_parameter(this);
    if (this == parameter_function_handler->ptr_request_) {
      parameter_function_handler->coroutine_yield();
    }
  } else if (it->first != parameter_function_handler->instance_) {
    // Retracing. Ideally would check a hash of the arguments
    assert(this != parameter_function_handler->ptr_request_);
    // Otherwise would have been no need to run this instance
  } else {
    assert(false && "ParameterPtr already instantiated");
  }
}

template <class T>
void ParameterPtr<T>::set_null() {
  // BAD CODE DUPLICATION with the one above

  // Check the parameter thread being ran is one associated with this parameter
  assert(computation_type != MAIN &&
         "Cannot assign a ParameterPtr from the main function");

  assert(!parameter_function_handler->is_calling_new());
  assert(parameter_function_handler->owner_->owns(this));
  // Are we retracing ?
  auto it = find_instance(parameter_function_handler->instance_);
  if (it == this->instances_.end()) {
    // First time instantiating
    instantiate_internal(std::unique_ptr<T>(nullptr),
                         parameter_function_handler->instance_);
    parameter_function_handler->instance_->add_parameter(this);
    if (this == parameter_function_handler->ptr_request_) {
      parameter_function_handler->coroutine_yield();
    }
  } else if (it->first != parameter_function_handler->instance_) {
    // Retracing. Ideally would check a hash of the arguments as well
    // But for now just check the type is correct
    assert(this != parameter_function_handler->ptr_request_);  // Otherwise
                                                               // would have
                                                               // been no need
                                                               // to run this
                                                               // instance
  } else {
    assert(false && "ParameterPtr already instantiated");
  }
}

template <class T>
CategoricalParameter<T>::CategoricalParameter(std::vector<T> values)
    : values_(std::move(values)) {}

template <class T>
bool CategoricalParameter<T>::check_value(const T& val) {
  return std::find(values_.begin(), values_.end(), val) != values_.end();
}

template <class T>
RangeParameter<T>::RangeParameter(T first, T last)
    : lower_(std::move(first)), upper_(std::move(last)) {}

template <class T>
bool RangeParameter<T>::check_value(const T& val) {
  return lower_ <= val && val <= upper_;
}

template <class T>
void RangeParameter<T>::instantiate_in_range(const T& val) {
  if (val < lower_) {
    this->assign(lower_);
  } else if (val > upper_) {
    this->assign(upper_);
  } else {
    this->assign(val);
  }
}

template <class T>
const T& RangeParameter<T>::get_lower() const {
  return lower_;
}

template <class T>
const T& RangeParameter<T>::get_upper() const {
  return upper_;
}

template <class T>
AnyParameter<T>::AnyParameter() {}

template <class T>
bool AnyParameter<T>::check_value(const T& val) {
  return true;
}
}
#endif  // PARAMETERS_T_HPP_INCLUDED
