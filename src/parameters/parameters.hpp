#ifndef PARAMETERS_HPP_INCLUDED
#define PARAMETERS_HPP_INCLUDED

#include <mutex>
#include <thread>
#include <condition_variable>
#include <functional>
#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <assert.h>
#include <algorithm>
#include <memory>
#define BOOST_COROUTINES_NO_DEPRECATION_WARNING
#include <boost/coroutine/all.hpp>

namespace boat{
class ParameterFunction;
class ParameterSpace;
class ParameterInterface;
typedef enum { MAIN, PARAMETER } computation_type_t;
extern computation_type_t computation_type;
extern ParameterFunction* parameter_function_handler;

template <class T>
class TopParameterPtr;

typedef boost::coroutines::symmetric_coroutine<void>::call_type
    ParameterCoroutine;
typedef boost::coroutines::symmetric_coroutine<void>::yield_type YieldCoroutine;

// Main classes
class ParameterInstance {
  template <class T>
  friend class TopParameterPtr;

 public:
  ~ParameterInstance();

  // Modifiers
  ParameterInstance* create_child();
  void promote(ParameterInstance* child);
  void remove_children();
  void set_fixed();
  void unset_fixed();

  // Queries
  ParameterInstance* parent();
  void check_structure();
  bool is_child_of(ParameterInstance* instance);
  bool has_children();

  bool has_parameter_space(ParameterSpace* parameter);
  bool has_parameter(ParameterInterface* parameter);

  // New ones

  void add_parameter_space(ParameterSpace* parameter_space);
  void remove_parameter_space(ParameterSpace* parameter);

  void add_parameter(ParameterInterface* parameter);
  void remove_parameter(ParameterInterface* parameter);

  static ParameterInstance* current_instance;
  static std::unique_ptr<ParameterInstance> root_parameter_instance_;

  // private:
  ParameterInstance();
  ParameterInstance(ParameterInstance* parent);

  // The tree structure fields
  ParameterInstance* parent_;
  std::vector<std::unique_ptr<ParameterInstance>> children_;
  ParameterInstance* current_child_;

  // The parameters (including ParameterPtrs) which have been instantiated with
  // this instance
  std::unordered_set<ParameterInterface*> parameters_;

  // Parameters which hold a functions of this instance
  // The parent parameters are before in the ordering
  std::unordered_set<ParameterSpace*> parameter_spaces_;

  // Used to know whether all functions have blocked
};

template <class T>
class TopParameterPtr {
 public:
  template <class... Args>
  TopParameterPtr(Args&&... args);
  ~TopParameterPtr();

  // Sets the root instance fixed, and the previous instance current
  void set_previous_instance_current();

  T& deref() { return *param_; }

 private:
  ParameterInstance* init_previous();
  std::unique_ptr<ParameterInstance> root_instance_;
  ParameterInstance* previous_instance_;
  std::unique_ptr<T> param_;
};

class ParameterSpace {
  friend class ParameterFunction;
  friend class ParameterInterface;
  friend class ParameterInstance;
  template <class T>
  friend class TopParameterPtr;
  template <class T>
  friend class ParameterPtr;
  template <class T>
  friend class Parameter;
  friend void spawn_parameter_space_funcs(std::vector<ParameterSpace*>& in,
                                          ParameterInstance* instance);

 public:
  typedef std::unordered_map<ParameterInstance*,
                             std::unique_ptr<ParameterFunction>>
      parameter_function_map;
  ParameterSpace();
  virtual ~ParameterSpace() {}

 protected:
  // The user implements parameter_function(). If they need to spawn some new
  // parameters, they use the new_parameter function.
  virtual void parameter_function() = 0;

 private:
  void promote_instance(ParameterInstance* parent, ParameterInstance* child);
  void add_parameter_function_instance(ParameterInstance* instance);
  void remove_parameter_function_instance(ParameterInstance* instance);

  typename ParameterSpace::parameter_function_map::iterator
  find_lowest_level_param_func(ParameterInstance* inst);

  void request(const ParameterInterface* parameter_ptr);
  void register_paramptr(ParameterInterface* param_ptr);
  bool owns(const ParameterInterface* param_ptr) const;

  parameter_function_map parameter_functions_;
  std::vector<ParameterInterface*> owned_paramptrs_;  // Only used fo debug

  static bool constructing_top_;  // Equivalent of the same fields for param
                                  // funcs but for the top parameter
  static ParameterSpace* last_constructed_parameter_space_;
};

class ParameterFunction {
  friend void ParameterSpace::promote_instance(ParameterInstance* parent,
                                               ParameterInstance* child);
  friend ParameterInstance;
  template <class T>
  friend class ParameterPtr;

 public:
  ParameterFunction(ParameterInstance* instance, ParameterSpace* owner);
  ParameterFunction(const ParameterFunction& other,
                    ParameterInstance* instance);
  ~ParameterFunction();

  void spawn_lower_instance_function(ParameterInstance* instance);
  void execute_up_to(const ParameterInterface* parameter_ptr);
  bool is_calling_new();
  void set_fixed();
  void unset_fixed();

  ParameterInstance* instance_;

 private:
  static void promote_from(ParameterFunction& parent_function,
                           ParameterFunction& child_function);
  void promote(ParameterInstance* parent);

  void coroutine_yield();

  ParameterSpace* const owner_;

  ParameterCoroutine parameter_coroutine_;
  YieldCoroutine* yield_coroutine_;
  ParameterFunction* calling_function_;
  // The pointer for which we're trying to get a value
  const ParameterInterface* ptr_request_;

  bool calling_new_;  // Whether we're currently constructing a new object
};

class ParameterInterface {
 public:
  virtual void remove_instance(ParameterInstance* instance) = 0;
  virtual void promote_instance(ParameterInstance* parent,
                                ParameterInstance* child) = 0;
  // Return nullptr if not instantiated, the corresponding instance otherwise
  //(Either inst or a parent)
  virtual ParameterInstance* is_instantiated(ParameterInstance* inst) const = 0;

  virtual void set_fixed() = 0;
  virtual void unset_fixed() = 0;
};

template <class T>
class Parameter : public ParameterInterface {
  template <class U>
  friend class ParameterPtr;
  friend class ParameterSpace;

 public:
  typedef std::unordered_map<ParameterInstance*, T> instance_map;
  Parameter();
  ~Parameter();
  Parameter(const Parameter&) = delete;

  const T& value() const;
  const T& value(ParameterInstance* instance) const;
  bool is_set() const;
  bool is_set(ParameterInstance* instance) const;
  void assign(T val);
  typename instance_map::const_iterator find_instance(
      ParameterInstance* instance) const;
  void remove_instance(ParameterInstance* instance) override;
  void promote_instance(ParameterInstance* parent,
                        ParameterInstance* child) override;
  ParameterInstance* is_instantiated(ParameterInstance* inst) const override;
  void set_fixed() override;
  void unset_fixed() override;

  virtual bool check_value(const T& val) = 0;

  // protected:
  // private:
  void remove_instance_internal(ParameterInstance* instance, bool deteting);
  void instantiate_internal(T arg, ParameterInstance* instance);
  const T& value_main(ParameterInstance* instance) const;
  const T& value_param() const;

  instance_map instances_;

  static const instance_map dummy_;  // Used so we can have an end() iterator
                                     // which definitely doesn't get invalidated

  ParameterInstance* fixed_instance_;
  mutable typename instance_map::const_iterator cached_iterator_;  // New, needs
                                                                   // checking
};

template <class T>
class ParameterPtr final : protected Parameter<std::unique_ptr<T>> {
  using Parameter<std::unique_ptr<T>>::value;
  using Parameter<std::unique_ptr<T>>::find_instance;
  using Parameter<std::unique_ptr<T>>::instantiate_internal;

 public:
  using Parameter<std::unique_ptr<T>>::is_set;
  ParameterPtr();  //Deprecated
  ParameterPtr(ParameterSpace* owner);
  bool check_value(const std::unique_ptr<T>& val) override;
  T& deref() const;
  T& deref(ParameterInstance* instance) const;
  bool is_null(
      ParameterInstance* instance = ParameterInstance::current_instance) const;
  template <class U = T, class... Args>
  void new_parameter(Args&&... args);
  void set_null();

 //private:
  void request_instantiated(ParameterInstance* instance) const;

  ParameterSpace* owner_;
};

// Useful instances
class DummyInstance {
 public:
  DummyInstance() : instance_(ParameterInstance::current_instance) {
    instance_->create_child();
  }
  ~DummyInstance() {
    assert(instance_ == ParameterInstance::current_instance->parent());
    instance_->remove_children();
  }

 private:
  ParameterInstance* instance_;
};

class InstanceView {
  public:
  InstanceView(ParameterInstance* instance)
      : true_current_(ParameterInstance::current_instance),
        other_(instance) {
    set();
  }

  ~InstanceView() {
    unset();
  }

  void set(){
    ParameterInstance::current_instance = other_;
  }

  void unset() {
    ParameterInstance::current_instance = true_current_;
  }
 private:
  ParameterInstance* true_current_;
  ParameterInstance* other_;
};

// Useful Parameters
template <class T>
class CategoricalParameter : public Parameter<T> {
 public:
  CategoricalParameter(std::vector<T> values);
  bool check_value(const T& val) override;

 private:
  std::vector<T> values_;
};

template <class T>
class RangeParameter : public Parameter<T> {
 public:
  RangeParameter(T lower, T upper);
  bool check_value(const T& val) override;
  const T& get_lower() const;
  const T& get_upper() const;
  void instantiate_in_range(const T& val);

 private:
  T lower_;
  T upper_;
};

template <class T>
class AnyParameter : public Parameter<T> {
 public:
  AnyParameter();
  bool check_value(const T& val) override;
};

class BoolParameter : public CategoricalParameter<bool> {
 public:
  BoolParameter(std::vector<bool> v = {true, false})
      : CategoricalParameter(v) {}
};

}
// template implementations
#include "parameters_t.hpp"

#endif  // PARAMETERS_HPP_INCLUDED
