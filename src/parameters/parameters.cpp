#include "parameters.hpp"

// Static Data
namespace boat{
std::unique_ptr<ParameterInstance> ParameterInstance::root_parameter_instance_(
    new ParameterInstance());
ParameterInstance* ParameterInstance::current_instance(
  root_parameter_instance_.get());
computation_type_t computation_type = MAIN;
ParameterFunction* parameter_function_handler;
ParameterSpace* ParameterSpace::last_constructed_parameter_space_;

// local helper functions

ParameterInstance::ParameterInstance()
    : parent_(nullptr), current_child_(nullptr) {}

ParameterInstance::ParameterInstance(ParameterInstance* parent)
    : parent_(parent), current_child_(nullptr) {
  current_instance = this;
  // Should make everything current: allow for values associated with
  // instances of the parents only
}

void ParameterInstance::promote(ParameterInstance* promoted_child) {
  // No lock threads should be blocked
  assert(computation_type == MAIN);

  // this becomes the current instance
  assert(current_instance == children_.back().get());
  current_instance = this;
  // Other children die
  // The child must have no children

  for (auto& child : children_) {
    assert(!child->has_children());
    if (child.get() == promoted_child) {
      std::swap(child, *children_.begin());
    }
  }
  assert(children_.begin()->get() == promoted_child);

  children_.resize(1);  // Will destruct all other children instances

  // We promote the instantiables intantiated by the child
  size_t old_size = promoted_child->parameters_.size();
  while (!promoted_child->parameters_.empty()) {
    ParameterInterface* parameter = *promoted_child->parameters_.begin();
    parameter->promote_instance(this, promoted_child);

    assert(old_size > promoted_child->parameters_.size());
    old_size = promoted_child->parameters_.size();
  }

  // We promote the parameter functions (which were genereated by reading
  // these instantiables)
  old_size = promoted_child->parameter_spaces_.size();  // Debug
  while (promoted_child->parameter_spaces_.size() > 0) {
    (*promoted_child->parameter_spaces_.begin())
        ->promote_instance(this, promoted_child);
    assert(old_size > promoted_child->parameter_spaces_.size());
    old_size = promoted_child->parameter_spaces_.size();
  }

  assert(children_.size() == 1);
  assert(children_.begin()->get() == promoted_child);
  assert(promoted_child->parameters_.empty());
  assert(promoted_child->parameter_spaces_.empty());
  // Finally delete the promoted child, this should yield no further destruction
  children_.clear();
}

void ParameterInstance::remove_children() {
  // Function used to simply remove all children instances without promoting any
  assert(computation_type == MAIN);

  // this becomes the current instance
  assert(current_instance == children_.back().get());
  current_instance = this;

  children_.clear();
}

void ParameterInstance::check_structure() {
  // Can only run on main branch
  for (auto& child : children_) {
    if (child.get() == current_child_) {
      child->check_structure();
    } else {
      assert(!child->has_children());
    }
  }
}

ParameterInstance* ParameterInstance::parent() { return parent_; }

bool ParameterInstance::is_child_of(ParameterInstance* instance) {
  ParameterInstance* to_test = parent_;
  // We go up the tree
  while (to_test != nullptr) {
    if (instance == to_test) {
      return true;
    }
    to_test = to_test->parent_;
  }
  return false;
}

bool ParameterInstance::has_children() { return !children_.empty(); }

void ParameterInstance::add_parameter_space(ParameterSpace* parameter) {
  parameter_spaces_.insert(parameter);
}
void ParameterInstance::remove_parameter_space(
    ParameterSpace* parameter_space) {
  assert(computation_type == MAIN);

  auto it = parameter_spaces_.find(parameter_space);
  assert(it != parameter_spaces_.end());
  parameter_spaces_.erase(it);
}

bool ParameterInstance::has_parameter_space(ParameterSpace* parameter_space) {
  return parameter_spaces_.find(parameter_space) != parameter_spaces_.end();
}

bool ParameterInstance::has_parameter(ParameterInterface* parameter) {
  return parameters_.find(parameter) != parameters_.end();
}

ParameterInstance::~ParameterInstance() {
  assert(computation_type == MAIN);

  size_t old_size = parameter_spaces_.size();
  while (!parameter_spaces_.empty()) {
    (*parameter_spaces_.begin())->remove_parameter_function_instance(this);

    // Check it worked
    assert(parameter_spaces_.size() < old_size);
    old_size = parameter_spaces_.size();
  }

  // old_size for debug
  // delete the parameters, this may include some parameter_ptrs
  // which will yield futher destructions
  old_size = parameters_.size();
  while (!parameters_.empty()) {
    (*parameters_.begin())->remove_instance(this);

    // Check it worked
    assert(parameters_.size() < old_size);
    old_size = parameters_.size();
  }
}

ParameterInstance* ParameterInstance::create_child() {
  // No lock threads should be blocked
  assert(computation_type == MAIN);

  for (auto& child : children_) {
    assert(!child->has_children());
  }
  children_.emplace_back(new ParameterInstance(this));
  current_child_ = children_.back().get();
  return children_.back().get();
}

void ParameterInstance::add_parameter(ParameterInterface* parameter) {
  parameters_.insert(parameter);
}

void ParameterInstance::remove_parameter(ParameterInterface* parameter) {
  assert(computation_type == MAIN);
  parameters_.erase(parameter);
}

void ParameterInstance::set_fixed() {
  if (parent_ != nullptr) {
    parent_->set_fixed();
  }
  for (auto& param : parameters_) {
    param->set_fixed();
  }
}

void ParameterInstance::unset_fixed() {
  // Calling on the parent first is important so wehn the function returns,
  // the bottom most instance is current_instance
  if (parent_ != nullptr) {
    parent_->unset_fixed();
  }
  for (auto& param : parameters_) {
    param->unset_fixed();
  }
  current_instance = this;
}

// ParameterFunction
ParameterFunction::ParameterFunction(ParameterInstance* instance,
                                     ParameterSpace* owner)
    : instance_(instance),
      owner_(owner),
      parameter_coroutine_([&](YieldCoroutine& coro) {
        yield_coroutine_ = &coro;
        owner_->parameter_function();
      }),
      calling_function_(nullptr),
      ptr_request_(nullptr),
      calling_new_(false) {
  instance_->add_parameter_space(owner_);
}

ParameterFunction::~ParameterFunction() {
  assert(computation_type == MAIN);

  // Notify the instance we have died
  assert(instance_ != nullptr);  // TODO simplify this after debug
  if (instance_ != nullptr) {
    instance_->remove_parameter_space(owner_);
  }
}

void ParameterFunction::spawn_lower_instance_function(
    ParameterInstance* instance) {
  assert(instance->is_child_of(instance_));
  owner_->add_parameter_function_instance(instance);
  coroutine_yield();
}

void ParameterFunction::promote(ParameterInstance* parent) {
  // Called when there is no parameter function already on the parent instance

  assert(parent == instance_->parent());

  // We remove ourselves from the child instance (which is about to be emptied)
  // and add ourselves to the parent
  assert(instance_->has_parameter_space(owner_));
  assert(!parent->has_parameter_space(owner_));

  instance_->remove_parameter_space(owner_);
  instance_ = parent;
  instance_->add_parameter_space(owner_);
}

void ParameterFunction::promote_from(ParameterFunction& parent_function,
                                     ParameterFunction& child_function) {
  ParameterInstance* child_instance = child_function.instance_;
  ParameterInstance* parent_instance = parent_function.instance_;

  // Check the parent has the correct properties
  assert(parent_instance == child_instance->parent());
  assert(parent_function.owner_ == child_function.owner_);

  // Check both instances know about this parameter
  assert(child_instance->has_parameter_space(child_function.owner_));
  assert(parent_instance->has_parameter_space(parent_function.owner_));

  // Finally give the child function the correct instance
  // the parent function will delete the parameter from the instance
  // once it is deleted, which should be just after
  std::swap(child_function.instance_, parent_function.instance_);
}

bool ParameterFunction::is_calling_new() { return calling_new_; }

void ParameterFunction::execute_up_to(const ParameterInterface* parameter_ptr) {
  assert(ptr_request_ ==
         nullptr);  // The coroutine isn't currently part of the "stack"
  if (computation_type == MAIN) {
    calling_function_ = nullptr;
  } else {
    calling_function_ = parameter_function_handler;
  }
  computation_type = PARAMETER;
  parameter_function_handler = this;
  ptr_request_ = parameter_ptr;
  parameter_coroutine_();
}

void ParameterFunction::coroutine_yield() {
  // We just instantiated the requested ptr
  assert(computation_type == PARAMETER);
  assert(ptr_request_ != nullptr);

  ptr_request_ = nullptr;
  if (calling_function_ == nullptr) {
    // We were called by MAIN
    computation_type = MAIN;
    parameter_function_handler = nullptr;
    (*yield_coroutine_)();
  } else {
    // We were called by another PARAMETRER
    parameter_function_handler = calling_function_;
    (*yield_coroutine_)(calling_function_->parameter_coroutine_);
  }

  // When we resume here we've been called again
}

// ParameterSpace

ParameterSpace::ParameterSpace() {
  if (computation_type == MAIN) {
    add_parameter_function_instance(ParameterInstance::current_instance);
  } else {
    assert(computation_type == PARAMETER);
    assert(parameter_function_handler->is_calling_new());
    add_parameter_function_instance(parameter_function_handler->instance_);
  }
  last_constructed_parameter_space_ = this;
}

void ParameterSpace::add_parameter_function_instance(
    ParameterInstance* instance) {
  parameter_functions_.emplace(instance,
                               std::unique_ptr<ParameterFunction>(
                                   new ParameterFunction(instance, this)));
}

void ParameterSpace::promote_instance(ParameterInstance* parent,
                                      ParameterInstance* child) {
  // No lock, threads should be blocked
  assert(parent == child->parent());
  assert(child->has_parameter_space(this));
  auto child_function = parameter_functions_.find(child);
  assert(child_function != parameter_functions_.end());
  auto parent_function = parameter_functions_.find(parent);

  // Did we already have a function for that instance ?
  if (parent_function == parameter_functions_.end()) {
    // Either:
    //  -This is the only parameter function, OR
    //  -The parent is at a higher level
    assert(!parent->has_parameter_space(this));

    child_function->second->promote(parent);

    parameter_functions_.emplace(parent, std::move(child_function->second));
    parameter_functions_.erase(child_function);
  } else {
    // There is a parent function here as well
    assert(parent->has_parameter_space(this));
    // We make the child acquire the ownership of the parameters spawed by its
    // parent
    ParameterFunction::promote_from(*parent_function->second.get(),
                                    *child_function->second.get());

    // We then kill the parent and replace it with the child
    std::swap(parent_function->second, child_function->second);
    parameter_functions_.erase(child_function);
  }
}

void ParameterSpace::remove_parameter_function_instance(
    ParameterInstance* instance) {
  auto it = parameter_functions_.find(instance);
  assert(it != parameter_functions_.end());
  parameter_functions_.erase(it);
}

void ParameterSpace::register_paramptr(ParameterInterface* param_ptr) {
  owned_paramptrs_.push_back(param_ptr);
}

bool ParameterSpace::owns(const ParameterInterface* param_ptr) const {
  return std::find(owned_paramptrs_.begin(), owned_paramptrs_.end(),
                   param_ptr) != owned_paramptrs_.end();
}

void ParameterSpace::request(const ParameterInterface* parameter_ptr) {
  assert(owns(parameter_ptr));
  auto prev_it = parameter_functions_.end();
  while (parameter_ptr->is_instantiated(ParameterInstance::current_instance) ==
         nullptr) {
    auto it = find_lowest_level_param_func(ParameterInstance::current_instance);

    // Check we're making progress
    assert(it != prev_it);
    prev_it = it;

    it->second->execute_up_to(parameter_ptr);
    // The paramfunc creates a new paramfunc if it is
  }
}

typename ParameterSpace::parameter_function_map::iterator
ParameterSpace::find_lowest_level_param_func(ParameterInstance* inst) {
  parameter_function_map::iterator it = parameter_functions_.end();
  while (it == parameter_functions_.end()) {
    assert(inst != nullptr);
    it = parameter_functions_.find(inst);
    inst = inst->parent();
  }
  return it;
}

}
