# BOAT: A framework to build bespoke auto-tuners

Most computer systems expose many parameters or flags to users. These flags can be tuned by users to maximise the system's performance in a particular context. Auto-tuners remove this burden from users. They repeatedly evaluate flag values to find ones with good performance. Unfortunately, traditional auto-tuners often require thousands of iterations to converge. This makes them unsuitable for systems which are expensive to evaluate such as distributed systems.

**BOAT** is a framework to build auto-tuners that are specific to a particular system. This can be done by a system developer who has a good understanding of the underlying system behaviour. BOAT allows them to expose the system's behaviour in the form of a probabilistic model. The resulting auto-tuner will be able to exploit this information and quickly converge to good flag values, typically in tens of iterations.

The file [examples/branin_hoo/branin_hoo.cpp](examples/branin_hoo/branin_hoo.cpp) shows an example of this. It implements multiple auto-tuners for the [Branin-Hoo function](https://www.sfu.ca/~ssurjano/branin.html). The first auto-tuner implements a traditional Bayesian optimisation, which makes little assumptions about the shape of the function. The second auto-tuner implements a **structured Bayesian optimisation** by using an informed prior about the shape of the function.

# Dependencies

BOAT has dependencies to:
* [Boost](http://www.boost.org/)
* [Eigen](http://eigen.tuxfamily.org) 
* [NLOpt](http://ab-initio.mit.edu/wiki/index.php/NLopt)

