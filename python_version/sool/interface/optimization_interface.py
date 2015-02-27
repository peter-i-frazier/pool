# -*- coding: utf-8 -*-
"""Interfaces for optimization (MINIMIZATION):
1. OptimizableInterface: things that can be optimized (objective function)
2. OptimizerInterface: perform the optimization.

"""
from abc import ABCMeta, abstractmethod, abstractproperty


class OptimizableInterface(object):

    """Interface that an object must fulfill to be optimized by an implementation
    of OptimizerInterface.

    Think of it as objective function f(x), and OptimizerInterface implementation
    finds maximum value of f(x). Then we need OptimizableInterface to
    a. compute problem size (how many independent parameters to optimize over, i.e
       size of vector x)
    b. compute objective value f(x) given x
    c. compute special form of f(x): if f(x) is convex, then we can use convex
       optimization solver to find min f(x). This requires that we write f(x)
       in a specific form, and this function returns matrices, vectors that
       construct such f(x)

    """

    __metaclass__ = ABCMeta

    @abstractproperty
    def problem_size(self):
        """Return the number of independent parameters to optimize."""
        pass

    @abstractproperty
    def current_point(self):
        """Return the current_point (array of float64 with shape (problem_size))
        """
        pass

    @current_point.setter
    def current_point(self, current_point):
        """Set current_point to the specified point; ordering must match.

        :param current_point: current_point at which to evaluate the objective
                              function
        :type current_point: array of float64 with shape (problem_size)

        """
        pass

    @abstractproperty
    def special_form(self):
        """Return special form that constructs objective function."""
        pass

    @abstractmethod
    def compute_objective_function(self, **kwargs):
        """Compute f(current_point).

        :return: value of objective function evaluated at ``current_point``
        :rtype: float64

        """
        pass

class OptimizerInterface(object):

    """Interface to minimize any object implementing OptimizableInterface."""

    __metaclass__ = ABCMeta

    @abstractproperty
    def performance_value(self):
        pass
    
    @abstractmethod
    def optimize(self, **kwargs):
        """Minimize a function f(x), represented by an implementation of
        OptimizableInterface.

        """
        pass