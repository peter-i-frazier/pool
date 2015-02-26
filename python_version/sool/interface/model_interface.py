from abc import ABCMeta, abstractmethod, abstractproperty

class ModelParameterInterface(object):
    """Interface of tunable parameters for a model"""

    __metaclass__ = ABCMeta

    @abstractproperty
    def param(self):
        """Return tunable parameters for a model"""
        pass

    @param.setter
    def param(self, **kwargs):
        pass

class ModelInterface(object):
    """Interface of statistical model"""

    __metaclass__ = ABCMeta

    @abstractmethod
    def train(self, **kwargs):
        """Train the model given dataset"""
        pass

    @abstractmethod
    def predict(self, features, **kwargs):
        """Return predicted values of variates in feature form"""
        pass

    @abstractmethod
    def cross_validation(self, num_folds):
        pass

