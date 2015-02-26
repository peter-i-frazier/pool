from abc import ABCMeta, abstractmethod, abstractproperty

class DataContainerInterface(object):

    """Container that holds original sequence data, and also offers translation
    to model readable feature matrix.
    """

    __metaclass__ = ABCMeta

    @abstractproperty
    def dim(self):
        """Return size of feature vector"""
        pass

    @abstractproperty
    def num_sample(self):
        """Return number of samples in the dataset"""
        pass

    @abstractproperty
    def feature_data(self):
        """Return data translated to feature form, is a numpy.array with size
        (num_sample, dim)

        """
        pass

    @abstractproperty
    def prior(self):
        """Return parameters for prior distribution, is a numpy.array with size
        (dim)

        """
        pass

    @prior.setter
    def prior(self, prior_in):
        """Set prior to prior_in

        :param prior_in: parameters for the new prior distribution
        :type prior_in: numpy.array of float64 with shape (dim)

        """
        pass

    @abstractmethod
    def seq_to_feature(self, sequences):
        """Translate an array of sequences to features.

        :param sequences: sequences to translate
        :type sequences: array of Sequence objects
        :return: features
        :rtype: numpy.array with shape (num_sample, dim)

        """
        pass

    @abstractmethod
    def feature_to_seq(self, features):
        """Translate features to an array of sequences.

        :param features: features to translate
        :type features: numpy.array with shape (num_sample, dim)
        :return: array of sequences
        :rtype: array of Sequence objects

        """
        pass
