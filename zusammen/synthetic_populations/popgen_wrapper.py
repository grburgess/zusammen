import abc


class PopGenWrapper(object, metaclass=abc.ABCMeta):
    def __init__(self, name):

        self._name = name

    @abc.abstractmethod
    def _build_population_generator(self):

        pass

    def generate_population(self):

        pass

    @property
    def population(self):
        return self._population
