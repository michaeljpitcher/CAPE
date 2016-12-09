import numpy as np


class Automaton:

    def __init__(self, shape, attributes, formats):

        self.grid = np.zeros(shape, dtype={'names':attributes, 'formats':formats})


if __name__ == '__main__':
    atts = ['oxygen', 'chemotherapy', 'chemokine', 'contents']
    formats = ['float', 'float', 'float', 'object']
    a = Automaton((10,10), atts, formats)
    print a.grid