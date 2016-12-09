import numpy as np
import time

# --------------------------------------- CAPE ------------------------------------------------------
class Automaton:

    def __init__(self, shape, attributes, formats, parameters):
        self.parameters = parameters
        self.grid = np.zeros(shape, dtype={'names':attributes, 'formats':formats})
        self.work_grid = np.zeros(self.grid.shape, dtype={'names': attributes, 'formats': formats})

class Agent:
    def __init__(self):
        pass

# ---------------------------------------- TB -------------------------------------------------------

class TBAutomaton(Automaton):

    def __init__(self, shape, parameters):
        attributes = ['oxygen', 'chemotherapy', 'chemokine', 'contents', 'oxygen_diffusion_rate',
                      'chemotherapy_diffusion_rate']
        formats = ['float', 'float', 'float', 'object', 'float', 'float']
        Automaton.__init__(self, shape, attributes, formats, parameters)

    def diffusion(self, chemo):
        cell = self.grid[1:-1, 1:-1]
        above = self.grid[:-2, 1:-1]
        below = self.grid[2:, 1:-1]
        left = self.grid[1:-1, :-2]
        right = self.grid[1:-1, 2:]

        # oxygen
        self.work_grid['oxygen'][1:-1,1:-1] = cell['oxygen'] + self.parameters['time_step'] * \
                                        (((((cell['oxygen_diffusion_rate'] + below['oxygen_diffusion_rate'])/2) *
                                        (below['oxygen'] - cell['oxygen'])
                                        - ((cell['oxygen_diffusion_rate'] + above['oxygen_diffusion_rate'])/2) *
                                        (cell['oxygen'] - above['oxygen']))
                                        / self.parameters['spatial_step']**2) +
                                        ((((cell['oxygen_diffusion_rate'] + right['oxygen_diffusion_rate'])/2) *
                                        (right['oxygen'] - cell['oxygen'])
                                        - ((cell['oxygen_diffusion_rate'] + left['oxygen_diffusion_rate'])/2) *
                                        (cell['oxygen'] - left['oxygen']))
                                        / self.parameters['spatial_step']**2) +
                                        (self.parameters['oxygen_from_source'] * isinstance(cell, BloodVessel) *
                                         self.parameters['blood_vessel_excretion_rate']) +
                                        (self.parameters['oxygen_uptake_from_bacteria'] * cell['oxygen'] *
                                         isinstance(cell['contents'], Bacterium)))

        # chemotherapy
        if chemo:
            self.work_grid['chemotherapy'][1:-1, 1:-1] = cell['chemotherapy'] + self.parameters['time_step'] * \
                                (((((cell['chemotherapy_diffusion_rate'] + below['chemotherapy_diffusion_rate']) / 2) *
                                    (below['chemotherapy'] - cell['chemotherapy']) -
                                    ((cell['chemotherapy_diffusion_rate'] +above['chemotherapy_diffusion_rate']) / 2) *
                                    (cell['chemotherapy'] - above['chemotherapy'])) /
                                   self.parameters['spatial_step'] ** 2) + \
                                  ((((cell['chemotherapy_diffusion_rate'] + right['chemotherapy_diffusion_rate']) / 2) *
                                    (right['chemotherapy'] - cell['chemotherapy']) -
                                    ((cell['chemotherapy_diffusion_rate'] + left['chemotherapy_diffusion_rate']) / 2) *
                                    (cell['chemotherapy'] -left['chemotherapy'])) /
                                   self.parameters['spatial_step'] ** 2) +
                                  (self.parameters['chemotherapy_from_source'] * isinstance(cell, BloodVessel) *
                                   self.parameters['blood_vessel_excretion_rate']) +
                                  (self.parameters['chemotherapy_decay'] * cell['chemotherapy']))


class Bacterium(Agent):

    def __init__(self):
        Agent.__init__(self)


class BloodVessel(Agent):

    def __init__(self):
        Agent.__init__(self)


if __name__ == '__main__':

    parameters = {}
    parameters['spatial_step'] = 0.1
    parameters['time_step'] = 0.001
    parameters['oxygen_from_source'] = 0.1
    parameters['blood_vessel_excretion_rate'] = 0.1
    parameters['oxygen_uptake_from_bacteria'] = 0.1
    parameters['chemotherapy_from_source'] = 0.1
    parameters['chemotherapy_decay'] = 0.1

    shape = (100,100)
    limit = 1000
    start_time = time.time()
    tba = TBAutomaton((100, 100), parameters)
    for t in range(limit):
        print t * 0.001
        tba.diffusion(True)
    end_time = time.time()
    print "DURATION:", end_time - start_time