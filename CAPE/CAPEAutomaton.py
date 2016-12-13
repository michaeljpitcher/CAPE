import numpy as np
import time
from collections import Counter

# --------------------------------------- CAPE ------------------------------------------------------


class Automaton:

    def __init__(self, shape, attributes, formats, time_parameters, model_parameters, initialisation):
        self.attributes = attributes
        self.model_parameters = model_parameters

        assert ('initial_time' in time_parameters.keys()), "Time parameter 'initial_time' must be defined"
        self.time = time_parameters['initial_time']
        assert ('time_step' in time_parameters.keys()), "Time parameter 'time_step' must be defined"
        self.time_step = time_parameters['time_step']
        assert ('time_limit' in time_parameters.keys()), "Time parameter 'time_limit' must be defined"
        self.time_limit = time_parameters['time_limit'] / self.time_step

        self.grid = np.zeros(shape, dtype={'names':attributes, 'formats':formats})
        self.work_grid = np.zeros(shape, dtype={'names': attributes, 'formats': formats})
        self.agents = []

        # Grid initialisation
        for attribute in initialisation.keys():
            assert attribute in self.attributes, "Invalid initialisation: key[{0}] is not valid".format(attribute)
            values = initialisation[attribute]
            for address in values:
                self.grid[address][attribute] = values[address]
                if isinstance(values[address], Agent):
                    self.agents.append(values[address])

    def run(self):

        while self.time < self.time_limit:
            self.timestep_output()
            self.update_cells()
            self.update_agents()
            self.grid, self.work_grid = self.work_grid, self.grid
            self.time += 1

    def timestep_output(self):
        print "t = ", self.time

    def update_cells(self):
        raise NotImplementedError

    def update_agents(self):
        raise NotImplementedError

    def moore_neighbours(self, address, depth):
        x,y = address
        neighbours = self.grid[x - depth:x + depth + 1, y - depth:y + depth + 1].flatten()
        neighbours = np.hstack((neighbours[:len(neighbours) // 2], neighbours[len(neighbours) // 2 + 1:]))
        return neighbours

    def von_neumann_neighbours(self, address, depth):
        x, y = address
        neighbours = None
        # TODO
        return neighbours


class Agent:
    def __init__(self):
        self.age = 0.0


class Event:
    def __init__(self):
        pass

