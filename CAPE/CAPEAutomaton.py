import numpy as np
import time
from collections import Counter
import itertools
import math
import csv

# --------------------------------------- CAPE ------------------------------------------------------


class Automaton:

    def __init__(self, shape, attributes, formats, time_parameters, model_parameters, output_location, values_to_record,
                 initialisation):
        self.attributes = attributes
        self.model_parameters = model_parameters
        self.time_parameters = time_parameters
        self.dimensions = 2
        self.max_depth = 3

        # Output files
        if output_location != '':
            output_location += '/'
        self.output_location = output_location
        self.grid_file = open(output_location + 'grid.csv', 'w')
        self.count_file = open(output_location + 'counts.csv', 'w')

        writer = csv.writer(self.count_file, delimiter=',')
        row = ['timestep']
        row += values_to_record
        writer.writerow(row)

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

        # Initialise empty dictionaries
        self.moore_relative = dict()
        self.von_neumann_relative = dict()
        # Add an entry for each depth
        for d in range(1, self.max_depth + 1):
            self.von_neumann_relative[d] = []

        for depth in range(1, self.max_depth + 1):
            # Get truth table values (e.g. depth 2 gives [-2,-1,0,1,2] for range_)
            range_ = range(-depth, depth + 1)
            # Use product to find all combinations for given depth and number of dimensions
            row = list(itertools.product(range_, repeat=self.dimensions))
            # Remove the 0 entry e.g. (0,0) for 2 dimensions
            row.remove((0,) * self.dimensions)

            reduced_row_moore = []
            self.von_neumann_relative[depth] = []
            for neighbour in row:
                # Calculate Manhattan distance and add to appropriate von Neumann table row
                manhattan_distance = int(sum([math.fabs(x) for x in neighbour]))
                if manhattan_distance <= self.max_depth and neighbour not in \
                        self.von_neumann_relative[manhattan_distance]:
                    self.von_neumann_relative[manhattan_distance].append(neighbour)
                # Check if one of coordinates = depth, if so then use for moore at this depth
                for x in neighbour:
                    if int(math.fabs(x)) == depth:
                        reduced_row_moore.append(neighbour)
                        break
            self.moore_relative[depth] = reduced_row_moore

    def run(self):
        while self.time < self.time_limit:
            # Output to console
            self.timestep_output()
            # Run cellular automaton update
            self.update_cells()
            # Run agent-based model update
            self.update_agents()
            # Swap grids
            self.grid, self.work_grid = self.work_grid, self.grid
            # Recording
            if self.time % self.time_parameters['interval_to_record_grid'] == 0:
                self.record_grid()
            if self.time % self.time_parameters['interval_to_record_counts'] == 0:
                self.record_counts()
            # Increment time
            self.time += 1

    def timestep_output(self):
        print "t = ", self.time * self.time_step

    def update_cells(self):
        raise NotImplementedError

    def update_agents(self):
        raise NotImplementedError

    def is_on_grid(self, address):
        for i in range(len(address)):
            if address[i] < 0 or address[i] >= self.grid.shape[i]:
                return False
        return True

    def neighbours(self, address, depth, type='moore'):
        neighbours = {}
        if type == 'moore':
            relative_addresses = self.moore_relative[depth]
        elif type == 'von_neumann':
            relative_addresses = self.von_neumann_relative[depth]
        else:
            raise Exception, "Invalid neighbourhood type"
        for n in relative_addresses:
            neighbour_address = tuple([address[i] + n[i] for i in range(len(address))])
            if self.is_on_grid(neighbour_address):
                neighbours[neighbour_address] = self.grid[neighbour_address]
        return neighbours

    def moore_neighbours(self, address, depth):
        return self.neighbours(address, depth, 'moore')

    def von_neumann_neighbours(self, address, depth):
        return self.neighbours(address, depth, 'von_neumann')

    def record_grid(self):
        writer = csv.writer(self.grid_file, delimiter=',')
        for i in range(self.grid.shape[0]):
            row = []
            for j in range(self.grid.shape[1]):
                contents = self.grid[(i,j)]['contents']
                if isinstance(contents, Agent):
                    row.append(contents.output_code())
                else:
                    row.append(0.0)
            writer.writerow(row)

    def record_counts(self):
        raise NotImplementedError


class Agent:
    def __init__(self):
        self.age = 0.0

    def output_code(self):
        raise NotImplementedError


class Event:
    def __init__(self):
        pass

