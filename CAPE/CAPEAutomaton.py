from Event import *
from Agent import *
import numpy as np
import itertools
import math
import csv

# --------------------------------------- CAPE ------------------------------------------------------


class Automaton:

    def __init__(self, shape, attributes, formats, time_parameters, model_parameters, output_location, values_to_record,
                 attribute_grids_to_record, initialisation):
        """
        Hybrid cellular automaton and agent-based model.
        :param shape: Shape of cellular grid (lattice)
        :param attributes: Attributes of cells in grid
        :param formats: Formats of attributes (int, float, object, etc.)
        :param time_parameters: Time specific parameters of the system
        :param model_parameters: Model parameters
        :param output_location: OS location where output from system will be written to
        :param attribute_grids_to_record: list of attributes that require recording
        :param values_to_record: Column headers for csv file of record counts
        :param initialisation: Dictionary of objects/values to be assigned to attributes at beginning of run
        """

        self.attributes = attributes
        self.model_parameters = model_parameters
        self.time_parameters = time_parameters
        self.dimensions = 2
        self.max_depth = int(self.model_parameters['max_depth'])

        # Output files
        if output_location != '':
            output_location += '/'
        self.output_location = output_location

        self.grid_files = {}
        for attribute in attribute_grids_to_record:
            assert attribute in self.attributes, "Incorrect attribute to record: {0} is not in attribute list"\
                .format(attribute)
            self.grid_files[attribute] = open(output_location + attribute + '.csv', 'w')

        self.count_file = open(output_location + 'counts.csv', 'w')
        # Initialise the count csv file with the column headers
        writer = csv.writer(self.count_file, delimiter=',')
        row = ['timestep']
        row += values_to_record
        writer.writerow(row)
        # Ensure necessary time parameters are present
        assert ('initial_time' in time_parameters.keys()), "Time parameter 'initial_time' must be defined"
        self.time = time_parameters['initial_time']
        assert ('time_step' in time_parameters.keys()), "Time parameter 'time_step' must be defined"
        self.time_step = time_parameters['time_step']
        assert ('time_limit' in time_parameters.keys()), "Time parameter 'time_limit' must be defined"
        self.time_limit = time_parameters['time_limit'] / self.time_step
        # Create the grids
        self.grid = np.zeros(shape, dtype={'names':attributes, 'formats':formats})
        self.work_grid = np.zeros(shape, dtype={'names': attributes, 'formats': formats})
        # List of agents TODO - may be redundant
        self.agents = []
        # Grid initialisation
        for attribute in initialisation.keys():
            # Make sure all attributes in initialisation are valid
            assert attribute in self.attributes, "Invalid initialisation: key[{0}] is not valid".format(attribute)
            values = initialisation[attribute]
            # Set the required values. If it's an agent, add it to the list
            for address in values:
                self.grid[address][attribute] = values[address]
                if isinstance(values[address], Agent):
                    self.agents.append(values[address])

        # NEIGHBOURHOODS
        # Builds dictionaries of neighbour cells for use with neighbour functions
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

        # Event lists
        self.potential_events = []
        self.acceptable_events = []

    def close_files(self):
        """
        Close all files
        :return:
        """
        for grid_file in self.grid_files.values():
            grid_file.close()
        self.count_file.close()

    def run(self):
        """
        Run the automaton from start time to time limit
        :return:
        """
        # Record the initial state
        self.record()
        while self.time < self.time_limit:

            # Increment time
            self.time += 1

            # Output to console
            self.timestep_output()
            # Run cellular automaton update
            self.update_cells()

            self.potential_events = []
            self.acceptable_events = []

            # Run agent-based model update
            new_events = self.generate_events_from_agents()
            self.potential_events += new_events

            self.acceptable_events = self.conflict_resolve_events()
            #
            # self.perform_events()

            # Set the main grid
            self.grid = self.work_grid.copy()

            # Record if necessary
            self.record()

    def record(self):
        # Recording
        if self.time % self.time_parameters['interval_to_record_grid'] == 0:
            self.record_grids()
        if self.time % self.time_parameters['interval_to_record_counts'] == 0:
            self.record_counts()

    def timestep_output(self):
        """
        Output to the console
        :return:
        """
        print "t = ", self.time * self.time_step

    def update_cells(self):
        """
        Cellular automaton update. Specific rules must be overriden by subclass
        :return:
        """
        raise NotImplementedError

    def generate_events_from_agents(self):
        """
        Agent-based model update. Specific rules must be overriden by subclass
        :return:
        """
        raise NotImplementedError

    def is_on_grid(self, address):
        """
        Check if a given address is actually on the grid
        :param address:
        :return:
        """
        for i in range(len(address)):
            if address[i] < 0 or address[i] >= self.grid.shape[i]:
                return False
        return True

    def neighbours(self, address, depth, type='moore'):
        """
        Get the neighbours for a given neighbourhood type and depth
        :param address: The address which requires neighbours
        :param depth: The depth to search for
        :param type: The type of neighbourhood (moore or von_neumann)
        :return:
        """
        neighbours = {}
        # Check type
        if type == 'moore':
            relative_addresses = self.moore_relative[depth]
        elif type == 'von_neumann':
            relative_addresses = self.von_neumann_relative[depth]
        else:
            raise Exception, "Invalid neighbourhood type"
        # For each relative address, apply the coordinates to address to get new addresses
        for n in relative_addresses:
            neighbour_address = tuple([address[i] + n[i] for i in range(len(address))])
            if self.is_on_grid(neighbour_address):
                neighbours[neighbour_address] = self.grid[neighbour_address]
        return neighbours

    def moore_neighbours(self, address, depth):
        """
        Get all moore neighbours (N, E, S, W, NE, NW, SE, SW) to a given address
        :param address: Address of cell which requires neighbours
        :param depth: Depth to search for
        :return:
        """
        return self.neighbours(address, depth, 'moore')

    def von_neumann_neighbours(self, address, depth):
        """
        Get all von Neumann neighbours (N, E, S, W) to a given address
        :param address: Address of cell which requires neighbours
        :param depth: Depth to search for
        :return:
        """
        return self.neighbours(address, depth, 'von_neumann')

    def record_grids(self):
        """
        Write the contents of the grid to the output file (based on specified agent codes)
        :return:
        """
        writers= {}
        for a in self.grid_files:
            writers[a] = csv.writer(self.grid_files[a], delimiter=',')
        # Loop through every cell
        for i in range(self.grid.shape[0]):
            rows = {}
            for a in writers:
                rows[a] = []
            for j in range(self.grid.shape[1]):
                # Write the value (or code if it's an agent)
                for a in writers:
                    value = self.grid[(i,j)][a]
                    if isinstance(value, Agent):
                        rows[a].append(value.output_code())
                    else:
                        rows[a].append(value)
            for a in writers:
                writers[a].writerow(rows[a])

    def record_counts(self):
        """
        Records the counts of various items. Must be overriden by subclass.
        :return:
        """
        raise NotImplementedError

    def conflict_resolve_events(self):
        """
        Picks potential events in random order. Evaluates acceptability based on dependent addresses - if any dependent
        address for an event has already been processed, event is discarded.
        Can be overridden if a different resolution method is required (e.g. using a priority system)
        :return:
        """
        processed_addresses = []
        acceptable_events = []

        np.random.shuffle(self.potential_events)

        while len(self.potential_events) > 0:
            event = self.potential_events.pop()
            acceptable = True

            for address in event.dependent_addresses:
                if address in processed_addresses:
                    # Discard event as conflicts with a previous event
                    acceptable = False
                    break

            if acceptable:
                amended_impacted_addresses = []
                for address in event.impacted_addresses:
                    if address not in processed_addresses:
                        amended_impacted_addresses.append(address)
                        processed_addresses.append(address)
                event.impacted_addresses = amended_impacted_addresses
                acceptable_events.append(event)

        return acceptable_events

    def perform_events(self):
        raise NotImplementedError
