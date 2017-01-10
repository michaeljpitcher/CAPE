from CAPE.CAPEAutomaton import *
from TBAgents import *
from TBEvents import *
from collections import Counter
import cProfile


class TBAutomaton(Automaton):

    def __init__(self, shape, time_parameters, model_parameters, output_location,
                 blood_vessel_addresses, initial_macrophage_addresses,
                 initial_fast_bacteria_addresses, initial_slow_bacteria_addresses):
        """
        Specific model of CAPE Automaton to investigate TB infection. Grid is square of alveolar tissue, agents are
        bacteria and immune cells that act upon the tissue. Cellular automaton handles diffusion of oxygen,
        chemotherapy and chemokine
        :param shape: Shape of grid
        :param time_parameters: Time specific parameters
        :param model_parameters: Model parameters
        :param output_location: Location output files will be written to
        :param blood_vessel_addresses: Addresses to place blood vessels
        :param initial_macrophage_addresses: Addresses to place macrophages
        :param initial_fast_bacteria_addresses: Addresses to place fast bacteria
        :param initial_slow_bacteria_addresses: Addresses to place slow bacteria
        """
        # Hard-coded attributes and formats
        attributes = ['oxygen', 'chemotherapy', 'chemokine', 'contents', 'oxygen_diffusion_rate',
                      'chemotherapy_diffusion_rate', 'blood_vessel']
        formats = ['float', 'float', 'float', 'object', 'float', 'float', 'float']

        # Initialise list (blood vessels never change)
        self.blood_vessel_addresses = blood_vessel_addresses
        self.macrophages = []
        self.bacteria = []
        self.t_cells = []
        self.caseum_addresses = []

        # INITIALISE
        initialisation = {}
        initialisation['contents'] = {}
        initialisation['oxygen'] = {}
        initialisation['oxygen_diffusion_rate'] = {}
        initialisation['chemotherapy_diffusion_rate'] = {}
        initialisation['blood_vessel'] = {}
        # Blood vessels & oxygen
        self.blood_vessel_addresses = blood_vessel_addresses
        for bva in blood_vessel_addresses:
            # Not strictly necessary, but amend the contents to show the cell is not empty
            initialisation['contents'][bva] = model_parameters['blood_vessel_value']
            initialisation['blood_vessel'][bva] = model_parameters['blood_vessel_value']
            initialisation['oxygen'][bva] = model_parameters['blood_vessel_value'] * model_parameters['initial_oxygen']
        # Macrophages
        for ima in initial_macrophage_addresses:
            mac = Macrophage(ima, 'resting')
            initialisation['contents'][ima] = mac
            self.macrophages.append(mac)
        # Fast bacteria
        for ifba in initial_fast_bacteria_addresses:
            fbac = Bacterium(ifba, 'fast')
            initialisation['contents'][ifba] = fbac
            self.bacteria.append(fbac)
        # Fast bacteria
        for isba in initial_slow_bacteria_addresses:
            sbac = Bacterium(isba, 'slow')
            initialisation['contents'][isba] = sbac
            self.bacteria.append(sbac)

        # Set initial diffusion rates (will reduce with caseum)
        for x in range(shape[0]):
            for y in range(shape[1]):
                initialisation['oxygen_diffusion_rate'][(x,y)] = model_parameters['oxygen_diffusion']
                initialisation['chemotherapy_diffusion_rate'][(x, y)] = model_parameters['chemotherapy_diffusion']

        # Hard-coded column headers for recording
        self.values_to_record = ["fast_bacteria", "fast_bacteria_resting", "slow_bacteria", "slow_bacteria_resting",
                                 "intracellular_bac", "total_bacteria",
                                 "resting_macrophages", "active_macrophages", "infected_macrophages",
                                 "chron_infected_macrophages", "total_macrophages",
                                 "t_cells", "caseum"]

        self.grids_to_record = ['oxygen', 'chemotherapy', 'chemokine', 'contents']

        # Add an attribute for the maximum neighbourhood depth
        model_parameters['max_depth'] = 3

        # Super class initialisation
        Automaton.__init__(self, shape, attributes, formats, time_parameters, model_parameters, output_location,
                           self.values_to_record, self.grids_to_record, initialisation)

        # Maxima
        self.max_oxygen = 0.0
        self.max_chemotherapy = 0.0
        self.max_chemokine = 0.0

        # Chemotherapy scheduling
        self.chemo_schedule1_start = np.random.randint(self.model_parameters['chemotherapy_schedule1_start_lower'],
                                                       self.model_parameters['chemotherapy_schedule1_start_upper'])

    # OVERRIDE
    def timestep_output(self):
        """
        Output to console at each timestep
        :return:
        """
        print "t = ", self.time * self.time_step, "Bac = ", len(self.bacteria)

    # OVERRIDE
    def record_counts(self):
        # Count up the totals of each bacteria, macrophage, etc and write the to the file
        writer = csv.writer(self.count_file, delimiter=',')

        fast_bac_count = len([b for b in self.bacteria if b.metabolism == 'fast' and not b.resting])
        fast_bac_rest_count = len([b for b in self.bacteria if b.metabolism == 'fast' and b.resting])
        slow_bac_count = len([b for b in self.bacteria if b.metabolism == 'slow' and not b.resting])
        slow_bac_rest_count = len([b for b in self.bacteria if b.metabolism == 'slow' and b.resting])
        intracell_bac_count = sum([m.intracellular_bacteria for m in self.macrophages])
        total_bac_count = fast_bac_count + fast_bac_rest_count + slow_bac_count + slow_bac_rest_count + \
                          intracell_bac_count
        rest_mac_count = len([m for m in self.macrophages if m.state == 'resting'])
        active_mac_count = len([m for m in self.macrophages if m.state == 'active'])
        inf_mac_count = len([m for m in self.macrophages if m.state == 'infected'])
        chr_inf_mac_count = len([m for m in self.macrophages if m.state == 'chronically_infected'])
        total_mac_count = rest_mac_count + active_mac_count + inf_mac_count + chr_inf_mac_count
        t_cell_count = len(self.t_cells)
        caseum_count = len(self.caseum_addresses)

        row = [self.time * self.time_step, fast_bac_count, fast_bac_rest_count, slow_bac_count, slow_bac_rest_count,
               intracell_bac_count, total_bac_count, rest_mac_count, active_mac_count, inf_mac_count, chr_inf_mac_count,
               total_mac_count, t_cell_count, caseum_count]

        writer.writerow(row)

    # OVERRIDE
    def update_cells(self):
        """
        Run the cellular automaton update. Runs pre-process first to determine diffusion rates. Then runs diffusion
        to calculate new values (written to work grid)
        :return:
        """
        # Update the current maxima
        self.max_oxygen = self.grid['oxygen'].max()
        self.max_chemotherapy = self.grid['chemotherapy'].max()
        self.max_chemokine = self.grid['chemokine'].max()

        self.diffusion_pre_process()
        chemo = (self.chemo_schedule1_start / self.time_step) <= self.time < \
                (self.model_parameters['chemotherapy_schedule1_end'] / self.time_step) or \
                self.model_parameters['chemotherapy_schedule2_start'] / self.time_step <= self.time
        self.diffusion(chemo)

    def diffusion_pre_process(self):
        """
        Pre-processing to calculate diffusion rates. If a cell is too close to too much caseum (determined by model
        parameters caseum_distance_to_reduce_diffusion and caseum_threshold_to_reduce_diffusion, then it's diffusion
        rates (and excretion rate if a blood vessel is present) are decreased.
        :return:
        """

        # Loop through every caseum address and record the addresses that are within the required distance
        affected_addresses = []
        for caseum_address in self.caseum_addresses:
            for d in range(1, self.model_parameters['caseum_distance_to_reduce_diffusion']+1):
                neighbours = self.moore_neighbours(caseum_address,d)
                affected_addresses += neighbours

        # Affected addresses is now a list of all address within the range of cells with caseum. Multiple records of the
        # same address imply cell is too close to more than one caseum. So count how often each address appears.
        counted = Counter(affected_addresses)
        for address in counted:
            # If the count is greater than the threshold, reduce diffusion
            if counted[address] >= self.model_parameters['caseum_threshold_to_reduce_diffusion']:
                self.grid[address]['oxygen_diffusion_rate'] = self.model_parameters['oxygen_diffusion'] / \
                                                              self.model_parameters['diffusion_caseum_reduction']
                self.grid[address]['chemotherapy_diffusion_rate'] = self.model_parameters['chemotherapy_diffusion'] / \
                                                        self.model_parameters['diffusion_caseum_reduction']
                # Reduce excretion if blood vessel
                if self.grid[address]['blood_vessel'] > 0.0:
                    self.grid[address]['blood_vessel'] /= self.model_parameters['diffusion_caseum_reduction']

    def diffusion(self, chemo):
        """
        Calculate new values from diffusion of chemicals (oxygen, chemotherapy, chemokine). Finite difference scheme
        based on the differences between value in the cell and values in von Neumann neighbours to depth 1, and the
        contents of the cell. Different equations for center cells (who have all 4 neighbours) and edge cells (who have
        less than 4).
        :param chemo: Boolean to indicate if chemo is present.
        :return:
        """

        # TODO - veeeeery long (but based off TBModel.cpp). Have to cater for edge cases separately as they have
        # different equations

        # Diffusion works by taking a section of the grid, taking corresponding neighbour grids of same size by slicing
        # the original grid with an offset. Results in 5 grids of the same size. Then can apply the diffusion function
        # to all cells at once within that sub-grid

        # Need to check contents. Function returns 1 if bacterium, 2 if non-resting macrophage and 0 otherwise
        def check_contents(cell):
            if isinstance(cell, Bacterium):
                return 1
            elif (isinstance(cell, Macrophage) and cell.state != 'resting'):
                return 2
            else:
                return 0
        # Run the vectorised function against the main grid
        vfunction_contents = np.vectorize(check_contents)
        contents_check_grid = vfunction_contents(self.grid['contents'])

        # Checks values in grid to get grid of 0/1 to indicate presence of bacterium/non-resting macrophage
        bac_grid = contents_check_grid == 1
        non_resting_mac_grid = contents_check_grid == 2

        # Center grid (of size X-2 x Y-2)
        cell = self.grid[1:-1, 1:-1]
        above = self.grid[:-2, 1:-1]
        below = self.grid[2:, 1:-1]
        left = self.grid[1:-1, :-2]
        right = self.grid[1:-1, 2:]
        cell_has_bacteria = bac_grid[1:-1, 1:-1]
        # Multiply by -1 to give grid of 0s and 1s for macrophages
        cell_has_non_resting_macrophage = non_resting_mac_grid[1:-1, 1:-1]

        # oxygen
        self.work_grid['oxygen'][1:-1,1:-1] = cell['oxygen'] + self.time_step * \
            (((((cell['oxygen_diffusion_rate'] + below['oxygen_diffusion_rate'])/2) *
               (below['oxygen'] - cell['oxygen']) -
               ((cell['oxygen_diffusion_rate'] + above['oxygen_diffusion_rate'])/2) *
                (cell['oxygen'] - above['oxygen'])) /
            self.model_parameters['spatial_step'] ** 2) +
            ((((cell['oxygen_diffusion_rate'] + right['oxygen_diffusion_rate'])/2) *
                (right['oxygen'] - cell['oxygen']) -
              ((cell['oxygen_diffusion_rate'] + left['oxygen_diffusion_rate'])/2) *
                (cell['oxygen'] - left['oxygen'])) /
            self.model_parameters['spatial_step'] ** 2) +
            (self.model_parameters['oxygen_from_source'] * cell['blood_vessel']) -
            (self.model_parameters['oxygen_uptake_from_bacteria'] * cell['oxygen'] * cell_has_bacteria))

        # chemotherapy
        if chemo:
            self.work_grid['chemotherapy'][1:-1, 1:-1] = cell['chemotherapy'] + self.time_step * \
                (((((cell['chemotherapy_diffusion_rate'] + below['chemotherapy_diffusion_rate']) / 2) *
                    (below['chemotherapy'] - cell['chemotherapy']) -
                    ((cell['chemotherapy_diffusion_rate'] + above['chemotherapy_diffusion_rate']) / 2) *
                    (cell['chemotherapy'] - above['chemotherapy'])) /
                    self.model_parameters['spatial_step'] ** 2) +
                    ((((cell['chemotherapy_diffusion_rate'] + right['chemotherapy_diffusion_rate']) / 2) *
                    (right['chemotherapy'] - cell['chemotherapy']) -
                    ((cell['chemotherapy_diffusion_rate'] + left['chemotherapy_diffusion_rate']) / 2) *
                    (cell['chemotherapy'] -left['chemotherapy'])) /
                    self.model_parameters['spatial_step'] ** 2) +
                    (self.model_parameters['chemotherapy_from_source'] * cell['blood_vessel']) -
                    (self.model_parameters['chemotherapy_decay'] * cell['chemotherapy']))

        self.work_grid['chemokine'][1:-1, 1:-1] = cell['chemokine'] + self.time_step * \
            (((self.model_parameters['chemokine_diffusion'] * (below['chemokine'] - cell['chemokine']) -
               self.model_parameters['chemokine_diffusion'] * (cell['chemokine'] - above['chemokine'])) /
              self.model_parameters['spatial_step'] ** 2) +
            ((self.model_parameters['chemokine_diffusion'] * (right['chemokine'] - cell['chemokine']) -
              self.model_parameters['chemokine_diffusion'] * (cell['chemokine'] - left['chemokine'])) /
             self.model_parameters['spatial_step'] ** 2) +
             self.model_parameters['chemokine_from_bacteria'] * cell_has_bacteria +
             (self.model_parameters['chemokine_from_macrophage'] * cell_has_non_resting_macrophage) -
             self.model_parameters['chemokine_decay'] * cell['chemokine'])

        # Edges
        # Top row (no corners) - size (X-2 x Y)
        cell = self.grid[:1, 1:-1]
        below = self.grid[1:2, 1:-1]
        left = self.grid[:1, 0:-2]
        right = self.grid[:1, 2:]
        cell_has_bacteria = bac_grid[:1, 1:-1]
        cell_has_non_resting_macrophage = non_resting_mac_grid[:1, 1:-1] * -1
        work_grid = self.work_grid[:1, 1:-1]
        self.diffusion_3_neighbours(chemo, cell, [left, right], below, cell_has_bacteria,
                                    cell_has_non_resting_macrophage, work_grid)

        # LEFT
        cell = self.grid[1:-1, :1]
        above = self.grid[0:-2, :1]
        right = self.grid[1:-1, 1:2]
        below = self.grid[2:, :1]
        cell_has_bacteria = bac_grid[1:-1, :1]
        cell_has_non_resting_macrophage = non_resting_mac_grid[1:-1, :1] * -1
        work_grid = self.work_grid[1:-1, :1]
        self.diffusion_3_neighbours(chemo, cell, [above, below], right, cell_has_bacteria,
                                    cell_has_non_resting_macrophage, work_grid)

        # RIGHT
        cell = self.grid[1:-1, -1:]
        above = self.grid[0:-2, -1:]
        left = self.grid[1:-1, -2:-1]
        below = self.grid[2:, -1:]
        cell_has_bacteria = bac_grid[1:-1, -1:]
        cell_has_non_resting_macrophage = non_resting_mac_grid[1:-1, -1:] * -1
        work_grid = self.work_grid[1:-1, -1:]
        self.diffusion_3_neighbours(chemo, cell, [above, below], left, cell_has_bacteria,
                                    cell_has_non_resting_macrophage, work_grid)

        # BOTTOM ROW
        cell = self.grid[-1:, 1:-1]
        above = self.grid[-2:-1, 1:-1]
        left = self.grid[-1:, 0:-2]
        right = self.grid[-1:, 2:]
        cell_has_bacteria = bac_grid[-1:, 1:-1]
        cell_has_non_resting_macrophage = non_resting_mac_grid[-1:, 1:-1] * -1
        work_grid = self.work_grid[-1:, 1:-1]
        self.diffusion_3_neighbours(chemo, cell, [left, right], above, cell_has_bacteria,
                                    cell_has_non_resting_macrophage, work_grid)

        #CORNERS
        # Top left
        cell = self.grid[:1, :1]
        below = self.grid[:1, 1:2]
        right = self.grid[1:2, :1]
        cell_has_bacteria = bac_grid[:1, :1]
        cell_has_non_resting_macrophage = non_resting_mac_grid[:1, :1] * -1
        work_grid = self.work_grid[:1, :1]
        self.diffusion_2_neighbours(chemo, cell, [below, right], cell_has_bacteria,
                                    cell_has_non_resting_macrophage, work_grid)

        # Top right
        cell = self.grid[:1, -1:]
        below = self.grid[1:2, -1:]
        left = self.grid[:1, -2:-1]
        cell_has_bacteria = bac_grid[:1, :1]
        cell_has_non_resting_macrophage = non_resting_mac_grid[:1, -1:] * -1
        work_grid = self.work_grid[:1, -1:]
        self.diffusion_2_neighbours(chemo, cell, [below, left], cell_has_bacteria,
                                    cell_has_non_resting_macrophage, work_grid)

        # Bottom left
        cell = self.grid[-1:, :1]
        above = self.grid[-2:-1, :1]
        right = self.grid[-1:, 1:2]
        cell_has_bacteria = bac_grid[-1:, :1]
        cell_has_non_resting_macrophage = non_resting_mac_grid[-1:, :1] * -1
        work_grid = self.work_grid[-1:, :1]
        self.diffusion_2_neighbours(chemo, cell, [above, right], cell_has_bacteria,
                                    cell_has_non_resting_macrophage, work_grid)

        # Bottom right
        cell = self.grid[-1:, -1:]
        above = self.grid[-2:-1, -1:]
        left = self.grid[-1:, -2:-1]
        cell_has_bacteria = bac_grid[-1:, -1:]
        cell_has_non_resting_macrophage = non_resting_mac_grid[-1:, -1:] * -1
        work_grid = self.work_grid[-1:, -1:]
        self.diffusion_2_neighbours(chemo, cell, [above, left], cell_has_bacteria,
                                    cell_has_non_resting_macrophage, work_grid)

        if not chemo:
            self.work_grid['chemotherapy'] = np.zeros(self.grid.shape,dtype=float)

    def diffusion_3_neighbours(self, chemo, cell, paired_neighbours, non_paired_neighbour, cell_has_bacteria,
                               cell_has_non_resting_macrophage, work_grid):
        """
        For cells with 3 neighbours (4th would be off grid - i.e. top row, left column etc.) run diffusion rules
        :param chemo: 
        :param cell: 
        :param paired_neighbours: 
        :param non_paired_neighbour: 
        :param cell_has_bacteria: 
        :param cell_has_non_resting_macrophage: 
        :param work_grid: 
        :return: 
        """
        work_grid['oxygen'] = cell['oxygen'] + self.time_step * (
            ((cell['oxygen_diffusion_rate'] * (non_paired_neighbour['oxygen'] - 2 * cell['oxygen'] + non_paired_neighbour['oxygen'])) /
             self.model_parameters['spatial_step'] ** 2) +
            ((cell['oxygen_diffusion_rate'] * (paired_neighbours[0]['oxygen'] - 2 * cell['oxygen'] + paired_neighbours[1]['oxygen'])) /
             self.model_parameters['spatial_step'] ** 2) +
            (self.model_parameters['oxygen_from_source'] * cell['blood_vessel']) +
            (self.model_parameters['oxygen_uptake_from_bacteria'] * cell['oxygen'] *
             cell_has_bacteria)
        )
        if chemo:
            work_grid['chemotherapy'] = cell['chemotherapy'] + self.time_step * (
                ((cell['chemotherapy_diffusion_rate'] * (non_paired_neighbour['chemotherapy'] - 2 * cell['chemotherapy'] +
                                                         non_paired_neighbour['chemotherapy'])) /
                 self.model_parameters['spatial_step'] ** 2) +
                ((cell['chemotherapy_diffusion_rate'] * (paired_neighbours[0]['chemotherapy'] - 2 * cell['chemotherapy'] +
                                                         paired_neighbours[1]['chemotherapy'])) /
                 self.model_parameters['spatial_step'] ** 2) +
                (self.model_parameters['chemotherapy_from_source'] * cell['blood_vessel']) +
                (self.model_parameters['chemotherapy_decay'] * cell['chemotherapy'])
            )
        work_grid['chemokine'] = cell['chemokine'] + self.time_step * (
            ((self.model_parameters['chemokine_diffusion'] * (non_paired_neighbour['chemokine'] - 2 * cell['chemokine'] +
                                                              non_paired_neighbour['chemokine'])) /
             self.model_parameters['spatial_step'] ** 2) +
            ((self.model_parameters['chemokine_diffusion'] * (paired_neighbours[0]['chemokine'] - 2 * cell['chemokine'] +
                                                              paired_neighbours[1]['chemokine'])) /
             self.model_parameters['spatial_step'] ** 2) +
            self.model_parameters['chemokine_from_bacteria'] * cell_has_bacteria +
            self.model_parameters['chemokine_from_macrophage'] *
            (cell_has_non_resting_macrophage) +
            self.model_parameters['chemokine_decay'] * cell['chemokine']
        )

    def diffusion_2_neighbours(self, chemo, cell, neighbours, cell_has_bacteria,
                               cell_has_non_resting_macrophage, work_grid):
        """
        For cells with 2 neighbours (other 2 would be off grid - i.e. corners) run diffusion rules
        :param chemo:
        :param cell:
        :param neighbours:
        :param cell_has_bacteria:
        :param cell_has_non_resting_macrophage:
        :param work_grid:
        :return:
        """

        work_grid['oxygen'] = cell['oxygen'] + self.time_step * (
            ((cell['oxygen_diffusion_rate'] * (neighbours[0]['oxygen'] - 2 * cell['oxygen'] + neighbours[0]['oxygen'])) /
             self.model_parameters['spatial_step'] ** 2) +
            ((cell['oxygen_diffusion_rate'] * (neighbours[1]['oxygen'] - 2 * cell['oxygen'] + neighbours[1]['oxygen'])) /
             self.model_parameters['spatial_step'] ** 2) +
            (self.model_parameters['oxygen_from_source'] * cell['blood_vessel']) +
            (self.model_parameters['oxygen_uptake_from_bacteria'] * cell['oxygen'] *
             cell_has_bacteria)
        )
        if chemo:
            work_grid['chemotherapy'] = cell['chemotherapy'] + self.time_step * (
                ((cell['chemotherapy_diffusion_rate'] * (neighbours[0]['chemotherapy'] - 2 * cell['chemotherapy'] +
                                                         neighbours[0]['chemotherapy'])) /
                 self.model_parameters['spatial_step'] ** 2) +
                ((cell['chemotherapy_diffusion_rate'] * (neighbours[1]['chemotherapy'] - 2 * cell['chemotherapy'] +
                                                         neighbours[1]['chemotherapy'])) /
                 self.model_parameters['spatial_step'] ** 2) +
                (self.model_parameters['chemotherapy_from_source'] * cell['blood_vessel']) +
                (self.model_parameters['chemotherapy_decay'] * cell['chemotherapy'])
            )
        work_grid['chemokine'] = cell['chemokine'] + self.time_step * (
            ((self.model_parameters['chemokine_diffusion'] * (neighbours[0]['chemokine'] - 2 * cell['chemokine'] +
                                                              neighbours[0]['chemokine'])) /
             self.model_parameters['spatial_step'] ** 2) +
            ((self.model_parameters['chemokine_diffusion'] * (neighbours[1]['chemokine'] - 2 * cell['chemokine'] +
                                                              neighbours[1]['chemokine'])) /
             self.model_parameters['spatial_step'] ** 2) +
            self.model_parameters['chemokine_from_bacteria'] * cell_has_bacteria +
            self.model_parameters['chemokine_from_macrophage'] *
            (cell_has_non_resting_macrophage) +
            self.model_parameters['chemokine_decay'] * cell['chemokine']
        )

    # OVERRIDE
    def generate_events_from_agents(self):
        events = []
        events += self.bacteria_replication()
        events += self.t_cell_recruitment()
        events += self.macrophage_recruitment()
        events += self.chemotherapy_killing_bacteria()
        events += self.chemotherapy_killing_macrophages()
        events += self.t_cell_processes()
        events += self.macrophage_processes()
        events += self.macrophage_state_changes()
        events += self.bacteria_state_changes()
        return events

    def oxygen_scale(self, address):
        if self.max_oxygen == 0.0:
            return 0.0
        else:
            return self.grid[address]['oxygen'] / self.max_oxygen

    def chemotherapy_scale(self, address):
        if self.max_chemotherapy == 0.0:
            return 0.0
        else:
            return self.grid[address]['chemotherapy'] / self.max_chemotherapy

    def chemokine_scale(self, address):
        if self.max_chemokine == 0.0:
            return 0.0
        else:
            return self.grid[address]['chemokine'] / self.max_chemokine

    def total_bacteria(self):
        return len(self.bacteria) + sum([m.intracellular_bacteria for m in self.macrophages])

    def find_max_chemokine_neighbour(self, neighbours):
        """
        Given neighbour addresses, find the neighbour which has the highest level of chemokine
        :param neighbours:
        :return:
        """
        max_chemokine = 0.0
        highest_indices = []
        for index in neighbours:

            if neighbours[index]['chemokine'] > max_chemokine:
                max_chemokine = neighbours[index]['chemokine']
                highest_indices = [index]
            elif neighbours[index]['chemokine'] == max_chemokine:
                highest_indices.append(index)

        # Tie-breaking. If just one pick it, else pick any one index at random
        choice = np.random.randint(0, len(highest_indices))
        chosen_index = highest_indices[choice]

        return [chosen_index, self.chemokine_scale(chosen_index)]

    def bacteria_replication(self):
        """
        Bacteria replicate (produce a new bacterium agent) once they reach a certain age.
        :return:
        """
        replication_events = []
        # Loop through every bacteria, check age against a (stochastic) threshold, generate event if age is higher than
        # threshold
        for bacterium in self.bacteria:
            # Increment age
            bacterium.age += self.time_step

            # Skip if the bacterium is resting
            if bacterium.resting:
                continue

            if bacterium.metabolism == 'fast':
                maximum = self.model_parameters['bacteria_replication_fast_upper']
                minimum = self.model_parameters['bacteria_replication_fast_lower']
            else:  # Slow
                maximum = self.model_parameters['bacteria_replication_slow_upper']
                minimum = self.model_parameters['bacteria_replication_slow_lower']

            replication_time = np.random.randint(minimum, maximum) / self.time_step

            # TODO - MED - Does this really work as a modulo?
            # If the time is sufficient enough, bacteria can replicate
            if self.time % replication_time == 0:

                # Look for free neighbours
                free_neighbours = []
                # TODO - COMP - maybe 4 shouldn't be hard coded?
                for depth in range(1, 4):
                    # Pull the neighbours from the appropriate neighbourhood
                    if bacterium.division_neighbourhood == 'mo':
                        neighbours = self.moore_neighbours(bacterium.address, depth)
                    else:
                        neighbours = self.von_neumann_neighbours(bacterium.address, depth)
                    # Find a free neighbour (not a blood vessel and contents == 0.0)
                    for neighbour_address in neighbours:
                        neighbour = self.grid[neighbour_address]
                        if neighbour is not None and neighbour['contents'] == 0.0 and neighbour['blood_vessel'] == 0.0:
                            free_neighbours.append(neighbour_address)
                    # If a free neighbour found, don't look at greater depths
                    if len(free_neighbours) > 0:
                        break
                # A free neighbour has not been found anywhere
                if len(free_neighbours) == 0:
                    # Bacterium will change to resting state (quorum sensing)
                    new_event = BacteriumStateChange(bacterium.address, 'resting', True)
                    replication_events.append(new_event)
                else:  # Free space found
                    # Pick a free neighbour at random
                    neighbour_address = free_neighbours[np.random.randint(len(free_neighbours))]
                    # Create event and add to list of potential events
                    new_event = BacteriumReplication(bacterium.address, neighbour_address, bacterium.metabolism)
                    replication_events.append(new_event)
        return replication_events

    def t_cell_recruitment(self):
        """
        Once bacteria over entire system reach a threshold, t-cells enter the system. Creates an event to add a t-cell
        to a cell next to a blood vessel
        :return:
        """
        t_cell_recruitment_events = []
        # When global amount of bacteria exceeds threshold
        if self.total_bacteria() >= self.model_parameters['bacteria_threshold_for_t_cells']:
            # Each blood vessel
            for blood_vessel_address in self.blood_vessel_addresses:
                # Generate event if probability according to parameters
                r = np.random.randint(1, 101)
                if r <= self.model_parameters['t_cell_recruitment_probability']:
                    # Get von Neumann neighbours
                    neighbours = self.von_neumann_neighbours(blood_vessel_address,1)
                    # Loop through all neighbours to find suitable options
                    free_neighbours = []
                    for neighbour_address in neighbours:
                        neighbour = self.grid[neighbour_address]
                        # Check neighbour is on the grid, is empty and has a sufficiently high chemokine level
                        if neighbour is not None and neighbour['blood_vessel'] == 0.0 and neighbour['contents'] == 0.0 \
                                and self.chemokine_scale(neighbour_address) > \
                                self.model_parameters['chemokine_scale_for_t_cell_recruitment']:
                            free_neighbours.append(neighbour_address)
                    # Check there is at least one suitable neighbour
                    if len(free_neighbours) > 0:
                        # Pick one of the neighbours
                        neighbour_address = free_neighbours[np.random.randint(len(free_neighbours))]
                        # Create event
                        new_event = RecruitTCell(blood_vessel_address, neighbour_address)
                        t_cell_recruitment_events.append(new_event)
        return t_cell_recruitment_events

    def macrophage_recruitment(self):
        """
        Each step for each source vessel, there is a probability that macrophage will be recruited
        :return:
        """
        recruitment_events = []
        if self.total_bacteria() >= self.model_parameters['bacteria_threshold_for_macrophage_recruitment']:
            chemokine_threshold = self.model_parameters['chemokine_scale_for_macrophage_recruitment_above_threshold']
        else:
            chemokine_threshold = self.model_parameters['chemokine_scale_for_macrophage_recruitment_below_threshold']

        # Loop through each blood vessel
        for bv_address in self.blood_vessel_addresses:
            # Generate event with probability based on parameters
            r = np.random.randint(1, 101)
            if r <= self.model_parameters['macrophage_recruitment_probability']:
                # Get neighbours, then reduce to those that are free and have sufficient chemokine scale
                neighbours = self.von_neumann_neighbours(bv_address, 1)
                free_neighbours = []
                for neighbour_address in neighbours:
                    neighbour = self.grid[neighbour_address]
                    if neighbour is not None and neighbour['blood_vessel'] == 0.0 and neighbour['contents'] == 0.0 and \
                            self.chemokine_scale(neighbour_address) > chemokine_threshold:
                        free_neighbours.append(neighbour_address)

                if len(free_neighbours) > 0:
                    # Pick one of the neighbours
                    chosen_neighbour = free_neighbours[np.random.randint(len(free_neighbours))]
                    # Create event
                    new_event = RecruitMacrophage(bv_address, chosen_neighbour)
                    recruitment_events.append(new_event)
        return recruitment_events

    def chemotherapy_killing_bacteria(self):
        """
        Chemotherapy destroys bacterium if the level is high enough
        :return:
        """
        chemo_kill_bac_events = []
        # Loop through all bacteria
        for bacterium in self.bacteria:
            # Check chemotherapy scale against relevant parameter based on metabolism
            chemo_scale = self.chemotherapy_scale(bacterium.address)
            if (bacterium.metabolism == 'fast' and chemo_scale >
                    self.model_parameters['chemotherapy_scale_for_kill_fast_bacteria']) \
                    or (bacterium.metabolism == 'slow' and chemo_scale >
                    self.model_parameters['chemotherapy_scale_for_kill_slow_bacteria']):
                # Scale is high enough, so create event to destroy bacterium
                new_event = ChemoKillBacterium(bacterium.address)
                chemo_kill_bac_events.append(new_event)
        return chemo_kill_bac_events

    def chemotherapy_killing_macrophages(self):
        chemo_kill_mac_events = []
        # Loop through all macrophages
        for macrophage in self.macrophages:
            # Check chemotherapy scale against relevant parameter based on metabolism
            chemo_scale = self.chemotherapy_scale(macrophage.address)
            if (macrophage.state == 'infected' or macrophage.state == 'chronically_infected') \
                and chemo_scale > self.model_parameters['chemotherapy_scale_for_kill_macrophage']:
                # Scale is high enough, so create event to destroy bacterium
                new_event = ChemoKillMacrophage(macrophage.address)
                chemo_kill_mac_events.append(new_event)
        return chemo_kill_mac_events

    def t_cell_processes(self):
        """
        T-cells movement, death and apoptosis of other agents
        :return:
        """
        t_cell_events = []

        # T-cells only move after set period of time
        if self.time % self.model_parameters['t_cell_movement_time'] == 0:

            # Loop through all T-cells
            for t_cell in self.t_cells:
                # Increment age
                t_cell.age += self.time_step
                # Stochastic age threshold
                age_threshold = np.random.randint(0, self.model_parameters['t_cell_age_threshold'])
                # T-CELL DEATH
                # If age > threshold, t-cell dies
                if t_cell.age >= age_threshold:
                    new_event = TCellDeath(t_cell.address)
                    t_cell_events.append(new_event)
                else:  # T-CELL MOVE
                    # T-cells move in biased random walk. Determine if move will be random based on probability in
                    # parameters
                    random_move = False
                    prob_random_move = np.random.randint(1, 101)
                    if prob_random_move <= self.model_parameters['t_cell_random_move_probability']:
                        random_move = True
                    # Get neighbours
                    neighbours = self.moore_neighbours(t_cell.address, 1)
                    # If a random move, pick a neighbour at random
                    if random_move:
                        index = np.random.randint(0, len(neighbours))
                        chosen_neighbour_address = neighbours.keys()[index]
                    else: # Pick the neighbour with the highest chemokine level
                        chosen_neighbour_address = self.find_max_chemokine_neighbour(neighbours)[0]

                    # Get neighbour
                    neighbour = self.grid[chosen_neighbour_address]
                    # Check neighbour is empty, then move T-cell there
                    if neighbour['contents'] == 0.0 and neighbour['blood_vessel'] == 0.0:
                        new_event = TCellMovement(t_cell.address, chosen_neighbour_address)
                        t_cell_events.append(new_event)
                    # Else if the address contains an infected macrophage, then t-cell may kill it
                    elif isinstance(neighbour['contents'], Macrophage) and (neighbour['contents'].state == 'infected'
                            or neighbour['contents'].state == 'chronically_infected'):
                        # T-cell killing based on parameter probability
                        prob_t_cell_killing = np.random.randint(1, 101)
                        if prob_t_cell_killing <= self.model_parameters['t_cell_kills_macrophage_probability']:
                            new_event = TCellKillsMacrophage(t_cell.address, chosen_neighbour_address)
                            t_cell_events.append(new_event)
        return t_cell_events

    def macrophage_processes(self):
        """
        Macrophages move, die and ingest bacteria
        :return:
        """
        mac_events = []
        # Loop through macrophages
        for macrophage in self.macrophages:
            death = False
            move = False
            ingest = False
            # Increment age
            macrophage.age += self.time_step
            # Different events/movement rates/death rates depending on state
            # TODO - MED - time > 1/dt added to match TBModel.cpp - but what is significance of this?
            if macrophage.state == 'resting' and self.time > 1 / self.time_step:
                # Death by age is stochastic
                random_macrophage_age = np.random.randint(0, self.model_parameters['resting_macrophage_age_limit'])
                if macrophage.age >= random_macrophage_age:
                    death = True
                # Within a set time for movement
                if (not death) and self.time % self.model_parameters['resting_macrophage_movement_time'] == 0:
                    # Chemokine moves on random biased walk. Random move with probability based on parameters, if
                    # highest chemokine scale at neighbours does not exceed threshold, then also random move
                    neighbours = self.moore_neighbours(macrophage.address, 1)
                    max_chemokine_address, max_chemokine_scale = self.find_max_chemokine_neighbour(neighbours)
                    # Generate random number for probability of random move
                    prob_random_move = np.random.randint(1, 101)
                    random_move = False
                    if prob_random_move <= self.model_parameters['prob_resting_macrophage_random_move'] \
                            or max_chemokine_scale <= \
                            self.model_parameters['minimum_chemokine_for_resting_macrophage_movement']:
                        random_move = True
                    # Pick the neighbour to move to, either random or highest chemokine scale
                    if random_move:
                        chosen_neighbour_address = neighbours.keys()[np.random.randint(0, len(neighbours))]
                    else:
                        chosen_neighbour_address = max_chemokine_address
                    # Check if leaving the grid
                    neighbour = self.grid[chosen_neighbour_address]
                    # If neighbour is empty, create a move event
                    if neighbour['contents'] == 0.0 and neighbour['blood_vessel'] == 0.0:
                        move = True
                    # If neighbour contains a bacterium, ingest it
                    elif isinstance(neighbour['contents'], Bacterium):
                        ingest = True
            # Active macrophage processes
            elif macrophage.state == 'active':
                # Active macrophages die after a set time (not stochastic)
                if macrophage.age > self.model_parameters['active_macrophage_age_limit']:
                    death = True
                # Set time for macrophage movement
                if (not death) and self.time % self.model_parameters['active_macrophage_movement_time'] == 0:
                    # Active macrophages always move to highest chemokine neighbour
                    neighbours = self.moore_neighbours(macrophage.address, 1)
                    chosen_neighbour_address = self.find_max_chemokine_neighbour(neighbours)[0]
                    neighbour = self.grid[chosen_neighbour_address]
                    # If cell to move to has a bacterium
                    if isinstance(neighbour['contents'], Bacterium):
                        # Macrophages ingests with set probability (active macrophages will destroy)
                        prob_macrophage_ingest = np.random.randint(1, 101)
                        # Probabilities differ based on bacterium metabolism
                        if (neighbour['contents'].metabolism == 'fast' and prob_macrophage_ingest <=
                            self.model_parameters['prob_active_macrophage_kill_fast_bacteria']) or (
                                neighbour['contents'].metabolism == 'slow' and prob_macrophage_ingest <=
                                self.model_parameters['prob_active_macrophage_kill_slow_bacteria']):
                            ingest = True
                    # Cell is empty so create a move event
                    elif neighbour['contents'] == 0.0 and neighbour['blood_vessel'] == 0.0:
                        move = True
            # Infected Macrophage processes
            elif macrophage.state == 'infected':
                # Death is stochastic
                random_macrophage_age = np.random.randint(0, self.model_parameters['infected_macrophage_age_limit'])
                if macrophage.age >= random_macrophage_age:
                    death = True
                # Move after certain time
                if (not death) and self.time % self.model_parameters['infected_macrophage_movement_time'] == 0:
                    # Infected move to highest chemokine neighbour
                    neighbours = self.moore_neighbours(macrophage.address, 1)
                    chosen_neighbour_address = self.find_max_chemokine_neighbour(neighbours)[0]
                    neighbour = self.grid[chosen_neighbour_address]
                    # Neighbour is empty, so move event
                    if neighbour['contents'] == 0.0 and neighbour['blood_vessel'] == 0.0:
                        move = True
                    # Neighbour has a bacterium, so kill event
                    elif isinstance(neighbour['contents'], Bacterium):
                        ingest = True
            # Chronically infected macrophage processes
            elif macrophage.state == 'chronically_infected':
                # Stochastic death
                random_macrophage_age = np.random.randint(0,
                                                    self.model_parameters['chronically_infected_macrophage_age_limit'])
                if macrophage.age >= random_macrophage_age:
                    death = True
                # Movement at set times
                if (not death) and self.time % \
                        self.model_parameters['chronically_infected_macrophage_movement_time'] == 0:
                    # Move to highest chemokine scale neighbour
                    neighbours = self.moore_neighbours(macrophage.address, 1)
                    chosen_neighbour_address = self.find_max_chemokine_neighbour(neighbours)[0]
                    neighbour = self.grid[chosen_neighbour_address]
                    # Neighbour is empty, so move event
                    if neighbour['contents'] == 0.0 and neighbour['blood_vessel'] == 0.0:
                        move = True
                    # Neighbour has bacterium, so kill event
                    elif isinstance(neighbour['contents'], Bacterium):
                        ingest = True

            if death:
                new_event = MacrophageDeath(macrophage.address)
                mac_events.append(new_event)
            elif move:
                new_event = MacrophageMovement(macrophage.address, chosen_neighbour_address)
                mac_events.append(new_event)
            elif ingest:
                new_event = MacrophageIngestsBacterium(macrophage.address, chosen_neighbour_address)
                mac_events.append(new_event)

        return mac_events

    def macrophage_state_changes(self):
        """
        Macrophages change state based on chemokine levels and intracellular bacteria
        :return:
        """
        mac_state_change_events = []

        # Loop through macrophages
        for macrophage in self.macrophages:
            new_event = None

            # Resting macrophages
            if macrophage.state == 'resting':
                # RESTING -> ACTIVE
                if self.chemokine_scale(macrophage.address) > \
                        self.model_parameters['chemokine_scale_for_macrophage_activation'] \
                        and macrophage.intracellular_bacteria == 0:
                    new_event = MacrophageChangesState(macrophage.address, 'active')
                # RESTING -> INFECTED
                elif macrophage.intracellular_bacteria == 1:
                    new_event = MacrophageChangesState(macrophage.address, 'infected')
            elif macrophage.state == 'active':
                # ACTIVE -> RESTING
                if self.chemokine_scale(macrophage.address) < \
                        self.model_parameters['chemokine_scale_for_macrophage_deactivation']:
                    new_event = MacrophageChangesState(macrophage.address, 'resting')
            elif macrophage.state == 'infected':
                # INFECTED -> CHRONICALLY INFECTED
                if macrophage.intracellular_bacteria >= self.model_parameters['bacteria_to_turn_chronically_infected']:
                    new_event = MacrophageChangesState(macrophage.address, 'chronically_infected')
            elif macrophage.state == 'chronically_infected':
                # MACROPHAGE BURSTS
                if macrophage.intracellular_bacteria == self.model_parameters['bacteria_to_burst_macrophage']:
                    # Loop through all neighbours (up to depth 3) and try to find enough to distribute bacteria on to
                    bacteria_addresses = []
                    for depth in range(1, 4):
                        neighbours = self.moore_neighbours(macrophage.address, depth).keys()
                        # Shuffle the neighbours so we don't give priority
                        np.random.shuffle(neighbours)
                        for n in neighbours:
                            # Find empty neighbours
                            neighbour = self.grid[n]
                            if neighbour['contents'] == 0.0 and neighbour['blood_vessel'] == 0.0:
                                bacteria_addresses.append(n)
                            # Limit reached - break here stops checking other neighbours at this depth, also need to
                            # stop searching further depths
                            if len(bacteria_addresses) == self.model_parameters['bacteria_to_burst_macrophage']:
                                # Break neighbour loop
                                break
                        # Limit reached earlier so don't check other depths
                        if len(bacteria_addresses) == self.model_parameters['bacteria_to_burst_macrophage']:
                            # Break depth loop
                            break
                    # Macrophage bursting event
                    new_event = MacrophageBursts(macrophage.address, bacteria_addresses)

            if new_event is not None:
                mac_state_change_events.append(new_event)

        return mac_state_change_events

    def bacteria_state_changes(self):
        """
        Bacteria switch between metabolism based on oxygen and switch resting true to false based on space
        :return:
        """
        bac_state_events = []
        # Loop through bacteria
        for bacterium in self.bacteria:
            # Metabolism change only happens later in process (after 2 hours)
            if self.time > 2 / self.time_step:
                # Check if state change - different scales based on metabolism
                if bacterium.metabolism == 'fast' and self.oxygen_scale(bacterium.address) <= \
                        self.model_parameters['oxygen_scale_for_metabolism_change_to_slow']:
                    new_event = BacteriumStateChange(bacterium.address, 'metabolism', 'slow')
                    bac_state_events.append(new_event)
                elif bacterium.metabolism == 'slow' and self.oxygen_scale(bacterium.address) > \
                        self.model_parameters['oxygen_scale_for_metabolism_change_to_fast']:
                    new_event = BacteriumStateChange(bacterium.address, 'metabolism', 'fast')
                    bac_state_events.append(new_event)
            # If bacteria is resting, check if there is now space in neighbourhood, if so, revert to non-resting
            if bacterium.resting:
                space_found = False
                for depth in range(1, 4):
                    # Get neighbours
                    neighbours = self.moore_neighbours(bacterium.address, depth)
                    for n in neighbours:
                        # Is neighbour empty?
                        neighbour = self.grid[n]
                        if neighbour is not None and neighbour['blood_vessel'] == 0.0 and neighbour['contents'] == 0.0:
                            new_event = BacteriumStateChange(bacterium.address, 'resting', False)
                            bac_state_events.append(new_event)
                            space_found = True
                            # Don't check other neighbours
                            break
                    # Space found so don't check further depths
                    if space_found:
                        break
        return bac_state_events
