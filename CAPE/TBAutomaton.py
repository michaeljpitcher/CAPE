from CAPEAutomaton import *
from collections import Counter
# ---------------------------------------- TB -------------------------------------------------------


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
        self.tcells = []
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
            mac = Macrophage('resting')
            initialisation['contents'][ima] = mac
            self.macrophages.append(mac)
        # Fast bacteria
        for ifba in initial_fast_bacteria_addresses:
            fbac = Bacterium('fast')
            initialisation['contents'][ifba] = fbac
            self.bacteria.append(fbac)
        # Fast bacteria
        for isba in initial_slow_bacteria_addresses:
            sbac = Bacterium('slow')
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

        # Chemotherapy scheduling
        self.chemo_schedule1_start = np.random.randint(self.model_parameters['chemotherapy_schedule1_start_lower'],
                                                       self.model_parameters['chemotherapy_schedule1_start_upper'])

    def timestep_output(self):
        """
        Output to console at each timestep
        :return:
        """
        print "t = ", self.time * self.time_step, "Bac = ", len(self.bacteria)

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
        chr_inf_mac_count = len([m for m in self.macrophages if m.state == 'chr_inf'])
        total_mac_count = rest_mac_count + active_mac_count + inf_mac_count + chr_inf_mac_count
        t_cell_count = len(self.tcells)
        caseum_count = len(self.caseum_addresses)

        row = [self.time * self.time_step, fast_bac_count, fast_bac_rest_count, slow_bac_count, slow_bac_rest_count,
               intracell_bac_count, total_bac_count, rest_mac_count, active_mac_count, inf_mac_count, chr_inf_mac_count,
               total_mac_count, t_cell_count, caseum_count]

        writer.writerow(row)

    def update_cells(self):
        """
        Run the cellular automaton update. Runs pre-process first to determine diffusion rates. Then runs diffusion
        to calculate new values (written to work grid)
        :return:
        """
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
            for d in range(self.model_parameters['caseum_distance_to_reduce_diffusion']):
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

        # Center grid (of size X-2 x Y-2)
        cell = self.grid[1:-1, 1:-1]
        above = self.grid[:-2, 1:-1]
        below = self.grid[2:, 1:-1]
        left = self.grid[1:-1, :-2]
        right = self.grid[1:-1, 2:]

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
            (self.model_parameters['oxygen_from_source'] * cell['blood_vessel']) +
            (self.model_parameters['oxygen_uptake_from_bacteria'] * cell['oxygen'] *
             isinstance(cell['contents'], Bacterium)))

        # chemotherapy
        if chemo:
            self.work_grid['chemotherapy'][1:-1, 1:-1] = cell['chemotherapy'] + self.time_step * \
                (((((cell['chemotherapy_diffusion_rate'] + below['chemotherapy_diffusion_rate']) / 2) *
                    (below['chemotherapy'] - cell['chemotherapy']) -
                    ((cell['chemotherapy_diffusion_rate'] +above['chemotherapy_diffusion_rate']) / 2) *
                    (cell['chemotherapy'] - above['chemotherapy'])) /
                    self.model_parameters['spatial_step'] ** 2) +
                    ((((cell['chemotherapy_diffusion_rate'] + right['chemotherapy_diffusion_rate']) / 2) *
                    (right['chemotherapy'] - cell['chemotherapy']) -
                    ((cell['chemotherapy_diffusion_rate'] + left['chemotherapy_diffusion_rate']) / 2) *
                    (cell['chemotherapy'] -left['chemotherapy'])) /
                    self.model_parameters['spatial_step'] ** 2) +
                    (self.model_parameters['chemotherapy_from_source'] * cell['blood_vessel']) +
                    (self.model_parameters['chemotherapy_decay'] * cell['chemotherapy']))

        self.work_grid['chemokine'][1:-1, 1:-1] = cell['chemokine'] + self.time_step * \
            (((self.model_parameters['chemokine_diffusion'] * (below['chemokine'] - cell['chemokine']) -
               self.model_parameters['chemokine_diffusion'] * (cell['chemokine'] - above['chemokine'])) /
              self.model_parameters['spatial_step'] ** 2) +
            ((self.model_parameters['chemokine_diffusion'] * (right['chemokine'] - cell['chemokine']) -
              self.model_parameters['chemokine_diffusion'] * (cell['chemokine'] - left['chemokine'])) /
             self.model_parameters['spatial_step'] ** 2) +
             self.model_parameters['chemokine_from_bacteria'] * isinstance(cell['contents'], Bacterium) +
             self.model_parameters['chemokine_from_macrophage'] *
             (isinstance(cell['contents'], Macrophage) and cell['contents'].state != 'resting') +
             self.model_parameters['chemokine_decay'] * cell['chemokine'])

        # Edges
        # Top row (no corners) - size (X-2 x Y)
        cell = self.grid[:1, 1:-1]
        below = self.grid[1:2, 1:-1]
        left = self.grid[:1, 0:-2]
        right = self.grid[:1, 2:]
        self.work_grid['oxygen'][:1, 1:-1] = cell['oxygen'] + self.time_step * (
            ((cell['oxygen_diffusion_rate'] * (below['oxygen'] - 2 * cell['oxygen'] + below['oxygen'])) /
            self.model_parameters['spatial_step'] ** 2) +
            ((cell['oxygen_diffusion_rate'] * (left['oxygen'] - 2 * cell['oxygen'] + right['oxygen'])) /
            self.model_parameters['spatial_step'] ** 2)+
            (self.model_parameters['oxygen_from_source'] * cell['blood_vessel']) +
            (self.model_parameters['oxygen_uptake_from_bacteria'] * cell['oxygen'] *
             isinstance(cell['contents'], Bacterium))
        )
        if chemo:
            self.work_grid['chemotherapy'][:1, 1:-1] = cell['chemotherapy'] + self.time_step * (
                ((cell['chemotherapy_diffusion_rate'] * (below['chemotherapy'] - 2 * cell['chemotherapy'] +
                below['chemotherapy'])) /
                self.model_parameters['spatial_step'] ** 2) +
                ((cell['chemotherapy_diffusion_rate'] * (left['chemotherapy'] - 2 * cell['chemotherapy'] +
                right['chemotherapy'])) /
                self.model_parameters['spatial_step'] ** 2) +
                (self.model_parameters['chemotherapy_from_source'] * cell['blood_vessel']) +
                (self.model_parameters['chemotherapy_decay'] * cell['chemotherapy'])
        )
        self.work_grid['chemokine'][:1, 1:-1] = cell['chemokine'] + self.time_step * (
            ((self.model_parameters['chemokine_diffusion'] * (below['chemokine'] - 2 * cell['chemokine'] +
              below['chemokine'])) /
             self.model_parameters['spatial_step'] ** 2) +
            ((self.model_parameters['chemokine_diffusion'] * (left['chemokine'] - 2 * cell['chemokine'] +
              right['chemokine'])) /
             self.model_parameters['spatial_step'] ** 2) +
            self.model_parameters['chemokine_from_bacteria'] * isinstance(cell['contents'], Bacterium) +
            self.model_parameters['chemokine_from_macrophage'] *
            (isinstance(cell['contents'], Macrophage) and cell['contents'].state != 'resting') +
            self.model_parameters['chemokine_decay'] * cell['chemokine']
        )

        # LEFT
        cell = self.grid[1:-1, :1]
        above = self.grid[0:-2, :1]
        right = self.grid[1:-1, 1:2]
        below = self.grid[2:, :1]
        self.work_grid['oxygen'][1:-1, :1] = cell['oxygen'] + self.time_step * (
            ((cell['oxygen_diffusion_rate'] * (right['oxygen'] - 2 * cell['oxygen'] + right['oxygen'])) /
             self.model_parameters['spatial_step'] ** 2) +
            ((cell['oxygen_diffusion_rate'] * (above['oxygen'] - 2 * cell['oxygen'] + below['oxygen'])) /
             self.model_parameters['spatial_step'] ** 2) +
            (self.model_parameters['oxygen_from_source'] * cell['blood_vessel']) +
            (self.model_parameters['oxygen_uptake_from_bacteria'] * cell['oxygen'] *
             isinstance(cell['contents'], Bacterium))
        )
        if chemo:
            self.work_grid['chemotherapy'][1:-1, :1] = cell['chemotherapy'] + self.time_step * (
                ((cell['chemotherapy_diffusion_rate'] * (right['chemotherapy'] - 2 * cell['chemotherapy'] +
                                                         right['chemotherapy'])) /
                 self.model_parameters['spatial_step'] ** 2) +
                ((cell['chemotherapy_diffusion_rate'] * (above['chemotherapy'] - 2 * cell['chemotherapy'] +
                                                         above['chemotherapy'])) /
                 self.model_parameters['spatial_step'] ** 2) +
                (self.model_parameters['chemotherapy_from_source'] * cell['blood_vessel']) +
                (self.model_parameters['chemotherapy_decay'] * cell['chemotherapy'])
            )
        self.work_grid['chemokine'][1:-1, :1] = cell['chemokine'] + self.time_step * (
            ((self.model_parameters['chemokine_diffusion'] * (right['chemokine'] - 2 * cell['chemokine'] +
                                                              right['chemokine'])) /
             self.model_parameters['spatial_step'] ** 2) +
            ((self.model_parameters['chemokine_diffusion'] * (above['chemokine'] - 2 * cell['chemokine'] +
                                                              below['chemokine'])) /
             self.model_parameters['spatial_step'] ** 2) +
            self.model_parameters['chemokine_from_bacteria'] * isinstance(cell['contents'], Bacterium) +
            self.model_parameters['chemokine_from_macrophage'] *
            (isinstance(cell['contents'], Macrophage) and cell['contents'].state != 'resting') +
            self.model_parameters['chemokine_decay'] * cell['chemokine']
        )
        # RIGHT
        cell = self.grid[1:-1, -1:]
        above = self.grid[0:-2, -1:]
        left = self.grid[1:-1, -2:-1]
        below = self.grid[2:, -1:]
        self.work_grid['oxygen'][1:-1, -1:] = cell['oxygen'] + self.time_step * (
            ((cell['oxygen_diffusion_rate'] * (left['oxygen'] - 2 * cell['oxygen'] + left['oxygen'])) /
             self.model_parameters['spatial_step'] ** 2) +
            ((cell['oxygen_diffusion_rate'] * (above['oxygen'] - 2 * cell['oxygen'] + below['oxygen'])) /
             self.model_parameters['spatial_step'] ** 2) +
            (self.model_parameters['oxygen_from_source'] * cell['blood_vessel']) +
            (self.model_parameters['oxygen_uptake_from_bacteria'] * cell['oxygen'] *
             isinstance(cell['contents'], Bacterium))
        )
        if chemo:
            self.work_grid['chemotherapy'][1:-1, -1:] = cell['chemotherapy'] + self.time_step * (
                ((cell['chemotherapy_diffusion_rate'] * (left['chemotherapy'] - 2 * cell['chemotherapy'] +
                                                         left['chemotherapy'])) /
                 self.model_parameters['spatial_step'] ** 2) +
                ((cell['chemotherapy_diffusion_rate'] * (above['chemotherapy'] - 2 * cell['chemotherapy'] +
                                                         above['chemotherapy'])) /
                 self.model_parameters['spatial_step'] ** 2) +
                (self.model_parameters['chemotherapy_from_source'] * cell['blood_vessel']) +
                (self.model_parameters['chemotherapy_decay'] * cell['chemotherapy'])
            )
        self.work_grid['chemokine'][1:-1, -1:] = cell['chemokine'] + self.time_step * (
            ((self.model_parameters['chemokine_diffusion'] * (left['chemokine'] - 2 * cell['chemokine'] +
                                                              left['chemokine'])) /
             self.model_parameters['spatial_step'] ** 2) +
            ((self.model_parameters['chemokine_diffusion'] * (above['chemokine'] - 2 * cell['chemokine'] +
                                                              below['chemokine'])) /
             self.model_parameters['spatial_step'] ** 2) +
            self.model_parameters['chemokine_from_bacteria'] * isinstance(cell['contents'], Bacterium) +
            self.model_parameters['chemokine_from_macrophage'] *
            (isinstance(cell['contents'], Macrophage) and cell['contents'].state != 'resting') +
            self.model_parameters['chemokine_decay'] * cell['chemokine']
        )

        # BOTTOM ROW
        cell = self.grid[-1:, 1:-1]
        above = self.grid[-2:-1, 1:-1]
        left = self.grid[-1:, 0:-2]
        right = self.grid[-1:, 2:]
        self.work_grid['oxygen'][-1:, 1:-1] = cell['oxygen'] + self.time_step * (
            ((cell['oxygen_diffusion_rate'] * (above['oxygen'] - 2 * cell['oxygen'] + above['oxygen'])) /
             self.model_parameters['spatial_step'] ** 2) +
            ((cell['oxygen_diffusion_rate'] * (left['oxygen'] - 2 * cell['oxygen'] + right['oxygen'])) /
             self.model_parameters['spatial_step'] ** 2) +
            (self.model_parameters['oxygen_from_source'] * cell['blood_vessel']) +
            (self.model_parameters['oxygen_uptake_from_bacteria'] * cell['oxygen'] *
             isinstance(cell['contents'], Bacterium))
        )
        if chemo:
            self.work_grid['chemotherapy'][-1:, 1:-1] = cell['chemotherapy'] + self.time_step * (
                ((cell['chemotherapy_diffusion_rate'] * (above['chemotherapy'] - 2 * cell['chemotherapy'] +
                                                         above['chemotherapy'])) /
                 self.model_parameters['spatial_step'] ** 2) +
                ((cell['chemotherapy_diffusion_rate'] * (left['chemotherapy'] - 2 * cell['chemotherapy'] +
                                                         right['chemotherapy'])) /
                 self.model_parameters['spatial_step'] ** 2) +
                (self.model_parameters['chemotherapy_from_source'] * cell['blood_vessel']) +
                (self.model_parameters['chemotherapy_decay'] * cell['chemotherapy'])
            )
        self.work_grid['chemokine'][-1:, 1:-1] = cell['chemokine'] + self.time_step * (
            ((self.model_parameters['chemokine_diffusion'] * (above['chemokine'] - 2 * cell['chemokine'] +
                                                              above['chemokine'])) /
             self.model_parameters['spatial_step'] ** 2) +
            ((self.model_parameters['chemokine_diffusion'] * (left['chemokine'] - 2 * cell['chemokine'] +
                                                              right['chemokine'])) /
             self.model_parameters['spatial_step'] ** 2) +
            self.model_parameters['chemokine_from_bacteria'] * isinstance(cell['contents'], Bacterium) +
            self.model_parameters['chemokine_from_macrophage'] *
            (isinstance(cell['contents'], Macrophage) and cell['contents'].state != 'resting') +
            self.model_parameters['chemokine_decay'] * cell['chemokine']
        )

        #CORNERS
        # Top left
        cell = self.grid[:1, :1]
        below = self.grid[:1, 1:2]
        right = self.grid[1:2, :1]
        self.work_grid['oxygen'][:1, :1] = cell['oxygen'] + self.time_step * (
            ((cell['oxygen_diffusion_rate'] * (below['oxygen'] - 2 * cell['oxygen'] + below['oxygen'])) /
             self.model_parameters['spatial_step'] ** 2) +
            ((cell['oxygen_diffusion_rate'] * (right['oxygen'] - 2 * cell['oxygen'] + right['oxygen'])) /
             self.model_parameters['spatial_step'] ** 2) +
            (self.model_parameters['oxygen_from_source'] * cell['blood_vessel']) +
            (self.model_parameters['oxygen_uptake_from_bacteria'] * cell['oxygen'] *
             isinstance(cell['contents'], Bacterium))
        )
        if chemo:
            self.work_grid['chemotherapy'][:1, 1:-1] = cell['chemotherapy'] + self.time_step * (
                ((cell['chemotherapy_diffusion_rate'] * (below['chemotherapy'] - 2 * cell['chemotherapy'] +
                                                         below['chemotherapy'])) /
                 self.model_parameters['spatial_step'] ** 2) +
                ((cell['chemotherapy_diffusion_rate'] * (right['chemotherapy'] - 2 * cell['chemotherapy'] +
                                                         right['chemotherapy'])) /
                 self.model_parameters['spatial_step'] ** 2) +
                (self.model_parameters['chemotherapy_from_source'] * cell['blood_vessel']) +
                (self.model_parameters['chemotherapy_decay'] * cell['chemotherapy'])
            )
        self.work_grid['chemokine'][:1, 1:-1] = cell['chemokine'] + self.time_step * (
            ((self.model_parameters['chemokine_diffusion'] * (below['chemokine'] - 2 * cell['chemokine'] +
                                                              below['chemokine'])) /
             self.model_parameters['spatial_step'] ** 2) +
            ((self.model_parameters['chemokine_diffusion'] * (right['chemokine'] - 2 * cell['chemokine'] +
                                                              right['chemokine'])) /
             self.model_parameters['spatial_step'] ** 2) +
            self.model_parameters['chemokine_from_bacteria'] * isinstance(cell['contents'], Bacterium) +
            self.model_parameters['chemokine_from_macrophage'] *
            (isinstance(cell['contents'], Macrophage) and cell['contents'].state != 'resting') +
            self.model_parameters['chemokine_decay'] * cell['chemokine']
        )
        # Top right
        cell = self.grid[:1, -1:]
        below = self.grid[1:2, -1:]
        left = self.grid[:1, -2:-1]
        self.work_grid['oxygen'][:1, :1] = cell['oxygen'] + self.time_step * (
            ((cell['oxygen_diffusion_rate'] * (below['oxygen'] - 2 * cell['oxygen'] + below['oxygen'])) /
             self.model_parameters['spatial_step'] ** 2) +
            ((cell['oxygen_diffusion_rate'] * (left['oxygen'] - 2 * cell['oxygen'] + left['oxygen'])) /
             self.model_parameters['spatial_step'] ** 2) +
            (self.model_parameters['oxygen_from_source'] * cell['blood_vessel']) +
            (self.model_parameters['oxygen_uptake_from_bacteria'] * cell['oxygen'] *
             isinstance(cell['contents'], Bacterium))
        )
        if chemo:
            self.work_grid['chemotherapy'][:1, 1:-1] = cell['chemotherapy'] + self.time_step * (
                ((cell['chemotherapy_diffusion_rate'] * (below['chemotherapy'] - 2 * cell['chemotherapy'] +
                                                         below['chemotherapy'])) /
                 self.model_parameters['spatial_step'] ** 2) +
                ((cell['chemotherapy_diffusion_rate'] * (left['chemotherapy'] - 2 * cell['chemotherapy'] +
                                                         left['chemotherapy'])) /
                 self.model_parameters['spatial_step'] ** 2) +
                (self.model_parameters['chemotherapy_from_source'] * cell['blood_vessel']) +
                (self.model_parameters['chemotherapy_decay'] * cell['chemotherapy'])
            )
        self.work_grid['chemokine'][:1, 1:-1] = cell['chemokine'] + self.time_step * (
            ((self.model_parameters['chemokine_diffusion'] * (below['chemokine'] - 2 * cell['chemokine'] +
                                                              below['chemokine'])) /
             self.model_parameters['spatial_step'] ** 2) +
            ((self.model_parameters['chemokine_diffusion'] * (left['chemokine'] - 2 * cell['chemokine'] +
                                                              left['chemokine'])) /
             self.model_parameters['spatial_step'] ** 2) +
            self.model_parameters['chemokine_from_bacteria'] * isinstance(cell['contents'], Bacterium) +
            self.model_parameters['chemokine_from_macrophage'] *
            (isinstance(cell['contents'], Macrophage) and cell['contents'].state != 'resting') +
            self.model_parameters['chemokine_decay'] * cell['chemokine']
        )
        # Bottom left
        cell = self.grid[-1:, :1]
        above = self.grid[-2:-1, :1]
        right = self.grid[-1:, 1:2]
        self.work_grid['oxygen'][:1, :1] = cell['oxygen'] + self.time_step * (
            ((cell['oxygen_diffusion_rate'] * (above['oxygen'] - 2 * cell['oxygen'] + above['oxygen'])) /
             self.model_parameters['spatial_step'] ** 2) +
            ((cell['oxygen_diffusion_rate'] * (right['oxygen'] - 2 * cell['oxygen'] + right['oxygen'])) /
             self.model_parameters['spatial_step'] ** 2) +
            (self.model_parameters['oxygen_from_source'] * cell['blood_vessel']) +
            (self.model_parameters['oxygen_uptake_from_bacteria'] * cell['oxygen'] *
             isinstance(cell['contents'], Bacterium))
        )
        if chemo:
            self.work_grid['chemotherapy'][:1, 1:-1] = cell['chemotherapy'] + self.time_step * (
                ((cell['chemotherapy_diffusion_rate'] * (above['chemotherapy'] - 2 * cell['chemotherapy'] +
                                                         above['chemotherapy'])) /
                 self.model_parameters['spatial_step'] ** 2) +
                ((cell['chemotherapy_diffusion_rate'] * (right['chemotherapy'] - 2 * cell['chemotherapy'] +
                                                         right['chemotherapy'])) /
                 self.model_parameters['spatial_step'] ** 2) +
                (self.model_parameters['chemotherapy_from_source'] * cell['blood_vessel']) +
                (self.model_parameters['chemotherapy_decay'] * cell['chemotherapy'])
            )
        self.work_grid['chemokine'][:1, 1:-1] = cell['chemokine'] + self.time_step * (
            ((self.model_parameters['chemokine_diffusion'] * (above['chemokine'] - 2 * cell['chemokine'] +
                                                              above['chemokine'])) /
             self.model_parameters['spatial_step'] ** 2) +
            ((self.model_parameters['chemokine_diffusion'] * (right['chemokine'] - 2 * cell['chemokine'] +
                                                              right['chemokine'])) /
             self.model_parameters['spatial_step'] ** 2) +
            self.model_parameters['chemokine_from_bacteria'] * isinstance(cell['contents'], Bacterium) +
            self.model_parameters['chemokine_from_macrophage'] *
            (isinstance(cell['contents'], Macrophage) and cell['contents'].state != 'resting') +
            self.model_parameters['chemokine_decay'] * cell['chemokine']
        )
        # Bottom right
        cell = self.grid[-1:, -1:]
        above = self.grid[-2:-1, -1:]
        left = self.grid[-1:, -2:-1]
        self.work_grid['oxygen'][:1, :1] = cell['oxygen'] + self.time_step * (
            ((cell['oxygen_diffusion_rate'] * (above['oxygen'] - 2 * cell['oxygen'] + above['oxygen'])) /
             self.model_parameters['spatial_step'] ** 2) +
            ((cell['oxygen_diffusion_rate'] * (left['oxygen'] - 2 * cell['oxygen'] + left['oxygen'])) /
             self.model_parameters['spatial_step'] ** 2) +
            (self.model_parameters['oxygen_from_source'] * cell['blood_vessel']) +
            (self.model_parameters['oxygen_uptake_from_bacteria'] * cell['oxygen'] *
             isinstance(cell['contents'], Bacterium))
        )
        if chemo:
            self.work_grid['chemotherapy'][:1, 1:-1] = cell['chemotherapy'] + self.time_step * (
                ((cell['chemotherapy_diffusion_rate'] * (above['chemotherapy'] - 2 * cell['chemotherapy'] +
                                                         above['chemotherapy'])) /
                 self.model_parameters['spatial_step'] ** 2) +
                ((cell['chemotherapy_diffusion_rate'] * (left['chemotherapy'] - 2 * cell['chemotherapy'] +
                                                         left['chemotherapy'])) /
                 self.model_parameters['spatial_step'] ** 2) +
                (self.model_parameters['chemotherapy_from_source'] * cell['blood_vessel']) +
                (self.model_parameters['chemotherapy_decay'] * cell['chemotherapy'])
            )
        self.work_grid['chemokine'][:1, 1:-1] = cell['chemokine'] + self.time_step * (
            ((self.model_parameters['chemokine_diffusion'] * (above['chemokine'] - 2 * cell['chemokine'] +
                                                              above['chemokine'])) /
             self.model_parameters['spatial_step'] ** 2) +
            ((self.model_parameters['chemokine_diffusion'] * (left['chemokine'] - 2 * cell['chemokine'] +
                                                              left['chemokine'])) /
             self.model_parameters['spatial_step'] ** 2) +
            self.model_parameters['chemokine_from_bacteria'] * isinstance(cell['contents'], Bacterium) +
            self.model_parameters['chemokine_from_macrophage'] *
            (isinstance(cell['contents'], Macrophage) and cell['contents'].state != 'resting') +
            self.model_parameters['chemokine_decay'] * cell['chemokine']
        )

        if not chemo:
            self.work_grid['chemotherapy'] = np.zeros(self.grid.shape,dtype=float)

    def update_agents(self):
        pass


class Bacterium(Agent):

    def __init__(self, metabolism):
        Agent.__init__(self)
        self.metabolism = metabolism
        self.resting = False

    def output_code(self):
        # 1.0 == Fast, 2.0 == slow
        code = 1.0 + (self.metabolism == 'slow')
        # Add .5 if resting
        if self.resting:
            code += 0.5
        return code


class Caseum(Agent):

    def __init__(self):
        Agent.__init__(self)

    def output_code(self):
        return 100.0


class TCell(Agent):

    def __init__(self):
        Agent.__init__(self)

    def output_code(self):
        return 3.0


class Macrophage(Agent):

    def __init__(self, state):
        Agent.__init__(self)
        self.state = state
        self.intracellular_bacteria = 0

    def output_code(self):
        code = 4.0
        if self.state == 'active':
            code += 1.0
        elif self.state == 'infected':
            code += 2.0
        elif self.state == 'chronically_infected':
            code += 3.0
        return code


# ---------------------------------------- Runner -------------------------------------------------------


def initialise(config, total_shape):
    """
    Convert the defined lists/settings into a series of coordinates for initialisation of the automaton
    :param config:
    :param total_shape:
    :return:
    """

    available_addresses = []
    for a in itertools.product(range(total_shape[0]), range(total_shape[1])):
        available_addresses.append(a)

    # BLOOD_VESSELS
    blood_vessels_method = config.get("InitialiseSection", "blood_vessels")
    blood_vessel_addresses = []

    if blood_vessels_method == 'hard_code':
        bv_list = config.get("InitialiseSection", 'blood_vessels_hard_code').split('/')
        for b in bv_list:
            address = tuple(int(c) for c in b.split(","))
            available_addresses.remove(address)
            blood_vessel_addresses.append(address)
    elif blood_vessels_method == 'from_file':
        path = config.get("InitialiseSection", 'blood_vessels_from_file')
        # Add the values to the list
        list_of_vessels = [line.rstrip('\n') for line in open(path)]
        for index in range(len(list_of_vessels)):
            if float(list_of_vessels[index]) > 0.0:
                try:
                    address = np.unravel_index(index, total_shape)
                except ValueError:
                    raise Exception(index, total_shape)
                blood_vessel_addresses.append(address)
                available_addresses.remove(address)
    elif blood_vessels_method == 'random':
        number = config.getint("InitialiseSection", "blood_vessels_random_number")
        assert len(available_addresses) > number
        for i in range(number):
            address = available_addresses.pop(np.random.randint(0, len(available_addresses)))
            blood_vessel_addresses.append(address)

    # FAST BACTERIA
    bacteria_method = config.get("InitialiseSection", "bacteria")
    fast_addresses = []
    slow_addresses = []

    if bacteria_method == 'hard_code':
        fast_list = config.get("InitialiseSection", "bacteria_fast_hard_code").split("/")

        assert len(available_addresses) > len(fast_list)
        for a in fast_list:
            address = tuple(int(c) for c in a.split(","))
            if address in available_addresses:
                fast_addresses.append(address)
                available_addresses.remove(address)
            else:
                # TODO - avoid conflict
                pass

        slow_list = config.get("InitialiseSection", "bacteria_slow_hard_code").split("/")
        assert len(available_addresses) > len(slow_list)
        for a in slow_list:
            address = tuple(int(c) for c in a.split(","))
            if address not in blood_vessel_addresses:
                slow_addresses.append(address)
                available_addresses.remove(address)
            else:
                # TODO - avoid conflict
                pass
    elif bacteria_method == 'random':
        number_fast = config.getint("InitialiseSection", "bacteria_fast_random_number")
        assert len(available_addresses) > number_fast
        for i in range(number_fast):
            address = available_addresses.pop(np.random.randint(0, len(available_addresses)))
            fast_addresses.append(address)

        number_slow = config.getint("InitialiseSection", "bacteria_slow_random_number")
        assert len(available_addresses) > number_slow
        for i in range(number_slow):
            address = available_addresses.pop(np.random.randint(0, len(available_addresses)))
            slow_addresses.append(address)

    # MACROPHAGES
    macrophage_method = config.get("InitialiseSection", "macrophages")
    macrophage_addresses = []
    if macrophage_method == 'random':
        number = config.getint("InitialiseSection", "macrophages_random_number")

        # Make sure there's enough room
        assert len(available_addresses) > number

        for i in range(number):
            address = available_addresses.pop(np.random.randint(0, len(available_addresses)))
            macrophage_addresses.append(address)
    # TODO: hard-code macrophages?

    return blood_vessel_addresses, fast_addresses, slow_addresses, macrophage_addresses


if __name__ == '__main__':

    import ConfigParser
    import os
    import time

    print '------------------------'
    print 'TB Simulation Automaton'
    print '------------------------'
    whole_start_time = time.time()
    print "Begin:   {", whole_start_time, "}"

    config = ConfigParser.RawConfigParser()
    if not config.read('config.properties'):
        raise IOError("Config file (config.properties) not found")

    # LOAD PARAMETERS
    parameters = {}
    # Get all options in parameters section and add to the dictionary
    for i in config.options("ParametersSection"):
        parameters[i] = config.getfloat("ParametersSection", i)

    # TIME PARAMETERS
    time_parameters = {}
    # Get all options in time parameters section
    for i in config.options("TimeParametersSection"):
        time_parameters[i] = config.getfloat("TimeParametersSection", i)

    # LOAD GRID ATTRIBUTES
    total_shape = [int(a) for a in config.get("GridSection", "total_shape").split(",")]

    # LOAD RUN PARAMETERS
    output_location = config.get("RunParametersSection", "output_location")
    if not os.path.exists(output_location):
        os.makedirs(output_location)

    # LOAD INITIALISATION
    blood_vessels, fast_bacteria, slow_bacteria, macrophages = initialise(config, total_shape)

    automaton = TBAutomaton(total_shape, time_parameters, parameters, output_location,
                            blood_vessels, macrophages, fast_bacteria, slow_bacteria)

    automaton.run()

    whole_end_time = time.time()
    print "End:     {", whole_end_time, "}"
    print "Duration: ", whole_end_time - whole_start_time
