from CAPEAutomaton import *

# ---------------------------------------- TB -------------------------------------------------------


class TBAutomaton(Automaton):

    def __init__(self, shape, time_parameters, model_parameters, output_location,
                 blood_vessel_addresses, initial_macrophage_addresses,
                 initial_fast_bacteria_addresses, initial_slow_bacteria_addresses):
        attributes = ['oxygen', 'chemotherapy', 'chemokine', 'contents', 'oxygen_diffusion_rate',
                      'chemotherapy_diffusion_rate', 'blood_vessel']
        formats = ['float', 'float', 'float', 'object', 'float', 'float', 'float']

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

        self.values_to_record = ["fast_bacteria", "fast_bacteria_resting", "slow_bacteria", "slow_bacteria_resting",
                                 "intracellular_bac", "total_bacteria",
                                 "resting_macrophages", "active_macrophages", "infected_macrophages",
                                 "chron_infected_macrophages", "total_macrophages",
                                 "t_cells", "caseum"]

        Automaton.__init__(self, shape, attributes, formats, time_parameters, model_parameters, output_location,
                           self.values_to_record, initialisation)

        self.chemo_schedule1_start = np.random.randint(self.model_parameters['chemotherapy_schedule1_start_lower'],
                                                       self.model_parameters['chemotherapy_schedule1_start_upper'])

    def timestep_output(self):
        print "t = ", self.time * self.time_step, "Bac = ", len(self.bacteria)

    def record_counts(self):
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
        self.diffusion_pre_process()
        chemo = (self.chemo_schedule1_start / self.time_step) <= self.time < \
                (self.model_parameters['chemotherapy_schedule1_end'] / self.time_step) or \
                self.model_parameters['chemotherapy_schedule2_start'] / self.time_step <= self.time
        self.diffusion(chemo)

    def diffusion_pre_process(self):
        affected_addresses = []

        for caseum_address in self.caseum_addresses:
            for d in range(self.model_parameters['caseum_distance_to_reduce_diffusion']):
                neighbours = self.moore_neighbours(caseum_address,d)
                affected_addresses += neighbours

        # Affected addresses is now a list of all address within the range of cells with caseum
        counted = Counter(affected_addresses)
        for address in counted:
            if counted[address] >= self.model_parameters['caseum_threshold_to_reduce_diffusion']:
                self.grid[address]['oxygen_diffusion_rate'] = self.model_parameters['oxygen_diffusion'] / \
                                                              self.model_parameters['diffusion_caseum_reduction']
                self.grid[address]['chemotherapy_diffusion_rate'] = self.model_parameters['chemotherapy_diffusion'] / \
                                                        self.model_parameters['diffusion_caseum_reduction']
                if self.grid[address]['blood_vessel'] > 0.0:
                    self.grid[address]['blood_vessel'] /= self.model_parameters['diffusion_caseum_reduction']

    def diffusion(self, chemo):
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


if __name__ == '__main__':

    time_params={}
    time_params['initial_time'] = 0.0
    time_params['time_step'] = 0.001
    time_params['time_limit'] = 1.0
    time_params['interval_to_record_grid'] = 100.0
    time_params['interval_to_record_counts'] = 1.0

    model_params = {}
    model_params['spatial_step'] = 0.1
    model_params['oxygen_from_source'] = 0.1
    model_params['blood_vessel_value'] = 0.1
    model_params['oxygen_uptake_from_bacteria'] = 0.1
    model_params['chemotherapy_from_source'] = 0.1
    model_params['chemotherapy_decay'] = 0.1
    model_params['chemokine_diffusion'] = 0.1
    model_params['chemokine_from_bacteria'] = 0.1
    model_params['chemokine_from_macrophage'] = 0.1
    model_params['chemokine_decay'] = 0.1
    model_params['oxygen_diffusion'] = 0.1
    model_params['chemotherapy_diffusion'] = 0.1
    model_params['chemotherapy_schedule1_start_lower'] = 1.0
    model_params['chemotherapy_schedule1_start_upper'] = 2.0
    model_params['chemotherapy_schedule1_end'] = 2.0
    model_params['chemotherapy_schedule2_start'] = 2.0

    start_time = time.time()
    tba = TBAutomaton((100, 100), time_params, model_params, "output", [], [], [(1,1)], [(0,0)])
    tba.run()
