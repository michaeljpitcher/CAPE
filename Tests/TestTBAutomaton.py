import os
import shutil
import unittest

from TBAutomaton.TBAutomaton import *


class TBAutomatonTestCase(unittest.TestCase):
    def setUp(self):
        self.shape = (10, 10)
        self.time_params = {}
        self.time_params['initial_time'] = 0.0
        self.time_params['time_step'] = 0.1
        self.time_params['time_limit'] = 1.0
        self.model_params = {}
        self.model_params['chemotherapy_schedule1_start_lower'] = 1.0
        self.model_params['chemotherapy_schedule1_start_upper'] = 2.0
        self.model_params['blood_vessel_value'] = 1.0
        self.model_params['initial_oxygen'] = 1.0
        self.model_params['oxygen_diffusion'] = 1.0
        self.model_params['chemotherapy_diffusion'] = 1.0
        self.model_params['bacteria_replication_fast_upper'] = 10.0
        self.model_params['bacteria_replication_fast_lower'] = 9.0
        self.model_params['bacteria_replication_slow_upper'] = 20.0
        self.model_params['bacteria_replication_slow_lower'] = 19.0
        self.model_params['bacteria_threshold_for_t_cells'] = 100
        self.model_params['t_cell_recruitment_probability'] = 0
        self.model_params['chemokine_scale_for_t_cell_recruitment'] = 1.1
        self.model_params['bacteria_threshold_for_macrophage_recruitment'] = 100
        self.model_params['chemokine_scale_for_macrophage_recruitment_below_threshold'] = 1.01
        self.model_params['macrophage_recruitment_probability'] = 0.0
        self.model_params['chemotherapy_scale_for_kill_fast_bacteria'] = 1.01
        self.model_params['chemotherapy_scale_for_kill_slow_bacteria'] = 1.01
        self.model_params['chemotherapy_scale_for_kill_macrophage'] = 1.01
        self.model_params['t_cell_movement_time'] = 999999
        self.model_params['t_cell_age_threshold'] = 999999
        self.model_params['t_cell_random_move_probability'] = 0.0
        self.model_params['t_cell_kills_macrophage_probability'] = 0.0
        self.model_params['resting_macrophage_age_limit'] = 100000.0
        self.model_params['resting_macrophage_movement_time'] = 100000.0
        self.model_params['prob_resting_macrophage_random_move'] = 0.0
        self.model_params['minimum_chemokine_for_resting_macrophage_movement'] = 101.0
        self.model_params['active_macrophage_age_limit'] = 1000000.0
        self.model_params['active_macrophage_movement_time'] = 1000000.0
        self.model_params['prob_active_macrophage_kill_fast_bacteria'] = 0.0
        self.model_params['prob_active_macrophage_kill_slow_bacteria'] = 0.0
        self.model_params['infected_macrophage_age_limit'] = 1000000.0
        self.model_params['infected_macrophage_movement_time'] = 1000000.0
        self.model_params['chronically_infected_macrophage_age_limit'] = 1000000.0
        self.model_params['chronically_infected_macrophage_movement_time'] = 1000000.0
        self.model_params['chemokine_scale_for_macrophage_activation'] = 101.0
        self.model_params['chemokine_scale_for_macrophage_deactivation'] = 0.0
        self.model_params['bacteria_to_burst_macrophage'] = 999999
        self.model_params['oxygen_scale_for_metabolism_change_to_slow'] = -1.0
        self.model_params['oxygen_scale_for_metabolism_change_to_fast'] = 101.0


        self.bv = [(1, 1), (2, 3), (3, 5)]
        self.macs = [(9, 9), (8, 8), (7, 7), (6, 6)]
        self.fb = [(8, 1), (8, 2), (8, 3)]
        self.sb = [(1, 7), (2, 7)]

        self.output_loc = 'test_output'
        if not os.path.exists(self.output_loc):
            os.makedirs(self.output_loc)
        self.automaton = TBAutomaton(self.shape, self.time_params, self.model_params, self.output_loc,
                                     self.bv, self.macs, self.fb, self.sb)

    def tearDown(self):
        # Close output files and delete
        self.automaton.close_files()
        shutil.rmtree(self.output_loc)

    def test_initialise(self):
        atts = ['oxygen', 'chemotherapy', 'chemokine', 'contents', 'oxygen_diffusion_rate',
                'chemotherapy_diffusion_rate', 'blood_vessel']
        self.assertItemsEqual(atts, self.automaton.attributes)
        self.assertItemsEqual(self.model_params, self.automaton.model_parameters)
        self.assertEqual(self.automaton.time, self.time_params['initial_time'])
        self.assertEqual(self.automaton.time_step, self.time_params['time_step'])
        self.assertEqual(self.automaton.time_limit, self.time_params['time_limit'] / self.time_params['time_step'])

        self.assertSequenceEqual(self.automaton.grid.shape, self.shape)
        self.assertSequenceEqual(self.automaton.work_grid.shape, self.shape)

        for x in range(self.shape[0]):
            for y in range(self.shape[1]):
                if (x,y) not in self.bv and (x,y) not in self.macs and (x,y) not in self.fb and (x,y) not in self.sb:
                    self.assertEqual(self.automaton.grid[(x,y)]['contents'], 0.0)
                elif (x,y) in self.bv:
                    self.assertEqual(self.automaton.grid[(x, y)]['contents'], self.model_params['blood_vessel_value'])
                    self.assertEqual(self.automaton.grid[(x, y)]['blood_vessel'], self.model_params['blood_vessel_value'])
                elif (x, y) in self.macs:
                    self.assertTrue(isinstance(self.automaton.grid[(x, y)]['contents'], Macrophage))
                    self.assertTrue(self.automaton.grid[(x, y)]['contents'] in self.automaton.agents)
                    self.assertTrue(self.automaton.grid[(x, y)]['contents'] in self.automaton.macrophages)
                elif (x,y) in self.fb:
                    self.assertTrue(isinstance(self.automaton.grid[(x,y)]['contents'], Bacterium))
                    self.assertEqual(self.automaton.grid[(x,y)]['contents'].metabolism, 'fast')
                    self.assertTrue(self.automaton.grid[(x, y)]['contents'] in self.automaton.agents)
                    self.assertTrue(self.automaton.grid[(x, y)]['contents'] in self.automaton.bacteria)
                elif (x,y) in self.sb:
                    self.assertTrue(isinstance(self.automaton.grid[(x,y)]['contents'], Bacterium))
                    self.assertEqual(self.automaton.grid[(x,y)]['contents'].metabolism, 'slow')
                    self.assertTrue(self.automaton.grid[(x, y)]['contents'] in self.automaton.agents)
                    self.assertTrue(self.automaton.grid[(x, y)]['contents'] in self.automaton.bacteria)

        self.assertEqual(self.automaton.chemo_schedule1_start, 1.0)
        self.assertEqual(len(self.automaton.agents), len(self.macs) + len(self.fb) + len(self.sb))
        self.assertItemsEqual(self.automaton.blood_vessel_addresses, self.bv)
        self.assertEqual(len(self.automaton.t_cells), 0.0)
        self.assertEqual(len(self.automaton.caseum_addresses), 0.0)

    def test_diffusion_pre_process_no_caseum(self):
        self.automaton.diffusion_pre_process()
        for x in range(10):
            for y in range(10):
                cell = self.automaton.grid[(x,y)]
                self.assertEqual(cell['oxygen_diffusion_rate'], self.model_params['oxygen_diffusion'])
                self.assertEqual(cell['chemotherapy_diffusion_rate'], self.model_params['chemotherapy_diffusion'])

    def test_diffusion_pre_process_caseum_no_reduction(self):
        self.automaton.model_parameters['caseum_distance_to_reduce_diffusion'] = 1
        self.automaton.model_parameters['caseum_threshold_to_reduce_diffusion'] = 2

        cas = Caseum((4,1))
        self.automaton.grid[(4,1)]['contents'] = cas
        self.automaton.caseum_addresses.append((4,1))

        self.automaton.diffusion_pre_process()
        for x in range(10):
            for y in range(10):
                cell = self.automaton.grid[(x,y)]
                self.assertEqual(cell['oxygen_diffusion_rate'], self.model_params['oxygen_diffusion'])
                self.assertEqual(cell['chemotherapy_diffusion_rate'], self.model_params['chemotherapy_diffusion'])

    def test_diffusion_pre_process_caseum_reduction(self):
        self.automaton.model_parameters['caseum_distance_to_reduce_diffusion'] = 1
        self.automaton.model_parameters['caseum_threshold_to_reduce_diffusion'] = 1

        self.automaton.model_parameters['diffusion_caseum_reduction'] = 2.0

        cas = Caseum((4,1))
        self.automaton.grid[(4,1)]['contents'] = cas
        self.automaton.caseum_addresses.append((4,1))

        expected_reductions = [(3,0),(3,1),(3,2),(4,0),(4,2),(5,0),(5,1),(5,2)]

        self.automaton.diffusion_pre_process()
        for x in range(10):
            for y in range(10):
                cell = self.automaton.grid[(x,y)]
                if (x,y) in expected_reductions:
                    self.assertEqual(cell['oxygen_diffusion_rate'], self.model_params['oxygen_diffusion'] / self.automaton.model_parameters['diffusion_caseum_reduction'])
                    self.assertEqual(cell['chemotherapy_diffusion_rate'], self.model_params['chemotherapy_diffusion'] / self.automaton.model_parameters['diffusion_caseum_reduction'])
                else:
                    self.assertEqual(cell['oxygen_diffusion_rate'], self.model_params['oxygen_diffusion'])
                    self.assertEqual(cell['chemotherapy_diffusion_rate'], self.model_params['chemotherapy_diffusion'])

    def test_record_counts(self):
        # BASE
        self.automaton.record_counts()
        # Increment all counts
        self.automaton.time = 10.0
        fast_bac = Bacterium((9,9), 'fast')
        self.automaton.bacteria.append(fast_bac)
        fast_bac_rest = Bacterium((8,9), 'fast')
        fast_bac_rest.resting = True
        self.automaton.bacteria.append(fast_bac_rest)
        slow_bac = Bacterium((7,9), 'slow')
        self.automaton.bacteria.append(slow_bac)
        slow_bac_rest = Bacterium((6,9), 'slow')
        slow_bac_rest.resting = True
        self.automaton.bacteria.append(slow_bac_rest)
        self.automaton.macrophages[0].intracellular_bacteria = 1
        rest_mac = Macrophage((5,9), 'resting')
        self.automaton.macrophages.append(rest_mac)
        act_mac = Macrophage((4, 9), 'active')
        self.automaton.macrophages.append(act_mac)
        inf_mac = Macrophage((3, 9), 'infected')
        self.automaton.macrophages.append(inf_mac)
        chrinf_mac = Macrophage((2, 9), 'chronically_infected')
        self.automaton.macrophages.append(chrinf_mac)
        tcell = TCell((1,9))
        self.automaton.t_cells.append(tcell)
        self.automaton.caseum_addresses.append((0,9))
        self.automaton.record_counts()

        self.automaton.close_files()
        count_file = open(self.output_loc + '/counts.csv', 'rb')
        reader = csv.reader(count_file, delimiter=',')
        counter = 0
        for row in reader:
            if counter == 0:
                self.assertSequenceEqual(row, ["timestep", "fast_bacteria", "fast_bacteria_resting", "slow_bacteria",
                                               "slow_bacteria_resting", "intracellular_bac", "total_bacteria",
                                               "resting_macrophages", "active_macrophages", "infected_macrophages",
                                               "chron_infected_macrophages", "total_macrophages", "t_cells", "caseum"])
            elif counter == 1:
                self.assertSequenceEqual(row, ['0.0', str(len(self.fb)), '0', str(len(self.sb)), '0', '0',
                                               str(len(self.fb) + len(self.sb)), str(len(self.macs)), '0', '0', '0',
                                               str(len(self.macs)), '0', '0'])
            elif counter == 2:
                self.assertSequenceEqual(row, [str(10.0 * 0.1), str(len(self.fb) + 1), str(1), str(len(self.sb)+1), str(1), str(1),
                                               str(len(self.fb) + len(self.sb) + 5), str(len(self.macs) + 1), str(1), str(1), str(1),
                                               str(len(self.macs) + 4), str(1), str(1)])
            counter += 1

    def test_bacterium_replicate_fast_not_slow(self):
        self.automaton.time = 50.0
        self.automaton.model_parameters['bacteria_replication_fast_upper'] = 6.0
        self.automaton.model_parameters['bacteria_replication_fast_lower'] = 5.0
        self.automaton.model_parameters['bacteria_replication_slow_upper'] = 100.0
        self.automaton.model_parameters['bacteria_replication_slow_lower'] = 99.0
        events = self.automaton.bacteria_replication()
        self.assertEqual(len(events), len(self.fb))
        for event in events:
            self.assertTrue(isinstance(event, BacteriumReplication))
            self.assertTrue(event.original_bac_address in self.fb)
            self.assertEqual(self.automaton.grid[event.new_bac_address]['contents'], 0.0)
            self.assertEqual(event.new_metabolism, 'fast')

    def test_bacterium_replicate_slow_not_fast(self):
        self.automaton.time = 50.0
        self.automaton.model_parameters['bacteria_replication_fast_upper'] = 3.0
        self.automaton.model_parameters['bacteria_replication_fast_lower'] = 2.0
        self.automaton.model_parameters['bacteria_replication_slow_upper'] = 6.0
        self.automaton.model_parameters['bacteria_replication_slow_lower'] = 5.0
        events = self.automaton.bacteria_replication()
        self.assertEqual(len(events), len(self.sb))
        for event in events:
            self.assertTrue(isinstance(event, BacteriumReplication))
            self.assertTrue(event.original_bac_address in self.sb)
            self.assertEqual(self.automaton.grid[event.new_bac_address]['contents'], 0.0)
            self.assertEqual(event.new_metabolism, 'slow')

    def test_replicate_no_room(self):
        self.automaton.bacteria = []
        bac = Bacterium((1,1), 'fast')
        self.automaton.bacteria.append(bac)

        for x in range(self.shape[0]):
            for y in range(self.shape[1]):
                if x == 1 and y == 1:
                    self.automaton.grid[(x, y)]['contents'] = bac
                else:
                    self.automaton.grid[(x, y)]['contents'] = Caseum((x, y))

        events = self.automaton.bacteria_replication()
        self.assertEqual(len(events), 1)
        self.assertTrue(isinstance(events[0], BacteriumStateChange))
        self.assertEqual(events[0].attribute, 'resting')
        self.assertEqual(events[0].value, True)

    def test_scales(self):
        self.automaton.max_oxygen = 0.0
        self.automaton.max_chemotherapy = 0.0
        self.automaton.max_chemokine = 0.0
        self.assertEqual(self.automaton.oxygen_scale((0, 0)), 0.0)
        self.assertEqual(self.automaton.chemotherapy_scale((0, 0)), 0.0)
        self.assertEqual(self.automaton.chemokine_scale((0, 0)), 0.0)

        self.automaton.max_oxygen = 100.0
        self.automaton.max_chemotherapy = 75.0
        self.automaton.max_chemokine = 40.0

        self.automaton.grid[(0, 0)]['oxygen'] = 56.0
        self.automaton.grid[(4, 7)]['chemotherapy'] = 36.0
        self.automaton.grid[(8, 8)]['chemokine'] = 1.0

        self.assertEqual(self.automaton.oxygen_scale((0, 0)), 0.56)
        self.assertAlmostEqual(self.automaton.chemotherapy_scale((4, 7)), 0.48)
        self.assertEqual(self.automaton.chemokine_scale((8, 8)), 0.025)

    def test_t_cell_recruitment_positive(self):
        self.automaton.model_parameters['bacteria_threshold_for_t_cells'] = 0
        self.automaton.model_parameters['t_cell_recruitment_probability'] = 100
        self.automaton.model_parameters['chemokine_scale_for_t_cell_recruitment'] = -0.1

        events = self.automaton.t_cell_recruitment()
        self.assertEqual(len(events), 3)

        self.automaton.model_parameters['bacteria_threshold_for_t_cells'] = 2
        self.automaton.model_parameters['t_cell_recruitment_probability'] = 100
        self.automaton.model_parameters['chemokine_scale_for_t_cell_recruitment'] = -0.1

        events = self.automaton.t_cell_recruitment()
        self.assertEqual(len(events), 3)

        self.automaton.model_parameters['bacteria_threshold_for_t_cells'] = 0
        self.automaton.model_parameters['t_cell_recruitment_probability'] = 100
        self.automaton.model_parameters['chemokine_scale_for_t_cell_recruitment'] = 0.5

        self.automaton.max_chemokine = 100.0
        self.automaton.grid[(0, 1)]['chemokine'] = 66.0
        self.automaton.grid[(1, 3)]['chemokine'] = 66.0
        self.automaton.grid[(4, 5)]['chemokine'] = 66.0

        events = self.automaton.t_cell_recruitment()
        self.assertEqual(len(events), 3)

        bv_addresses = [(e.blood_vessel_address, e.new_t_cell_address) for e in events]
        self.assertItemsEqual(bv_addresses, [((1,1),(0,1)), ((2,3),(1,3)), ((3,5), (4,5))])

    def test_t_cell_recruit_negative(self):
        events = self.automaton.t_cell_recruitment()
        self.assertEqual(len(events), 0)

    def test_total_bacteria(self):
        # 5 initial
        self.assertEqual(self.automaton.total_bacteria(), len(self.fb) + len(self.sb))
        # Add intracellular
        self.automaton.macrophages[0].intracellular_bacteria = 3
        self.automaton.macrophages[1].intracellular_bacteria = 2
        self.automaton.macrophages[2].intracellular_bacteria = 1
        self.assertEqual(self.automaton.total_bacteria(), len(self.fb) + len(self.sb) + 6)

    def test_macrophage_recruitment_positive_below_threshold(self):
        self.model_params['bacteria_threshold_for_macrophage_recruitment'] = 100
        self.model_params['chemokine_scale_for_macrophage_recruitment_below_threshold'] = -0.1
        self.model_params['macrophage_recruitment_probability'] = 100.0
        events = self.automaton.macrophage_recruitment()
        self.assertEqual(len(events), len(self.bv))

    def test_macrophage_recruitment_positive_above_threshold(self):
        self.model_params['bacteria_threshold_for_macrophage_recruitment'] = 0
        self.model_params['chemokine_scale_for_macrophage_recruitment_above_threshold'] = -0.1
        self.model_params['macrophage_recruitment_probability'] = 100.0
        events = self.automaton.macrophage_recruitment()
        self.assertEqual(len(events), len(self.bv))

    def test_macrophage_recruitment_negative_prob(self):
        self.model_params['bacteria_threshold_for_macrophage_recruitment'] = 100
        self.model_params['chemokine_scale_for_macrophage_recruitment_below_threshold'] = -0.1
        self.model_params['macrophage_recruitment_probability'] = 0.0
        events = self.automaton.macrophage_recruitment()
        self.assertEqual(len(events), 0)

        self.model_params['bacteria_threshold_for_macrophage_recruitment'] = 0
        self.model_params['chemokine_scale_for_macrophage_recruitment_above_threshold'] = -0.1
        self.model_params['macrophage_recruitment_probability'] = 0.0
        events = self.automaton.macrophage_recruitment()
        self.assertEqual(len(events), 0)

    def test_macrophage_recruitment_negative_scale(self):

        for x in range(self.shape[0]):
            for y in range(self.shape[1]):
                chemokine = np.random.randint(0,100)
                self.automaton.grid[(x,y)]['chemokine'] = chemokine
                self.automaton.max_chemokine = max(self.automaton.max_chemokine, chemokine)

        self.model_params['bacteria_threshold_for_macrophage_recruitment'] = 100
        self.model_params['chemokine_scale_for_macrophage_recruitment_below_threshold'] = 1.1
        self.model_params['macrophage_recruitment_probability'] = 100.0
        events = self.automaton.macrophage_recruitment()
        self.assertEqual(len(events), 0)

        self.model_params['bacteria_threshold_for_macrophage_recruitment'] = 0
        self.model_params['chemokine_scale_for_macrophage_recruitment_above_threshold'] = 1.1
        self.model_params['macrophage_recruitment_probability'] = 0.0
        events = self.automaton.macrophage_recruitment()
        self.assertEqual(len(events), 0)

    def test_chemo_kill_bacteria(self):
        # FAST
        self.automaton.model_parameters['chemotherapy_scale_for_kill_fast_bacteria'] = 0.5
        self.automaton.max_chemotherapy = 100.0
        bac_to_kill = self.fb[0]
        self.automaton.grid[bac_to_kill]['chemotherapy'] = 77.0
        bacs_to_save = self.fb[1:]
        for a in bacs_to_save:
            self.automaton.grid[a]['chemotherapy'] = 1.0
        events = self.automaton.chemotherapy_killing_bacteria()
        self.assertEqual(len(events), 1)
        self.assertEqual(events[0].bacterium_address, bac_to_kill)

        # SLOW
        self.automaton.model_parameters['chemotherapy_scale_for_kill_fast_bacteria'] = 1.01
        self.automaton.model_parameters['chemotherapy_scale_for_kill_slow_bacteria'] = 0.5
        self.automaton.max_chemotherapy = 100.0
        bac_to_kill = self.sb[0]
        self.automaton.grid[bac_to_kill]['chemotherapy'] = 77.0
        bacs_to_save = self.sb[1:]
        for a in bacs_to_save:
            self.automaton.grid[a]['chemotherapy'] = 1.0
        events = self.automaton.chemotherapy_killing_bacteria()
        self.assertEqual(len(events), 1)
        self.assertEqual(events[0].bacterium_address, bac_to_kill)

    def test_chemo_kill_bacteria_negative(self):
        # FAST
        self.automaton.model_parameters['chemotherapy_scale_for_kill_fast_bacteria'] = 0.5
        self.automaton.model_parameters['chemotherapy_scale_for_kill_slow_bacteria'] = 0.5
        self.automaton.max_chemotherapy = 100.0
        bacs_to_save = self.fb
        for a in bacs_to_save:
            self.automaton.grid[a]['chemotherapy'] = 1.0
        bacs_to_save = self.sb
        for a in bacs_to_save:
            self.automaton.grid[a]['chemotherapy'] = 1.0
        events = self.automaton.chemotherapy_killing_bacteria()
        self.assertEqual(len(events), 0)

    def test_chemo_kill_macrophages(self):
        self.automaton.model_parameters['chemotherapy_scale_for_kill_macrophage'] = -0.1

        self.automaton.max_chemotherapy = 100.0

        # Infected - killed
        self.automaton.grid[self.macs[0]]['contents'].state = 'infected'
        self.automaton.grid[self.macs[0]]['chemotherapy'] = 70.0
        # Chronically infected - killed
        self.automaton.grid[self.macs[1]]['contents'].state = 'chronically_infected'
        self.automaton.grid[self.macs[1]]['chemotherapy'] = 70.0

        events = self.automaton.chemotherapy_killing_macrophages()
        self.assertEqual(len(events), 2)
        for e in events:
            self.assertTrue(isinstance(e, ChemoKillMacrophage))
        addresses = [e.macrophage_address for e in events]
        self.assertItemsEqual(addresses, [self.macs[0], self.macs[1]])

    def test_chemo_kill_macrophages_negative_state(self):
        self.automaton.model_parameters['chemotherapy_scale_for_kill_macrophage'] = -0.1

        self.automaton.max_chemotherapy = 100.0

        # Resting - not killed
        self.automaton.grid[self.macs[0]]['contents'].state = 'resting'
        self.automaton.grid[self.macs[0]]['chemotherapy'] = 70.0
        # Active - not killed
        self.automaton.grid[self.macs[1]]['contents'].state = 'active'
        self.automaton.grid[self.macs[1]]['chemotherapy'] = 70.0

        events = self.automaton.chemotherapy_killing_macrophages()
        self.assertEqual(len(events), 0)

    def test_chemo_kill_macrophages_negative_scale(self):
        self.automaton.model_parameters['chemotherapy_scale_for_kill_macrophage'] = 0.8

        self.automaton.max_chemotherapy = 100.0

        # Infected
        self.automaton.grid[self.macs[0]]['contents'].state = 'infected'
        self.automaton.grid[self.macs[0]]['chemotherapy'] = 70.0
        # Chronically infected
        self.automaton.grid[self.macs[1]]['contents'].state = 'chronically_infected'
        self.automaton.grid[self.macs[1]]['chemotherapy'] = 70.0

        events = self.automaton.chemotherapy_killing_macrophages()
        self.assertEqual(len(events), 0)

    def test_t_cell_death(self):
        self.automaton.model_parameters['t_cell_age_threshold'] = 1.0
        self.automaton.time_step = 1

        # Add t-cells
        t_cell = TCell((7, 2))
        self.automaton.t_cells.append(t_cell)

        events = self.automaton.t_cell_processes()
        self.assertEqual(len(events), 1)
        self.assertTrue(isinstance(events[0], TCellDeath))
        self.assertEqual(events[0].t_cell_address, (7, 2))

    def test_find_max_chemokine_neighbour(self):

        self.automaton.grid[(0,0)]['chemokine'] = 1.0
        self.automaton.max_chemokine = 1.0

        neighbours = self.automaton.moore_neighbours((1,1), 1)
        address, scale = self.automaton.find_max_chemokine_neighbour(neighbours)
        self.assertEqual(address, (0,0))
        self.assertEqual(scale, 1.0)

        # tie-break
        self.automaton.grid[(2, 0)]['chemokine'] = 1.0
        np.random.seed(101) # Force pick of (2, 0)
        neighbours = self.automaton.moore_neighbours((1, 1), 1)
        address, scale = self.automaton.find_max_chemokine_neighbour(neighbours)
        self.assertEqual(address, (2, 0))
        self.assertEqual(scale, 1.0)

    def test_t_cell_moves_not_random(self):

        self.automaton.model_parameters['t_cell_movement_time'] = 1
        self.automaton.grid[(7,1)]['chemokine'] = 100.0
        self.automaton.max_chemokine = 100.0

        # Add t-cells
        t_cell = TCell((7, 2))
        self.automaton.t_cells.append(t_cell)

        events = self.automaton.t_cell_processes()
        self.assertEqual(len(events), 1)
        self.assertTrue(isinstance(events[0], TCellMovement))
        self.assertEqual(events[0].tcell_from_address, (7, 2))
        self.assertEqual(events[0].tcell_to_address, (7, 1))

    def test_t_cell_moves_random(self):

        self.automaton.model_parameters['t_cell_movement_time'] = 1
        self.automaton.model_parameters['t_cell_random_move_probability'] = 101.0
        self.automaton.grid[(7,1)]['chemokine'] = 100.0
        self.automaton.max_chemokine = 100.0

        # Add t-cells
        t_cell = TCell((5, 1))
        self.automaton.t_cells.append(t_cell)

        np.random.seed(101)

        events = self.automaton.t_cell_processes()
        self.assertEqual(len(events), 1)
        self.assertTrue(isinstance(events[0], TCellMovement))
        self.assertEqual(events[0].tcell_from_address, (5, 1))
        self.assertEqual(events[0].tcell_to_address, (6, 1))

    def test_t_cell_kills_mac(self):

        self.automaton.model_parameters['t_cell_kills_macrophage_probability'] = 100.0
        self.automaton.model_parameters['t_cell_movement_time'] = 1
        self.automaton.grid[(7,1)]['chemokine'] = 100.0
        self.automaton.max_chemokine = 100.0

        mac = Macrophage((7, 1), 'infected')
        self.automaton.macrophages.append(mac)
        self.automaton.grid[(7, 1)]['contents'] = mac

        # Add t-cells
        t_cell = TCell((7, 2))
        self.automaton.t_cells.append(t_cell)

        events = self.automaton.t_cell_processes()
        self.assertEqual(len(events), 1)
        self.assertTrue(isinstance(events[0], TCellKillsMacrophage))
        self.assertEqual(events[0].tcell_address, (7, 2))
        self.assertEqual(events[0].macrophage_address, (7, 1))

        mac.state = 'chronically_infected'
        events = self.automaton.t_cell_processes()
        self.assertEqual(len(events), 1)
        self.assertTrue(isinstance(events[0], TCellKillsMacrophage))
        self.assertEqual(events[0].tcell_address, (7, 2))
        self.assertEqual(events[0].macrophage_address, (7, 1))


    def test_t_cell_kills_mac_negative_prob(self):

        self.automaton.model_parameters['t_cell_kills_macrophage_probability'] = 0.0
        self.automaton.model_parameters['t_cell_movement_time'] = 1
        self.automaton.grid[(7,1)]['chemokine'] = 100.0
        self.automaton.max_chemokine = 100.0

        mac = Macrophage((7, 1), 'infected')
        self.automaton.macrophages.append(mac)
        self.automaton.grid[(7, 1)]['contents'] = mac

        # Add t-cells
        t_cell = TCell((7, 2))
        self.automaton.t_cells.append(t_cell)

        events = self.automaton.t_cell_processes()
        self.assertEqual(len(events), 0)

    def test_t_cell_kills_mac_negative_state(self):

        self.automaton.model_parameters['t_cell_kills_macrophage_probability'] = 100.0
        self.automaton.model_parameters['t_cell_movement_time'] = 1
        self.automaton.grid[(7,1)]['chemokine'] = 100.0
        self.automaton.max_chemokine = 100.0

        mac = Macrophage((7, 1), 'resting')
        self.automaton.macrophages.append(mac)
        self.automaton.grid[(7, 1)]['contents'] = mac

        # Add t-cells
        t_cell = TCell((7, 2))
        self.automaton.t_cells.append(t_cell)

        events = self.automaton.t_cell_processes()
        self.assertEqual(len(events), 0)

        mac.state = 'active'
        events = self.automaton.t_cell_processes()
        self.assertEqual(len(events), 0)

    def test_macrophage_resting_death(self):
        self.automaton.model_parameters['resting_macrophage_age_limit'] = 1.0

        self.automaton.macrophages = []
        for m in self.macs:
            self.automaton.grid[m]['contents'] = 0.0

        mac = Macrophage((8,8),'resting')
        self.automaton.macrophages.append(mac)
        self.automaton.grid[(8,8)]['contents'] = mac

        self.automaton.time = 100.0
        events = self.automaton.macrophage_processes()

        self.assertEqual(len(events), 1)
        self.assertTrue(isinstance(events[0], MacrophageDeath))
        self.assertEqual(events[0].macrophage_address, (8,8))

    def test_macrophage_resting_move_not_random(self):
        self.automaton.model_parameters['resting_macrophage_movement_time'] = 1.0
        self.automaton.model_parameters['minimum_chemokine_for_resting_macrophage_movement'] = 0.0

        self.automaton.grid[(7,8)]['chemokine'] = 100.0
        self.automaton.max_chemokine = 100.0

        self.automaton.macrophages = []
        for m in self.macs:
            self.automaton.grid[m]['contents'] = 0.0

        mac = Macrophage((8, 8), 'resting')
        self.automaton.macrophages.append(mac)
        self.automaton.grid[(8, 8)]['contents'] = mac

        self.automaton.time = 100.0
        events = self.automaton.macrophage_processes()

        self.assertEqual(len(events), 1)
        self.assertTrue(isinstance(events[0], MacrophageMovement))
        self.assertEqual(events[0].macrophage_from_address, (8, 8))
        self.assertEqual(events[0].macrophage_to_address, (7, 8))

    def test_macrophage_resting_move_random_through_prob(self):
        self.automaton.model_parameters['resting_macrophage_movement_time'] = 1.0
        self.automaton.model_parameters['prob_resting_macrophage_random_move'] = 100.0

        self.automaton.grid[(7,8)]['chemokine'] = 100.0
        self.automaton.max_chemokine = 100.0

        self.automaton.macrophages = []
        for m in self.macs:
            self.automaton.grid[m]['contents'] = 0.0

        mac = Macrophage((8, 8), 'resting')
        self.automaton.macrophages.append(mac)
        self.automaton.grid[(8, 8)]['contents'] = mac

        self.automaton.time = 100.0
        np.random.seed(101)
        events = self.automaton.macrophage_processes()

        self.assertEqual(len(events), 1)
        self.assertTrue(isinstance(events[0], MacrophageMovement))
        self.assertEqual(events[0].macrophage_from_address, (8, 8))
        self.assertEqual(events[0].macrophage_to_address, (7, 7))

    def test_macrophage_resting_move_random_through_not_enough_chemokine(self):
        self.automaton.model_parameters['resting_macrophage_movement_time'] = 1.0
        self.automaton.model_parameters['minimum_chemokine_for_resting_macrophage_movement'] = 101.0

        self.automaton.grid[(7, 8)]['chemokine'] = 100.0
        self.automaton.max_chemokine = 100.0

        self.automaton.macrophages = []
        for m in self.macs:
            self.automaton.grid[m]['contents'] = 0.0

        mac = Macrophage((8, 8), 'resting')
        self.automaton.macrophages.append(mac)
        self.automaton.grid[(8, 8)]['contents'] = mac

        self.automaton.time = 100.0
        np.random.seed(101)
        events = self.automaton.macrophage_processes()

        self.assertEqual(len(events), 1)
        self.assertTrue(isinstance(events[0], MacrophageMovement))
        self.assertEqual(events[0].macrophage_from_address, (8, 8))
        self.assertEqual(events[0].macrophage_to_address, (7, 7))

    def test_macrophage_resting_kill_bac(self):
        self.automaton.model_parameters['resting_macrophage_movement_time'] = 1.0
        self.automaton.model_parameters['minimum_chemokine_for_resting_macrophage_movement'] = 0.0

        self.automaton.grid[(7,8)]['chemokine'] = 100.0
        self.automaton.max_chemokine = 100.0

        self.automaton.macrophages = []
        for m in self.macs:
            self.automaton.grid[m]['contents'] = 0.0

        mac = Macrophage((8, 8), 'resting')
        self.automaton.macrophages.append(mac)
        self.automaton.grid[(8, 8)]['contents'] = mac

        self.automaton.time = 100.0
        events = self.automaton.macrophage_processes()

        self.assertEqual(len(events), 1)
        self.assertTrue(isinstance(events[0], MacrophageMovement))
        self.assertEqual(events[0].macrophage_from_address, (8, 8))
        self.assertEqual(events[0].macrophage_to_address, (7, 8))

    def test_macrophage_does_nothing(self):

        self.automaton.macrophages = []
        for m in self.macs:
            self.automaton.grid[m]['contents'] = 0.0

        mac = Macrophage((8,8),'resting')
        self.automaton.macrophages.append(mac)
        self.automaton.grid[(8,8)]['contents'] = mac

        self.automaton.time = 100.0
        events = self.automaton.macrophage_processes()

        self.assertEqual(len(events), 0)

    def test_macrophage_active_death(self):
        self.automaton.model_parameters['active_macrophage_age_limit'] = 0.0

        self.automaton.macrophages = []
        for m in self.macs:
            self.automaton.grid[m]['contents'] = 0.0

        mac = Macrophage((8,8),'active')
        self.automaton.macrophages.append(mac)
        self.automaton.grid[(8,8)]['contents'] = mac

        self.automaton.time = 100.0
        events = self.automaton.macrophage_processes()

        self.assertEqual(len(events), 1)
        self.assertTrue(isinstance(events[0], MacrophageDeath))
        self.assertEqual(events[0].macrophage_address, (8,8))

    def test_macrophage_active_move(self):
        self.automaton.model_parameters['active_macrophage_movement_time'] = 10.0
        self.automaton.time = 10.0
        self.automaton.macrophages = []
        for m in self.macs:
            self.automaton.grid[m]['contents'] = 0.0

        self.automaton.grid[(9,9)]['chemokine'] = 100.0
        self.automaton.max_chemokine = 100.0

        mac = Macrophage((8, 8), 'active')
        self.automaton.macrophages.append(mac)
        self.automaton.grid[(8, 8)]['contents'] = mac

        self.automaton.time = 100.0
        events = self.automaton.macrophage_processes()

        self.assertEqual(len(events), 1)
        self.assertTrue(isinstance(events[0], MacrophageMovement))
        self.assertEqual(events[0].macrophage_from_address, (8, 8))
        self.assertEqual(events[0].macrophage_to_address, (9, 9))

    def test_macrophage_kills_bacterium(self):
        self.automaton.model_parameters['active_macrophage_movement_time'] = 10.0
        self.automaton.model_parameters['prob_active_macrophage_kill_fast_bacteria'] = 100.0
        self.automaton.model_parameters['prob_active_macrophage_kill_slow_bacteria'] = 100.0
        self.automaton.time = 10.0

        self.automaton.macrophages = []
        for m in self.macs:
            self.automaton.grid[m]['contents'] = 0.0

        self.automaton.grid[(9, 9)]['chemokine'] = 100.0
        self.automaton.max_chemokine = 100.0

        bac = Bacterium((9, 9), 'fast')
        self.automaton.bacteria.append(bac)
        self.automaton.grid[(9, 9)]['contents'] = bac

        mac = Macrophage((8, 8), 'active')
        self.automaton.macrophages.append(mac)
        self.automaton.grid[(8, 8)]['contents'] = mac

        events = self.automaton.macrophage_processes()

        self.assertEqual(len(events), 1)
        self.assertTrue(isinstance(events[0], MacrophageIngestsBacterium))
        self.assertEqual(events[0].macrophage_address, (8, 8))
        self.assertEqual(events[0].bacterium_address, (9, 9))

        bac.metabolism = 'slow'
        events = self.automaton.macrophage_processes()

        self.assertEqual(len(events), 1)
        self.assertTrue(isinstance(events[0], MacrophageIngestsBacterium))
        self.assertEqual(events[0].macrophage_address, (8, 8))
        self.assertEqual(events[0].bacterium_address, (9, 9))

    def test_macrophage_kills_bacterium_negative(self):
        self.automaton.model_parameters['active_macrophage_movement_time'] = 10.0
        self.automaton.model_parameters['prob_active_macrophage_kill_fast_bacteria'] = 0.0
        self.automaton.model_parameters['prob_active_macrophage_kill_slow_bacteria'] = 0.0
        self.automaton.time = 10.0

        self.automaton.macrophages = []
        for m in self.macs:
            self.automaton.grid[m]['contents'] = 0.0

        self.automaton.grid[(9, 9)]['chemokine'] = 100.0
        self.automaton.max_chemokine = 100.0

        bac = Bacterium((9, 9), 'fast')
        self.automaton.bacteria.append(bac)
        self.automaton.grid[(9, 9)]['contents'] = bac

        mac = Macrophage((8, 8), 'active')
        self.automaton.macrophages.append(mac)
        self.automaton.grid[(8, 8)]['contents'] = mac

        events = self.automaton.macrophage_processes()
        self.assertEqual(len(events), 0)

        bac.metabolism = 'slow'
        events = self.automaton.macrophage_processes()
        self.assertEqual(len(events), 0)

    def test_macrophage_infected_death(self):
        self.automaton.model_parameters['infected_macrophage_age_limit'] = 1.0

        self.automaton.macrophages = []
        for m in self.macs:
            self.automaton.grid[m]['contents'] = 0.0

        mac = Macrophage((8,8),'infected')
        self.automaton.macrophages.append(mac)
        self.automaton.grid[(8,8)]['contents'] = mac

        self.automaton.time = 100.0
        events = self.automaton.macrophage_processes()

        self.assertEqual(len(events), 1)
        self.assertTrue(isinstance(events[0], MacrophageDeath))
        self.assertEqual(events[0].macrophage_address, (8,8))

    def test_macrophage_infected_move(self):
        self.automaton.model_parameters['infected_macrophage_movement_time'] = 10.0
        self.automaton.time = 10.0

        self.automaton.macrophages = []
        for m in self.macs:
            self.automaton.grid[m]['contents'] = 0.0

        mac = Macrophage((8,8),'infected')
        self.automaton.macrophages.append(mac)
        self.automaton.grid[(8,8)]['contents'] = mac

        self.automaton.grid[(9, 9)]['chemokine'] = 100.0
        self.automaton.max_chemokine = 100.0

        events = self.automaton.macrophage_processes()

        self.assertEqual(len(events), 1)
        self.assertTrue(isinstance(events[0], MacrophageMovement))
        self.assertEqual(events[0].macrophage_from_address, (8, 8))
        self.assertEqual(events[0].macrophage_to_address, (9, 9))

    def test_macrophage_infected_kill_bacteria(self):
        self.automaton.model_parameters['infected_macrophage_movement_time'] = 10.0
        self.automaton.time = 10.0

        self.automaton.macrophages = []
        for m in self.macs:
            self.automaton.grid[m]['contents'] = 0.0

        mac = Macrophage((8,8),'infected')
        self.automaton.macrophages.append(mac)
        self.automaton.grid[(8,8)]['contents'] = mac

        self.automaton.grid[(9, 9)]['chemokine'] = 100.0
        self.automaton.max_chemokine = 100.0

        bac = Bacterium((9, 9), 'fast')
        self.automaton.bacteria.append(bac)
        self.automaton.grid[(9, 9)]['contents'] = bac

        events = self.automaton.macrophage_processes()

        self.assertEqual(len(events), 1)
        self.assertTrue(isinstance(events[0], MacrophageIngestsBacterium))
        self.assertEqual(events[0].macrophage_address, (8, 8))
        self.assertEqual(events[0].bacterium_address, (9, 9))

    def test_macrophage_chr_infected_death(self):
        self.automaton.model_parameters['chronically_infected_macrophage_age_limit'] = 1.0

        self.automaton.macrophages = []
        for m in self.macs:
            self.automaton.grid[m]['contents'] = 0.0

        mac = Macrophage((8,8),'chronically_infected')
        self.automaton.macrophages.append(mac)
        self.automaton.grid[(8,8)]['contents'] = mac

        self.automaton.time = 100.0
        events = self.automaton.macrophage_processes()

        self.assertEqual(len(events), 1)
        self.assertTrue(isinstance(events[0], MacrophageDeath))
        self.assertEqual(events[0].macrophage_address, (8,8))

    def test_macrophage_chr_infected_move(self):
        self.automaton.model_parameters['chronically_infected_macrophage_movement_time'] = 10.0
        self.automaton.time = 10.0

        self.automaton.macrophages = []
        for m in self.macs:
            self.automaton.grid[m]['contents'] = 0.0

        mac = Macrophage((8,8),'chronically_infected')
        self.automaton.macrophages.append(mac)
        self.automaton.grid[(8,8)]['contents'] = mac

        self.automaton.grid[(9, 9)]['chemokine'] = 100.0
        self.automaton.max_chemokine = 100.0

        events = self.automaton.macrophage_processes()

        self.assertEqual(len(events), 1)
        self.assertTrue(isinstance(events[0], MacrophageMovement))
        self.assertEqual(events[0].macrophage_from_address, (8, 8))
        self.assertEqual(events[0].macrophage_to_address, (9, 9))

    def test_macrophage_chr_infected_kill_bacteria(self):
        self.automaton.model_parameters['chronically_infected_macrophage_movement_time'] = 10.0
        self.automaton.time = 10.0

        self.automaton.macrophages = []
        for m in self.macs:
            self.automaton.grid[m]['contents'] = 0.0

        mac = Macrophage((8,8),'chronically_infected')
        self.automaton.macrophages.append(mac)
        self.automaton.grid[(8,8)]['contents'] = mac

        self.automaton.grid[(9, 9)]['chemokine'] = 100.0
        self.automaton.max_chemokine = 100.0

        bac = Bacterium((9, 9), 'fast')
        self.automaton.bacteria.append(bac)
        self.automaton.grid[(9, 9)]['contents'] = bac

        events = self.automaton.macrophage_processes()

        self.assertEqual(len(events), 1)
        self.assertTrue(isinstance(events[0], MacrophageIngestsBacterium))
        self.assertEqual(events[0].macrophage_address, (8, 8))
        self.assertEqual(events[0].bacterium_address, (9, 9))

    def test_macrophage_activation(self):
        self.automaton.model_parameters['chemokine_scale_for_macrophage_activation'] = 0.0
        self.automaton.grid[self.macs[0]]['chemokine'] = 100.0
        self.automaton.max_chemokine = 100.0
        self.automaton.time = 88
        events = self.automaton.macrophage_processes()
        self.assertEqual(len(events), 1)
        self.assertTrue(isinstance(events[0], MacrophageActivation))
        self.assertEqual(events[0].macrophage_address, self.macs[0])
        self.assertEqual(events[0].new_state, 'active')

    def test_macrophage_deactivation(self):
        self.automaton.grid[self.macs[0]]['contents'].state = 'active'
        self.automaton.model_parameters['chemokine_scale_for_macrophage_deactivation'] = 101.0
        self.automaton.grid[self.macs[0]]['chemokine'] = 0.0
        self.automaton.max_chemokine = 100.0
        events = self.automaton.macrophage_processes()
        self.assertEqual(len(events), 1)
        self.assertTrue(isinstance(events[0], MacrophageActivation))
        self.assertEqual(events[0].macrophage_address, self.macs[0])
        self.assertEqual(events[0].new_state, 'resting')


    def test_macrophage_burst(self):
        self.automaton.model_parameters['bacteria_to_burst_macrophage'] = 20
        self.automaton.grid[self.macs[0]]['contents'].state = 'chronically_infected'
        self.automaton.grid[self.macs[0]]['contents'].intracellular_bacteria = 20

        np.random.seed(101)
        events = self.automaton.macrophage_processes()
        self.assertEqual(len(events), 1)
        self.assertTrue(isinstance(events[0], MacrophageBursts))
        self.assertEqual(events[0].macrophage_address, self.macs[0])
        self.assertItemsEqual(events[0].new_bacteria_addresses,
                    [(8, 9), (9, 8), (9, 7), (7, 8), (7, 9), (8, 7), (6, 8), (8, 6), (6, 7), (9, 6), (7, 6), (6, 9)])

    def test_bacteria_fast_to_slow(self):

        self.automaton.bacteria = []
        for b in self.fb:
            self.automaton.grid[b]['contents'] = 0.0
        for b in self.sb:
            self.automaton.grid[b]['contents'] = 0.0

        bac = Bacterium((8,8), 'fast')
        self.automaton.bacteria.append(bac)
        self.automaton.grid[(8,8)]['contents'] = bac

        self.model_params['oxygen_scale_for_metabolism_change_to_slow'] = 1.01
        self.automaton.time = 999.0
        events = self.automaton.bacteria_state_changes()
        self.assertEqual(len(events), 1)
        self.assertTrue(isinstance(events[0], BacteriumStateChange))
        self.assertEqual(events[0].attribute, "metabolism")
        self.assertEqual(events[0].value, "slow")

    def test_bacteria_fast_to_slow_negative(self):

        self.automaton.bacteria = []
        for b in self.fb:
            self.automaton.grid[b]['contents'] = 0.0
        for b in self.sb:
            self.automaton.grid[b]['contents'] = 0.0

        bac = Bacterium((8,8), 'fast')
        self.automaton.bacteria.append(bac)
        self.automaton.grid[(8,8)]['contents'] = bac

        self.model_params['oxygen_scale_for_metabolism_change_to_slow'] = 0.5

        self.automaton.max_oxygen = 100.0
        self.automaton.grid[(8,8)]['oxygen'] = 75.0

        self.automaton.time = 999.0
        events = self.automaton.bacteria_state_changes()
        self.assertEqual(len(events), 0)


    def test_bacteria_slow_to_fast(self):
        self.automaton.bacteria = []
        for b in self.fb:
            self.automaton.grid[b]['contents'] = 0.0
        for b in self.sb:
            self.automaton.grid[b]['contents'] = 0.0

        bac = Bacterium((8,8), 'slow')
        self.automaton.bacteria.append(bac)
        self.automaton.grid[(8,8)]['contents'] = bac

        self.model_params['oxygen_scale_for_metabolism_change_to_fast'] = -1
        self.automaton.time = 999.0
        events = self.automaton.bacteria_state_changes()
        self.assertEqual(len(events), 1)
        self.assertTrue(isinstance(events[0], BacteriumStateChange))
        self.assertEqual(events[0].attribute, "metabolism")
        self.assertEqual(events[0].value, "fast")

    def test_bacteria_slow_to_fast_negative(self):
        self.automaton.bacteria = []
        for b in self.fb:
            self.automaton.grid[b]['contents'] = 0.0
        for b in self.sb:
            self.automaton.grid[b]['contents'] = 0.0

        bac = Bacterium((8,8), 'slow')
        self.automaton.bacteria.append(bac)
        self.automaton.grid[(8,8)]['contents'] = bac

        self.model_params['oxygen_scale_for_metabolism_change_to_fast'] = 0.5

        self.automaton.max_oxygen = 100.0
        self.automaton.grid[(8, 8)]['oxygen'] = 25.0

        self.automaton.time = 999.0
        events = self.automaton.bacteria_state_changes()
        self.assertEqual(len(events), 0)

    def test_bacteria_resting_to_non_resting(self):
        self.automaton.bacteria = []
        for b in self.fb:
            self.automaton.grid[b]['contents'] = 0.0
        for b in self.sb:
            self.automaton.grid[b]['contents'] = 0.0

        bac = Bacterium((8, 8), 'fast')
        self.automaton.bacteria.append(bac)
        self.automaton.grid[(8, 8)]['contents'] = bac
        bac.resting = True
        events = self.automaton.bacteria_state_changes()
        self.assertEqual(len(events), 1)
        self.assertTrue(isinstance(events[0], BacteriumStateChange))
        self.assertEqual(events[0].attribute, "resting")
        self.assertEqual(events[0].value, False)

    def test_bacteria_resting_to_non_resting_negative(self):
        self.automaton.bacteria = []
        for b in self.fb:
            self.automaton.grid[b]['contents'] = 0.0
        for b in self.sb:
            self.automaton.grid[b]['contents'] = 0.0

        for x in range(self.shape[0]):
            for y in range(self.shape[1]):
                if not (x == 8 and y == 8):
                    self.automaton.grid[(x,y)]['contents'] = Caseum((x,y))

        bac = Bacterium((8, 8), 'fast')
        self.automaton.bacteria.append(bac)
        self.automaton.grid[(8, 8)]['contents'] = bac
        bac.resting = True
        events = self.automaton.bacteria_state_changes()
        self.assertEqual(len(events), 0)
        

if __name__ == '__main__':
    unittest.main()
