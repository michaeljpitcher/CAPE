import unittest
from CAPE.TBAutomaton import *
import shutil
import os


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

        self.bv = [(1, 1), (2, 2), (3, 3)]
        self.macs = [(9, 9), (8, 8), (7, 7)]
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
        self.assertEqual(len(self.automaton.tcells), 0.0)
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
        self.automaton.tcells.append(tcell)
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
        self.assertEqual(self.automaton.oxygen_scale((0,0)), 0.0)
        self.assertEqual(self.automaton.chemotherapy_scale((0,0)), 0.0)
        self.assertEqual(self.automaton.chemokine_scale((0,0)), 0.0)

        self.automaton.max_oxygen = 100.0
        self.automaton.max_chemotherapy = 75.0
        self.automaton.max_chemokine = 40.0

        self.automaton.grid[(0, 0)]['oxygen'] = 56.0
        self.automaton.grid[(4, 7)]['chemotherapy'] = 36.0
        self.automaton.grid[(8, 8)]['chemokine'] = 1.0

        self.assertEqual(self.automaton.oxygen_scale((0, 0)), 0.56)
        self.assertAlmostEqual(self.automaton.chemotherapy_scale((4, 7)), 0.48)
        self.assertEqual(self.automaton.chemokine_scale((8, 8)), 0.025)


if __name__ == '__main__':
    unittest.main()
