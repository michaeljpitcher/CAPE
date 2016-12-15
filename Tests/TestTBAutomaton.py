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

        self.bv = [(1, 1), (2, 2), (3, 3)]
        self.macs = [(9, 9), (8, 8), (7, 7)]
        self.fb = [(8, 1), (8, 2), (8, 3)]
        self.sb = [(1, 7), (2, 7), (3, 7)]

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

        cas = Caseum()
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

        cas = Caseum()
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

if __name__ == '__main__':
    unittest.main()
