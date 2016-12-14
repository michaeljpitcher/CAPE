import unittest
import os
from CAPE.TBAutomaton import *
import shutil


class DiffusionTestCase(unittest.TestCase):

    def setUp(self):
        self.shape = (10, 10)
        self.time_params = {}
        self.time_params['initial_time'] = 0.0
        self.time_params['time_step'] = 0.001
        self.time_params['time_limit'] = 1.0
        self.model_params = {}
        self.model_params['chemotherapy_schedule1_start_lower'] = 10000.0
        self.model_params['chemotherapy_schedule1_start_upper'] = 10001.0
        self.model_params['blood_vessel_value'] = 1.0
        self.model_params['initial_oxygen'] = 1.0
        self.model_params['oxygen_diffusion'] = 1.0
        self.model_params['chemotherapy_diffusion'] = 0.0
        self.model_params['spatial_step'] = 0.2
        self.model_params['oxygen_from_source'] = 0.0
        self.model_params['oxygen_uptake_from_bacteria'] = 0.0
        self.model_params['chemokine_diffusion'] = 0.0
        self.model_params['chemokine_from_bacteria'] = 0.0
        self.model_params['chemokine_from_macrophage'] = 0.0
        self.model_params['chemokine_decay'] = 0.0

        self.bv = [(4,4)]
        self.macs = []
        self.fb = []
        self.sb = []

        self.output_loc = 'test_output'
        if not os.path.exists(self.output_loc):
            os.makedirs(self.output_loc)
        self.automaton = TBAutomaton(self.shape, self.time_params, self.model_params, self.output_loc,
                                     self.bv, self.macs, self.fb, self.sb)
    def tearDown(self):
        # Close output files and delete
        self.automaton.close_files()
        shutil.rmtree(self.output_loc)

    def test_initial(self):
        for x in range(10):
            for y in range(10):
                if (x,y) in self.bv:
                    self.assertEqual(self.automaton.grid[(x,y)]['oxygen'],
                         self.model_params['initial_oxygen'] * self.model_params['blood_vessel_value'])
                else:
                    self.assertEqual(self.automaton.grid[(x, y)]['oxygen'], 0.0)
                self.assertEqual(self.automaton.grid[(x, y)]['chemotherapy'], 0.0)
                self.assertEqual(self.automaton.grid[(x, y)]['chemokine'], 0.0)

    def test_oxygen_diffusion(self):
        self.automaton.diffusion(False)
        previous_oxygen_at_source_cell = self.model_params['initial_oxygen'] * self.model_params['blood_vessel_value']
        expected_oxygen_at_cell = (previous_oxygen_at_source_cell) + self.time_params['time_step'] * (
            (((self.model_params['oxygen_diffusion'] + self.model_params['oxygen_diffusion']) / 2) *
             (0 - previous_oxygen_at_source_cell) -
            ((self.model_params['oxygen_diffusion'] + self.model_params['oxygen_diffusion']) / 2) *
             (previous_oxygen_at_source_cell - 0)) / self.model_params['spatial_step']**2
            + (((self.model_params['oxygen_diffusion'] + self.model_params['oxygen_diffusion']) / 2) *
               (0 - previous_oxygen_at_source_cell) -
            ((self.model_params['oxygen_diffusion'] + self.model_params['oxygen_diffusion']) / 2) *
               (previous_oxygen_at_source_cell - 0)) / self.model_params['spatial_step']**2
            + 0 + 0)
        self.assertEqual(self.automaton.work_grid[(4,4)]['oxygen'], expected_oxygen_at_cell)
        # Neighbours - check the neighbours get the oxygen
        previous_oxygen_at_neighb_cell = 0
        expected_oxygen_at_cell = (0) + self.time_params['time_step'] * (
            (((self.model_params['oxygen_diffusion'] + self.model_params['oxygen_diffusion']) / 2) *
             (previous_oxygen_at_source_cell - 0) -
             ((self.model_params['oxygen_diffusion'] + self.model_params['oxygen_diffusion']) / 2) *
             (0 - 0)) / self.model_params['spatial_step'] ** 2
            + (((self.model_params['oxygen_diffusion'] + self.model_params['oxygen_diffusion']) / 2) *
               (0 - 0) -
               ((self.model_params['oxygen_diffusion'] + self.model_params['oxygen_diffusion']) / 2) *
               (0 - 0)) / self.model_params['spatial_step'] ** 2
            + 0 + 0)
        self.assertEqual(self.automaton.work_grid[(3, 4)]['oxygen'], expected_oxygen_at_cell)
        self.assertEqual(self.automaton.work_grid[(4, 3)]['oxygen'], expected_oxygen_at_cell)
        self.assertEqual(self.automaton.work_grid[(4, 5)]['oxygen'], expected_oxygen_at_cell)
        self.assertEqual(self.automaton.work_grid[(5, 4)]['oxygen'], expected_oxygen_at_cell)

    def test_oxygen_from_source_no_diffusion(self):
        previous_oxygen_at_cell = self.model_params['initial_oxygen'] * self.model_params['blood_vessel_value']
        # Set from source
        self.automaton.model_parameters['oxygen_from_source'] = 1.0
        # Turn off diffusion
        self.automaton.grid['oxygen_diffusion_rate'] = np.zeros(self.shape,dtype=float)

        self.automaton.diffusion(False)

        expected_oxygen_at_cell = (previous_oxygen_at_cell) + self.time_params['time_step'] * (
            + self.automaton.model_parameters['oxygen_from_source'] *
            self.automaton.model_parameters['blood_vessel_value']
            + 0)
        self.assertEqual(self.automaton.work_grid[(4,4)]['oxygen'], expected_oxygen_at_cell)
        self.assertEqual(self.automaton.work_grid[(3, 4)]['oxygen'], 0.0)
        self.assertEqual(self.automaton.work_grid[(4, 3)]['oxygen'], 0.0)
        self.assertEqual(self.automaton.work_grid[(4, 5)]['oxygen'], 0.0)
        self.assertEqual(self.automaton.work_grid[(5, 4)]['oxygen'], 0.0)

    def test_oxygen_to_bacteria_no_diffusion(self):
        previous_oxygen_at_cell = self.model_params['initial_oxygen'] * self.model_params['blood_vessel_value']
        # Set from source
        self.automaton.model_parameters['oxygen_uptake_from_bacteria'] = 1.0
        # Turn off diffusion
        self.automaton.grid['oxygen_diffusion_rate'] = np.zeros(self.shape,dtype=float)
        b = Bacterium('fast')
        self.automaton.grid[(4,4)]['contents'] = b

        self.automaton.diffusion(False)

        expected_oxygen_at_cell = (previous_oxygen_at_cell) + self.time_params['time_step'] * (
            - self.automaton.model_parameters['oxygen_uptake_from_bacteria'] * previous_oxygen_at_cell
            + 0)
        self.assertEqual(self.automaton.work_grid[(4,4)]['oxygen'], expected_oxygen_at_cell)
        self.assertEqual(self.automaton.work_grid[(3, 4)]['oxygen'], 0.0)
        self.assertEqual(self.automaton.work_grid[(4, 3)]['oxygen'], 0.0)
        self.assertEqual(self.automaton.work_grid[(4, 5)]['oxygen'], 0.0)
        self.assertEqual(self.automaton.work_grid[(5, 4)]['oxygen'], 0.0)

    def test_oxygen_diffusion_differing_rates(self):

        self.automaton.grid['oxygen_diffusion_rate'][(4, 4)] = 1.0
        self.automaton.grid['oxygen_diffusion_rate'][(3, 4)] = 0.7
        self.automaton.grid['oxygen_diffusion_rate'][(4, 3)] = 0.6
        self.automaton.grid['oxygen_diffusion_rate'][(4, 5)] = 0.5
        self.automaton.grid['oxygen_diffusion_rate'][(5, 4)] = 0.4

        self.automaton.diffusion(False)
        previous_oxygen_at_source_cell = self.model_params['initial_oxygen'] * self.model_params['blood_vessel_value']
        expected_oxygen_at_cell = (previous_oxygen_at_source_cell) + self.time_params['time_step'] * (
            (((self.automaton.grid['oxygen_diffusion_rate'][(4, 4)] + self.automaton.grid['oxygen_diffusion_rate'][(3, 4)]) / 2) *
             (0 - previous_oxygen_at_source_cell) -
            ((self.automaton.grid['oxygen_diffusion_rate'][(4, 4)] + self.automaton.grid['oxygen_diffusion_rate'][(5, 4)]) / 2) *
             (previous_oxygen_at_source_cell - 0)) / self.model_params['spatial_step']**2
            + (((self.automaton.grid['oxygen_diffusion_rate'][(4, 4)] + self.automaton.grid['oxygen_diffusion_rate'][(4, 3)]) / 2) *
               (0 - previous_oxygen_at_source_cell) -
            ((self.automaton.grid['oxygen_diffusion_rate'][(4, 4)] + self.automaton.grid['oxygen_diffusion_rate'][(4, 5)]) / 2) *
               (previous_oxygen_at_source_cell - 0)) / self.model_params['spatial_step']**2
            + 0 + 0)
        self.assertEqual(self.automaton.work_grid[(4,4)]['oxygen'], expected_oxygen_at_cell)
        # Neighbours - check the neighbours get the oxygen
        expected_oxygen_at_cell = (0) + self.time_params['time_step'] * (
            (((self.automaton.grid['oxygen_diffusion_rate'][(3, 4)] + self.automaton.grid['oxygen_diffusion_rate'][(4, 4)]) / 2) *
             (previous_oxygen_at_source_cell - 0)) / self.model_params['spatial_step'] ** 2
            + 0 + 0)
        self.assertEqual(self.automaton.work_grid[(3, 4)]['oxygen'], expected_oxygen_at_cell)
        expected_oxygen_at_cell = (0) + self.time_params['time_step'] * (
            (((self.automaton.grid['oxygen_diffusion_rate'][(5, 4)] + self.automaton.grid['oxygen_diffusion_rate'][
                (4, 4)]) / 2) *
             (previous_oxygen_at_source_cell - 0)) / self.model_params['spatial_step'] ** 2
            + 0 + 0)
        self.assertEqual(self.automaton.work_grid[(5, 4)]['oxygen'], expected_oxygen_at_cell)
        expected_oxygen_at_cell = (0) + self.time_params['time_step'] * (
            (((self.automaton.grid['oxygen_diffusion_rate'][(4, 3)] + self.automaton.grid['oxygen_diffusion_rate'][
                (4, 4)]) / 2) *
             (previous_oxygen_at_source_cell - 0)) / self.model_params['spatial_step'] ** 2
            + 0 + 0)
        self.assertEqual(self.automaton.work_grid[(4, 3)]['oxygen'], expected_oxygen_at_cell)
        expected_oxygen_at_cell = (0) + self.time_params['time_step'] * (
            (((self.automaton.grid['oxygen_diffusion_rate'][(4, 5)] + self.automaton.grid['oxygen_diffusion_rate'][
                (4, 4)]) / 2) *
             (previous_oxygen_at_source_cell - 0)) / self.model_params['spatial_step'] ** 2
            + 0 + 0)
        self.assertEqual(self.automaton.work_grid[(4, 5)]['oxygen'], expected_oxygen_at_cell)


if __name__ == '__main__':
    unittest.main()
