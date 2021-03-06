import shutil
import unittest

from TBAutomaton.TBAutomaton import *


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
        self.model_params['chemotherapy_diffusion'] = 1.0
        self.model_params['spatial_step'] = 0.2
        self.model_params['oxygen_from_source'] = 0.0
        self.model_params['oxygen_uptake_from_bacteria'] = 0.0
        self.model_params['chemokine_diffusion'] = 1.0
        self.model_params['chemokine_from_bacteria'] = 0.0
        self.model_params['chemokine_from_macrophage'] = 0.0
        self.model_params['chemokine_decay'] = 0.0
        self.model_params['chemotherapy_from_source'] = 0.0
        self.model_params['chemotherapy_decay'] = 0.0

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
        b = Bacterium((4,4), 'fast')
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

    def test_chemotherapy_diffusion(self):

        self.automaton.grid[(4,4)]['chemotherapy'] = 10.0
        self.automaton.diffusion(True)

        previous_chemo_at_source_cell = 10.0
        expected_chemo_at_cell = (previous_chemo_at_source_cell) + self.time_params['time_step'] * (
            (((self.model_params['chemotherapy_diffusion'] + self.model_params['chemotherapy_diffusion']) / 2) *
             (0 - previous_chemo_at_source_cell) -
             ((self.model_params['chemotherapy_diffusion'] + self.model_params['chemotherapy_diffusion']) / 2) *
             (previous_chemo_at_source_cell - 0)) / self.model_params['spatial_step'] ** 2
            + (((self.model_params['chemotherapy_diffusion'] + self.model_params['chemotherapy_diffusion']) / 2) *
               (0 - previous_chemo_at_source_cell) -
               ((self.model_params['chemotherapy_diffusion'] + self.model_params['chemotherapy_diffusion']) / 2) *
               (previous_chemo_at_source_cell - 0)) / self.model_params['spatial_step'] ** 2
            + 0 + 0)
        self.assertEqual(self.automaton.work_grid[(4, 4)]['chemotherapy'], expected_chemo_at_cell)
        # Neighbours - check the neighbours get the chemotherapy
        previous_chemotherapy_at_neighb_cell = 0
        expected_chemo_at_cell = (0) + self.time_params['time_step'] * (
            (((self.model_params['chemotherapy_diffusion'] + self.model_params['chemotherapy_diffusion']) / 2) *
             (previous_chemo_at_source_cell - 0) -
             ((self.model_params['chemotherapy_diffusion'] + self.model_params['chemotherapy_diffusion']) / 2) *
             (0 - 0)) / self.model_params['spatial_step'] ** 2
            + (((self.model_params['chemotherapy_diffusion'] + self.model_params['chemotherapy_diffusion']) / 2) *
               (0 - 0) -
               ((self.model_params['chemotherapy_diffusion'] + self.model_params['chemotherapy_diffusion']) / 2) *
               (0 - 0)) / self.model_params['spatial_step'] ** 2
            + 0 + 0)
        self.assertEqual(self.automaton.work_grid[(3, 4)]['chemotherapy'], expected_chemo_at_cell)
        self.assertEqual(self.automaton.work_grid[(4, 3)]['chemotherapy'], expected_chemo_at_cell)
        self.assertEqual(self.automaton.work_grid[(4, 5)]['chemotherapy'], expected_chemo_at_cell)
        self.assertEqual(self.automaton.work_grid[(5, 4)]['chemotherapy'], expected_chemo_at_cell)

    def test_chemo_from_source_no_diffusion(self):
        previous_chemotherapy_at_cell = 0.0
        self.automaton.grid[(4, 4)]['chemotherapy'] = 0.0
        # Set from source
        self.automaton.model_parameters['chemotherapy_from_source'] = 1.0
        # Turn off diffusion
        self.automaton.grid['chemotherapy_diffusion_rate'] = np.zeros(self.shape,dtype=float)

        self.automaton.diffusion(True)

        expected_chemotherapy_at_cell = (previous_chemotherapy_at_cell) + self.time_params['time_step'] * (
            + self.automaton.model_parameters['chemotherapy_from_source'] *
            self.automaton.model_parameters['blood_vessel_value']
            + 0)
        self.assertEqual(self.automaton.work_grid[(4, 4)]['chemotherapy'], expected_chemotherapy_at_cell)
        self.assertEqual(self.automaton.work_grid[(3, 4)]['chemotherapy'], 0.0)
        self.assertEqual(self.automaton.work_grid[(4, 3)]['chemotherapy'], 0.0)
        self.assertEqual(self.automaton.work_grid[(4, 5)]['chemotherapy'], 0.0)
        self.assertEqual(self.automaton.work_grid[(5, 4)]['chemotherapy'], 0.0)

    def test_chemo_decay_no_diffusion(self):
        previous_chemotherapy_at_cell = 10.0
        self.automaton.grid[(4, 4)]['chemotherapy'] = 10.0
        # Set from source
        self.automaton.model_parameters['chemotherapy_decay'] = 1.0
        # Turn off diffusion
        self.automaton.grid['chemotherapy_diffusion_rate'] = np.zeros(self.shape,dtype=float)

        self.automaton.diffusion(True)

        expected_chemotherapy_at_cell = (previous_chemotherapy_at_cell) + self.time_params['time_step'] * (
            - self.automaton.model_parameters['chemotherapy_decay'] * previous_chemotherapy_at_cell
            + 0)
        self.assertEqual(self.automaton.work_grid[(4, 4)]['chemotherapy'], expected_chemotherapy_at_cell)
        self.assertEqual(self.automaton.work_grid[(3, 4)]['chemotherapy'], 0.0)
        self.assertEqual(self.automaton.work_grid[(4, 3)]['chemotherapy'], 0.0)
        self.assertEqual(self.automaton.work_grid[(4, 5)]['chemotherapy'], 0.0)
        self.assertEqual(self.automaton.work_grid[(5, 4)]['chemotherapy'], 0.0)

    def test_chemotherapy_diffusion_differing_rates(self):

        self.automaton.grid['chemotherapy_diffusion_rate'][(4, 4)] = 1.0
        self.automaton.grid['chemotherapy_diffusion_rate'][(3, 4)] = 0.7
        self.automaton.grid['chemotherapy_diffusion_rate'][(4, 3)] = 0.6
        self.automaton.grid['chemotherapy_diffusion_rate'][(4, 5)] = 0.5
        self.automaton.grid['chemotherapy_diffusion_rate'][(5, 4)] = 0.4

        previous_chemotherapy_at_source_cell = 10.0
        self.automaton.grid[(4, 4)]['chemotherapy'] = 10.0

        self.automaton.diffusion(True)

        expected_chemotherapy_at_cell = (previous_chemotherapy_at_source_cell) + self.time_params['time_step'] * (
            (((self.automaton.grid['chemotherapy_diffusion_rate'][(4, 4)] + self.automaton.grid['chemotherapy_diffusion_rate'][
                (3, 4)]) / 2) *
             (0 - previous_chemotherapy_at_source_cell) -
             ((self.automaton.grid['chemotherapy_diffusion_rate'][(4, 4)] + self.automaton.grid['chemotherapy_diffusion_rate'][
                 (5, 4)]) / 2) *
             (previous_chemotherapy_at_source_cell - 0)) / self.model_params['spatial_step'] ** 2
            + (((self.automaton.grid['chemotherapy_diffusion_rate'][(4, 4)] + self.automaton.grid['chemotherapy_diffusion_rate'][
                (4, 3)]) / 2) *
               (0 - previous_chemotherapy_at_source_cell) -
               ((self.automaton.grid['chemotherapy_diffusion_rate'][(4, 4)] + self.automaton.grid['chemotherapy_diffusion_rate'][
                   (4, 5)]) / 2) *
               (previous_chemotherapy_at_source_cell - 0)) / self.model_params['spatial_step'] ** 2
            + 0 + 0)
        self.assertEqual(self.automaton.work_grid[(4, 4)]['chemotherapy'], expected_chemotherapy_at_cell)
        # Neighbours - check the neighbours get the chemotherapy
        expected_chemotherapy_at_cell = (0) + self.time_params['time_step'] * (
            (((self.automaton.grid['chemotherapy_diffusion_rate'][(3, 4)] + self.automaton.grid['chemotherapy_diffusion_rate'][
                (4, 4)]) / 2) *
             (previous_chemotherapy_at_source_cell - 0)) / self.model_params['spatial_step'] ** 2
            + 0 + 0)
        self.assertEqual(self.automaton.work_grid[(3, 4)]['chemotherapy'], expected_chemotherapy_at_cell)
        expected_chemotherapy_at_cell = (0) + self.time_params['time_step'] * (
            (((self.automaton.grid['chemotherapy_diffusion_rate'][(5, 4)] + self.automaton.grid['chemotherapy_diffusion_rate'][
                (4, 4)]) / 2) *
             (previous_chemotherapy_at_source_cell - 0)) / self.model_params['spatial_step'] ** 2
            + 0 + 0)
        self.assertEqual(self.automaton.work_grid[(5, 4)]['chemotherapy'], expected_chemotherapy_at_cell)
        expected_chemotherapy_at_cell = (0) + self.time_params['time_step'] * (
            (((self.automaton.grid['chemotherapy_diffusion_rate'][(4, 3)] + self.automaton.grid['chemotherapy_diffusion_rate'][
                (4, 4)]) / 2) *
             (previous_chemotherapy_at_source_cell - 0)) / self.model_params['spatial_step'] ** 2
            + 0 + 0)
        self.assertEqual(self.automaton.work_grid[(4, 3)]['chemotherapy'], expected_chemotherapy_at_cell)
        expected_chemotherapy_at_cell = (0) + self.time_params['time_step'] * (
            (((self.automaton.grid['chemotherapy_diffusion_rate'][(4, 5)] + self.automaton.grid['chemotherapy_diffusion_rate'][
                (4, 4)]) / 2) *
             (previous_chemotherapy_at_source_cell - 0)) / self.model_params['spatial_step'] ** 2
            + 0 + 0)
        self.assertEqual(self.automaton.work_grid[(4, 5)]['chemotherapy'], expected_chemotherapy_at_cell)
        
    def test_chemokine_diffusion(self):

        self.automaton.grid[(4,4)]['chemokine'] = 10.0
        self.automaton.diffusion(True)

        previous_chemo_at_source_cell = 10.0
        expected_chemo_at_cell = (previous_chemo_at_source_cell) + self.time_params['time_step'] * (
            (((self.model_params['chemokine_diffusion'] + self.model_params['chemokine_diffusion']) / 2) *
             (0 - previous_chemo_at_source_cell) -
             ((self.model_params['chemokine_diffusion'] + self.model_params['chemokine_diffusion']) / 2) *
             (previous_chemo_at_source_cell - 0)) / self.model_params['spatial_step'] ** 2
            + (((self.model_params['chemokine_diffusion'] + self.model_params['chemokine_diffusion']) / 2) *
               (0 - previous_chemo_at_source_cell) -
               ((self.model_params['chemokine_diffusion'] + self.model_params['chemokine_diffusion']) / 2) *
               (previous_chemo_at_source_cell - 0)) / self.model_params['spatial_step'] ** 2
            + 0 + 0)
        self.assertEqual(self.automaton.work_grid[(4, 4)]['chemokine'], expected_chemo_at_cell)
        # Neighbours - check the neighbours get the chemokine
        previous_chemokine_at_neighb_cell = 0
        expected_chemo_at_cell = (0) + self.time_params['time_step'] * (
            (((self.model_params['chemokine_diffusion'] + self.model_params['chemokine_diffusion']) / 2) *
             (previous_chemo_at_source_cell - 0) -
             ((self.model_params['chemokine_diffusion'] + self.model_params['chemokine_diffusion']) / 2) *
             (0 - 0)) / self.model_params['spatial_step'] ** 2
            + (((self.model_params['chemokine_diffusion'] + self.model_params['chemokine_diffusion']) / 2) *
               (0 - 0) -
               ((self.model_params['chemokine_diffusion'] + self.model_params['chemokine_diffusion']) / 2) *
               (0 - 0)) / self.model_params['spatial_step'] ** 2
            + 0 + 0)
        self.assertEqual(self.automaton.work_grid[(3, 4)]['chemokine'], expected_chemo_at_cell)
        self.assertEqual(self.automaton.work_grid[(4, 3)]['chemokine'], expected_chemo_at_cell)
        self.assertEqual(self.automaton.work_grid[(4, 5)]['chemokine'], expected_chemo_at_cell)
        self.assertEqual(self.automaton.work_grid[(5, 4)]['chemokine'], expected_chemo_at_cell)

    def test_chemokine_no_diffusion(self):

        self.automaton.model_parameters['chemokine_diffusion'] = 0.0
        self.automaton.grid[(4,4)]['chemokine'] = 10.0
        self.automaton.diffusion(True)

        previous_chemo_at_source_cell = 10.0
        self.assertEqual(self.automaton.work_grid[(4, 4)]['chemokine'], 10.0)
        # Neighbours - check the neighbours get the chemokine
        self.assertEqual(self.automaton.work_grid[(3, 4)]['chemokine'], 0.0)
        self.assertEqual(self.automaton.work_grid[(4, 3)]['chemokine'], 0.0)
        self.assertEqual(self.automaton.work_grid[(4, 5)]['chemokine'], 0.0)
        self.assertEqual(self.automaton.work_grid[(5, 4)]['chemokine'], 0.0)

    def test_chemokine_decay_diffusion(self):

        self.automaton.model_parameters['chemokine_diffusion'] = 0.0
        self.automaton.model_parameters['chemokine_decay'] = 1.0
        self.automaton.grid[(4,4)]['chemokine'] = 10.0
        self.automaton.diffusion(True)

        previous_chemo_at_source_cell = 10.0
        expected_chemo_at_cell = (previous_chemo_at_source_cell) + self.time_params['time_step'] * (
            -previous_chemo_at_source_cell * self.automaton.model_parameters['chemokine_decay']
            + 0 + 0
            )
        self.assertEqual(self.automaton.work_grid[(4, 4)]['chemokine'], expected_chemo_at_cell)
        # Neighbours - check the neighbours get the chemokine
        self.assertEqual(self.automaton.work_grid[(3, 4)]['chemokine'], 0.0)
        self.assertEqual(self.automaton.work_grid[(4, 3)]['chemokine'], 0.0)
        self.assertEqual(self.automaton.work_grid[(4, 5)]['chemokine'], 0.0)
        self.assertEqual(self.automaton.work_grid[(5, 4)]['chemokine'], 0.0)

    def test_chemokine_from_nonresting_mac_diffusion(self):

        self.automaton.model_parameters['chemokine_diffusion'] = 0.0
        self.automaton.model_parameters['chemokine_from_macrophage'] = 1.0
        m = Macrophage((4,4), 'active')
        self.automaton.grid[(4,4)]['contents'] = m
        self.automaton.grid[(4,4)]['chemokine'] = 10.0

        self.automaton.diffusion(True)

        previous_chemo_at_source_cell = 10.0
        expected_chemo_at_cell = (previous_chemo_at_source_cell) + self.time_params['time_step'] * (
            + self.automaton.model_parameters['chemokine_from_macrophage']
            + 0 + 0
            )
        self.assertEqual(self.automaton.work_grid[(4, 4)]['chemokine'], expected_chemo_at_cell)
        # Neighbours - check the neighbours get the chemokine
        self.assertEqual(self.automaton.work_grid[(3, 4)]['chemokine'], 0.0)
        self.assertEqual(self.automaton.work_grid[(4, 3)]['chemokine'], 0.0)
        self.assertEqual(self.automaton.work_grid[(4, 5)]['chemokine'], 0.0)
        self.assertEqual(self.automaton.work_grid[(5, 4)]['chemokine'], 0.0)

    def test_chemokine_resting_mac_diffusion(self):
        # Resting mac so should be no change
        self.automaton.model_parameters['chemokine_diffusion'] = 0.0
        self.automaton.model_parameters['chemokine_from_macrophage'] = 1.0
        m = Macrophage((4,4), 'resting')
        self.automaton.grid[(4, 4)]['contents'] = m
        self.automaton.grid[(4, 4)]['chemokine'] = 10.0

        self.automaton.diffusion(True)

        self.assertEqual(self.automaton.work_grid[(4, 4)]['chemokine'], 10.0)
        # Neighbours - check the neighbours get the chemokine
        self.assertEqual(self.automaton.work_grid[(3, 4)]['chemokine'], 0.0)
        self.assertEqual(self.automaton.work_grid[(4, 3)]['chemokine'], 0.0)
        self.assertEqual(self.automaton.work_grid[(4, 5)]['chemokine'], 0.0)
        self.assertEqual(self.automaton.work_grid[(5, 4)]['chemokine'], 0.0)

    # TODO need to check edges cases as well :(

if __name__ == '__main__':
    unittest.main()
