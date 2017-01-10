import unittest
from TBAutomaton.TBAutomaton import *
from TBAutomaton.TBAgents import *
from TBAutomaton.TBEvents import *
import os
import shutil


class EventPerformsTestCase(unittest.TestCase):

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

        self.bv = [(1,1)]
        self.macs = [(1,8), (2,8), (3,8), (4,8)]
        self.fb = [(8, 1)]
        self.sb = [(8, 2)]

        self.output_loc = 'test_output'
        if not os.path.exists(self.output_loc):
            os.makedirs(self.output_loc)
        self.automaton = TBAutomaton(self.shape, self.time_params, self.model_params, self.output_loc,
                                     self.bv, self.macs, self.fb, self.sb)

    def tearDown(self):
        # Close output files and delete
        self.automaton.close_files()
        shutil.rmtree(self.output_loc)

    def test_bacteria_replication_perform(self):

        original_bac_number = len(self.automaton.bacteria)
        old_bac = self.automaton.grid[(8, 1)]['contents']

        bac_rep_event = BacteriumReplication((8, 1), (7, 1), 'fast')
        bac_rep_event.perform_event(self.automaton)

        self.assertTrue(len(self.automaton.bacteria), original_bac_number+1)
        self.assertTrue(isinstance(self.automaton.work_grid[(7, 1)]['contents'], Bacterium))
        new_bac = self.automaton.work_grid[(7, 1)]['contents']
        self.assertEqual(new_bac.metabolism, 'fast')
        self.assertEqual(old_bac.division_neighbourhood, 'vn')
        self.assertEqual(new_bac.division_neighbourhood, 'mo')

    def test_bacterium_state_change_perform(self):

        bac = self.automaton.grid[(8, 1)]['contents']
        self.assertEqual(bac.metabolism, 'fast')
        bac_sta_cha_event = BacteriumStateChange((8, 1), 'metabolism', 'slow')
        bac_sta_cha_event.perform_event(self.automaton)
        self.assertEqual(bac.metabolism, 'slow')

        self.assertEqual(bac.resting, False)
        bac_sta_cha_event = BacteriumStateChange((8, 1), 'resting', True)
        bac_sta_cha_event.perform_event(self.automaton)
        self.assertEqual(bac.resting, True)

    def test_t_cell_recruitment_perform(self):

        t_cell_rec_event = RecruitTCell((1,1), (1,2))
        t_cell_rec_event.perform_event(self.automaton)

        self.assertTrue(isinstance(self.automaton.work_grid[(1,2)]['contents'], TCell))
        tcell = self.automaton.work_grid[(1,2)]['contents']
        self.assertTrue(tcell in self.automaton.t_cells)


if __name__ == '__main__':
    unittest.main()
