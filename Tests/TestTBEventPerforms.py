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
        self.model_params['bacteria_to_turn_chronically_infected'] = 999

        self.bv = [(1,1)]
        self.macs = [(1,8), (2,8), (3,8), (4,8)]
        self.fb = [(8, 1)]
        self.sb = [(8, 2)]

        self.output_loc = 'test_output'
        if not os.path.exists(self.output_loc):
            os.makedirs(self.output_loc)
        self.automaton = TBAutomaton(self.shape, self.time_params, self.model_params, self.output_loc,
                                     self.bv, self.macs, self.fb, self.sb)

        # Add a t-cell
        tcell = TCell((5,5))
        self.automaton.grid[(5,5)]['contents'] = tcell
        self.automaton.t_cells.append(tcell)

        # Set macrophage states
        self.automaton.grid[(2, 8)]['contents'].state = 'active'
        self.automaton.grid[(3, 8)]['contents'].state = 'infected'
        self.automaton.grid[(3, 8)]['contents'].intracellular_bacteria = 7
        self.automaton.grid[(4, 8)]['contents'].state = 'chronically_infected'
        self.automaton.grid[(4, 8)]['contents'].intracellular_bacteria = 19

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

    def test_macrophage_recruitment_perform(self):
        mac_rec_event = RecruitMacrophage((1, 1), (1, 2))
        mac_rec_event.perform_event(self.automaton)
        self.assertTrue(isinstance(self.automaton.work_grid[(1, 2)]['contents'], Macrophage))
        macrophage = self.automaton.work_grid[(1, 2)]['contents']
        self.assertTrue(macrophage in self.automaton.t_cells)
        self.assertEqual(macrophage.state, 'resting')

    def test_chemo_kill_bacterium(self):
        chem_kill_bac_event = ChemoKillBacterium((8, 1))
        bac = self.automaton.grid[(8,1)]['contents']
        chem_kill_bac_event.perform_event(self.automaton)
        self.assertEqual(self.automaton.work_grid[(8,1)]['contents'], 0.0)
        self.assertTrue(bac not in self.automaton.bacteria)

    def test_chemo_kill_mac(self):
        chem_kill_mac_event = ChemoKillMacrophage((1, 8))
        mac = self.automaton.grid[(1, 8)]['contents']
        chem_kill_mac_event.perform_event(self.automaton)
        self.assertTrue(isinstance(self.automaton.work_grid[(1, 8)]['contents'], Caseum))
        self.assertTrue((1, 8) in self.automaton.caseum_addresses)
        self.assertTrue(mac not in self.automaton.bacteria)

    def test_t_cell_death_perform(self):
        t_cell_dea_event = TCellDeath((5,5))
        t_cell = self.automaton.grid[(5,5)]['contents']
        t_cell_dea_event.perform_event(self.automaton)
        self.assertEqual(self.automaton.work_grid[(5,5)]['contents'], 0.0)
        self.assertTrue(t_cell not in self.automaton.t_cells)

    def test_t_cell_move_perform(self):
        t_cell_move_event = TCellMovement((5,5),(5,4))
        t_cell = self.automaton.grid[(5,5)]['contents']
        t_cell_move_event.perform_event(self.automaton)
        self.assertEqual(self.automaton.work_grid[(5, 5)]['contents'], 0.0)
        self.assertEqual(self.automaton.work_grid[(5, 4)]['contents'], t_cell)

    def test_t_cell_kill_macrophage(self):
        mac = Macrophage((5, 4), 'infected')
        self.automaton.macrophages.append(mac)
        self.automaton.grid[(5, 4)]['contents'] = mac
        t_cell = self.automaton.grid[(5, 5)]['contents']
        tcell_kill_mac_event = TCellKillsMacrophage((5, 5), (5, 4))
        tcell_kill_mac_event.perform_event(self.automaton)
        self.assertTrue(isinstance(self.automaton.work_grid[(5,4)]['contents'], Caseum))
        self.assertTrue((5,4) in self.automaton.caseum_addresses)
        self.assertEqual(self.automaton.work_grid[(5,5)]['contents'], 0.0)
        self.assertTrue(mac not in self.automaton.macrophages)
        self.assertTrue(t_cell not in self.automaton.t_cells)

    def test_macrophage_death_perform(self):
        # resting
        mac = self.automaton.grid[(1, 8)]['contents']
        mac_death_event = MacrophageDeath((1,8))
        mac_death_event.perform_event(self.automaton)
        self.assertEqual(self.automaton.work_grid[(1, 8)]['contents'], 0.0)
        self.assertTrue(mac not in self.automaton.macrophages)
        # active
        mac_death_event = MacrophageDeath((2, 8))
        mac_death_event.perform_event(self.automaton)
        self.assertEqual(self.automaton.work_grid[(2, 8)]['contents'], 0.0)
        # Infected
        mac_death_event = MacrophageDeath((3, 8))
        mac_death_event.perform_event(self.automaton)
        self.assertTrue(isinstance(self.automaton.work_grid[(3, 8)]['contents'], Caseum))
        self.assertTrue((3,8) in self.automaton.caseum_addresses)
        # Chr Infected
        mac_death_event = MacrophageDeath((4, 8))
        mac_death_event.perform_event(self.automaton)
        self.assertTrue(isinstance(self.automaton.work_grid[(4, 8)]['contents'], Caseum))
        self.assertTrue((4, 8) in self.automaton.caseum_addresses)

    def test_macrophage_movement(self):
        mac = self.automaton.grid[(1, 8)]['contents']
        mac_move_event = MacrophageMovement((1,8), (0,8))
        mac_move_event.perform_event(self.automaton)
        self.assertEqual(self.automaton.work_grid[(1,8)]['contents'], 0.0)
        self.assertEqual(self.automaton.work_grid[(0,8)]['contents'], mac)

    def test_resting_macrophage_ingests_bacterium_perform(self):
        bac = Bacterium((1,7), 'fast')
        self.automaton.bacteria.append(bac)
        self.automaton.grid[(1,7)]['contents'] = bac
        mac = self.automaton.grid[(1,8)]['contents']
        mac_ing_bac_event = MacrophageIngestsBacterium((1,8),(1,7))
        mac_ing_bac_event.perform_event(self.automaton)

        self.assertTrue(bac not in self.automaton.bacteria)
        self.assertEqual(self.automaton.work_grid[(1, 8)]['contents'], 0.0)
        self.assertEqual(self.automaton.work_grid[(1, 7)]['contents'], mac)
        self.assertEqual(mac.state, 'infected')
        self.assertEqual(mac.intracellular_bacteria, 1)

    def test_active_macrophage_ingests_bacterium_perform(self):
        bac = Bacterium((2,7), 'fast')
        self.automaton.bacteria.append(bac)
        self.automaton.grid[(2,7)]['contents'] = bac
        mac = self.automaton.grid[(2,8)]['contents']
        mac_ing_bac_event = MacrophageIngestsBacterium((2,8),(2,7))
        mac_ing_bac_event.perform_event(self.automaton)

        self.assertTrue(bac not in self.automaton.bacteria)
        self.assertEqual(self.automaton.work_grid[(2, 8)]['contents'], 0.0)
        self.assertEqual(self.automaton.work_grid[(2, 7)]['contents'], mac)
        self.assertEqual(mac.state, 'active')
        self.assertEqual(mac.intracellular_bacteria, 0)

    def test_infected_macrophage_ingests_bacterium_perform_no_state_change(self):
        bac = Bacterium((3,7), 'fast')
        self.automaton.bacteria.append(bac)
        self.automaton.grid[(3,7)]['contents'] = bac
        mac = self.automaton.grid[(3,8)]['contents']
        orig_int_bac = mac.intracellular_bacteria
        mac_ing_bac_event = MacrophageIngestsBacterium((3,8),(3,7))
        mac_ing_bac_event.perform_event(self.automaton)

        self.assertTrue(bac not in self.automaton.bacteria)
        self.assertEqual(self.automaton.work_grid[(3, 8)]['contents'], 0.0)
        self.assertEqual(self.automaton.work_grid[(3, 7)]['contents'], mac)
        self.assertEqual(mac.state, 'infected')
        self.assertEqual(mac.intracellular_bacteria, orig_int_bac + 1)

    def test_infected_macrophage_ingests_bacterium_perform_change_to_chr_inf(self):

        bac = Bacterium((3,7), 'fast')
        self.automaton.bacteria.append(bac)
        self.automaton.grid[(3,7)]['contents'] = bac
        mac = self.automaton.grid[(3,8)]['contents']
        orig_int_bac = mac.intracellular_bacteria
        self.automaton.model_parameters['bacteria_to_turn_chronically_infected'] = orig_int_bac + 1
        mac_ing_bac_event = MacrophageIngestsBacterium((3,8),(3,7))
        mac_ing_bac_event.perform_event(self.automaton)

        self.assertTrue(bac not in self.automaton.bacteria)
        self.assertEqual(self.automaton.work_grid[(3, 8)]['contents'], 0.0)
        self.assertEqual(self.automaton.work_grid[(3, 7)]['contents'], mac)
        self.assertEqual(mac.state, 'chronically_infected')
        self.assertEqual(mac.intracellular_bacteria, orig_int_bac + 1)

    def test_chr_infected_macrophage_ingests_bacterium_perform(self):
        bac = Bacterium((4,7), 'fast')
        self.automaton.bacteria.append(bac)
        self.automaton.grid[(4,7)]['contents'] = bac
        mac = self.automaton.grid[(4,8)]['contents']
        orig_int_bac = mac.intracellular_bacteria
        mac_ing_bac_event = MacrophageIngestsBacterium((4,8),(4,7))
        mac_ing_bac_event.perform_event(self.automaton)

        self.assertTrue(bac not in self.automaton.bacteria)
        self.assertEqual(self.automaton.work_grid[(4, 8)]['contents'], 0.0)
        self.assertEqual(self.automaton.work_grid[(4, 7)]['contents'], mac)
        self.assertEqual(mac.state, 'chronically_infected')
        self.assertEqual(mac.intracellular_bacteria, orig_int_bac + 1)

    def test_macrophage_activation(self):
        # Activate
        mac = self.automaton.grid[(1,8)]['contents']
        mac_act_event = MacrophageActivation((1,8), 'active')
        mac_act_event.perform_event(self.automaton)
        self.assertEqual(mac.state, 'active')

        # Deactivate
        mac = self.automaton.grid[(2, 8)]['contents']
        mac_act_event = MacrophageActivation((2, 8), 'resting')
        mac_act_event.perform_event(self.automaton)
        self.assertEqual(mac.state, 'resting')

    def test_macrophage_bursting(self):
        mac = self.automaton.grid[(4,8)]['contents']
        mac_burst_event = MacrophageBursts((4,8), [(3,7), (3,9), (4,7), (4,9), (5,7), (5,8), (5,9)])

        # Remove an address from impacted - something else has happened here
        mac_burst_event.impacted_addresses.remove((4,9))

        mac_burst_event.perform_event(self.automaton)

        self.assertTrue(mac not in self.automaton.macrophages)
        self.assertTrue(isinstance(self.automaton.work_grid[(4,8)]['contents'], Caseum))
        self.assertTrue((4,8) in self.automaton.caseum_addresses)

        for address in [(3,7), (3,9), (4,7), (5,7), (5,8), (5,9)]:
            self.assertTrue(isinstance(self.automaton.work_grid[address]['contents'], Bacterium))
            bac = self.automaton.work_grid[address]['contents']
            self.assertTrue(bac in self.automaton.bacteria)
            self.assertEqual(bac.metabolism, 'slow')

        self.assertEqual(self.automaton.work_grid[(4,9)]['contents'], 0.0)


if __name__ == '__main__':
    unittest.main()
