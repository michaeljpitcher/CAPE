from CAPE.Event import *
from TBAgents import *


class BacteriumReplication(Event):
    def __init__(self, original_bac_address, new_bac_address, new_metabolism):
        Event.__init__(self, [original_bac_address, new_bac_address], [new_bac_address])
        self.original_bac_address = original_bac_address
        self.new_bac_address = new_bac_address
        self.new_metabolism = new_metabolism

    def perform_event(self, tb_automaton):
        new_bacterium = Bacterium(self.new_bac_address, self.new_metabolism)
        tb_automaton.bacteria.append(new_bacterium)
        tb_automaton.work_grid[(self.new_bac_address)]['contents'] = new_bacterium

        original_bacterium = tb_automaton.grid[self.original_bac_address]['contents']
        if original_bacterium.division_neighbourhood == 'mo':
            original_bacterium.division_neighbourhood = 'vn'
        else:
            original_bacterium.division_neighbourhood = 'mo'


class BacteriumStateChange(Event):
    def __init__(self, address, attribute, value):
        # No impacted addresses - changing state doesn't prevent other events (e.g. being ingested)
        Event.__init__(self, [address], [])
        self.bacterium_address = address
        self.attribute = attribute
        self.value = value

    def perform_event(self, automaton):
        bacterium = automaton.grid[self.bacterium_address]['contents']
        if self.attribute == 'metabolism':
            bacterium.metabolism = self.value
        elif self.attribute == 'resting':
            bacterium.resting = self.value


class RecruitTCell(Event):
    def __init__(self, bv_address, new_t_cell_address):
        # Vessel in impacted addresses to ensure only 1 thing recruited at a step
        Event.__init__(self, [bv_address, new_t_cell_address], [bv_address, new_t_cell_address])
        self.blood_vessel_address = bv_address
        self.new_t_cell_address = new_t_cell_address

    def perform_event(self, automaton):
        new_t_cell = TCell(self.new_t_cell_address)
        automaton.t_cells.append(new_t_cell)
        automaton.work_grid[self.new_t_cell_address]['contents'] = new_t_cell


class RecruitMacrophage(Event):
    def __init__(self, bv_address, new_macrophage_address):
        Event.__init__(self, [bv_address, new_macrophage_address], [bv_address, new_macrophage_address])
        self.blood_vessel_address = bv_address
        self.new_macrophage_address = new_macrophage_address

    def perform_event(self, automaton):
        new_macrophage = Macrophage(self.new_macrophage_address, 'resting')
        automaton.macrophages.append(new_macrophage)
        automaton.work_grid[self.new_macrophage_address]['contents'] = new_macrophage


class ChemoKillBacterium(Event):
    def __init__(self, bac_address):
        Event.__init__(self, [bac_address], [bac_address])
        self.bacterium_address = bac_address

    def perform_event(self, automaton):
        bacterium = automaton.grid[self.bacterium_address]['contents']
        automaton.bacteria.remove(bacterium)
        automaton.grid[self.bacterium_address]['contents'] = 0


class ChemoKillMacrophage(Event):
    def __init__(self, mac_address):
        Event.__init__(self, [mac_address], [mac_address])
        self.macrophage_address = mac_address

    def perform_event(self, automaton):
        macrophage = automaton.grid[self.macrophage_address]['contents']
        automaton.macrophages.remove(macrophage)
        caseum = Caseum(self.macrophage_address)
        automaton.caseum_addresses.append(self.macrophage_address)
        automaton.work_grid[self.macrophage_address]['contents'] = caseum


class TCellDeath(Event):
    def __init__(self, t_cell_address):
        Event.__init__(self, [t_cell_address], [t_cell_address])
        self.t_cell_address = t_cell_address

    def perform_event(self, automaton):
        t_cell = automaton.grid[self.t_cell_address]['contents']
        automaton.t_cells.remove(t_cell)
        automaton.work_grid[self.t_cell_address]['contents'] = 0


class TCellMovement(Event):
    def __init__(self, tcell_from_address, tcell_to_address):
        Event.__init__(self, [tcell_from_address, tcell_to_address], [tcell_from_address, tcell_to_address])
        self.tcell_from_address = tcell_from_address
        self.tcell_to_address = tcell_to_address

    def perform_event(self, automaton):
        t_cell = automaton.grid[self.tcell_from_address]['contents']
        t_cell.address = self.tcell_to_address
        automaton.work_grid[self.tcell_from_address]['contents'] = 0
        automaton.work_grid[self.tcell_to_address]['contents'] = t_cell


class TCellKillsMacrophage(Event):
    def __init__(self, tcell_address, macrophage_address):
        Event.__init__(self, [tcell_address, macrophage_address], [tcell_address, macrophage_address])
        self.tcell_address = tcell_address
        self.macrophage_address = macrophage_address

    def perform_event(self, automaton):
        t_cell = automaton.grid[self.tcell_address]['contents']
        macrophage = automaton.grid[self.macrophage_address]['contents']

        automaton.t_cells.remove(t_cell)
        automaton.work_grid[self.tcell_address]['contents'] = 0
        automaton.macrophages.remove(macrophage)
        caseum = Caseum(self.macrophage_address)
        automaton.caseum_addresses.append(self.macrophage_address)
        automaton.work_grid[self.macrophage_address]['contents'] = caseum


class MacrophageDeath(Event):
    def __init__(self, macrophage_address):
        Event.__init__(self, [macrophage_address], [macrophage_address])
        self.macrophage_address = macrophage_address

    def perform_event(self, automaton):
        macrophage = automaton.grid[self.macrophage_address]['contents']
        automaton.macrophages.remove(macrophage)
        if macrophage.state == 'infected' or macrophage.state == 'chronically_infected':
            caseum = Caseum(self.macrophage_address)
            automaton.caseum_addresses.append(self.macrophage_address)
            automaton.work_grid[self.macrophage_address]['contents'] = caseum
        else:
            automaton.work_grid[self.macrophage_address]['contents'] = 0


class MacrophageMovement(Event):
    def __init__(self, macrophage_from_address, macrophage_to_address):
        Event.__init__(self, [macrophage_from_address, macrophage_to_address],
                       [macrophage_from_address, macrophage_to_address])
        self.macrophage_from_address = macrophage_from_address
        self.macrophage_to_address = macrophage_to_address

    def perform_event(self, automaton):
        macrophage = automaton.grid[self.macrophage_from_address]['contents']
        macrophage.address = self.macrophage_to_address
        automaton.work_grid[self.macrophage_from_address]['contents'] = 0
        automaton.work_grid[self.macrophage_to_address]['contents'] = macrophage


class MacrophageIngestsBacterium(Event):
    def __init__(self, macrophage_address, bacterium_address):
        Event.__init__(self, [macrophage_address, bacterium_address], [macrophage_address, bacterium_address])
        self.macrophage_address = macrophage_address
        self.bacterium_address = bacterium_address

    def perform_event(self, automaton):
        macrophage = automaton.grid[self.macrophage_address]['contents']
        bacterium = automaton.grid[self.bacterium_address]['contents']

        macrophage.address = self.bacterium_address
        automaton.bacteria.remove(bacterium)
        automaton.work_grid[self.macrophage_address]['contents'] = 0
        automaton.work_grid[self.bacterium_address]['contents'] = macrophage

        # If not active, intracellular bacteria count increases by 1
        if macrophage.state != 'active':
            macrophage.intracellular_bacteria += 1
            # Resting macrophages become infected
            if macrophage.state == 'resting':
                macrophage.state = 'infected'
            # Infected macrophages become chronically infected if they breach threshold
            elif macrophage.state == 'infected' and macrophage.intracellular_bacteria == \
                automaton.model_parameters['bacteria_to_turn_chronically_infected']:
                macrophage.state = 'chronically_infected'


class MacrophageActivation(Event):
    def __init__(self, mac_address, state):
        Event.__init__(self, [mac_address], [mac_address])
        self.macrophage_address = mac_address
        self.new_state = state

    def perform_event(self, automaton):
        automaton.grid[self.macrophage_address]['contents'].state = self.new_state


class MacrophageBursts(Event):
    def __init__(self, mac_address, new_bacteria_addresses):
        # Bacteria addresses are impacted, but they're not dependent (if something else moves into a cell where a
        # bacterium would be deposited, this doesn't stop the macrophage bursting)
        Event.__init__(self, [mac_address], [mac_address] + new_bacteria_addresses)
        self.macrophage_address = mac_address
        self.new_bacteria_addresses = new_bacteria_addresses

    def perform_event(self, automaton):
        macrophage = automaton.grid[self.macrophage_address]['contents']
        automaton.macrophages.remove(macrophage)
        caseum = Caseum(self.macrophage_address)
        automaton.caseum_addresses.append(self.macrophage_address)
        automaton.work_grid[self.macrophage_address]['contents'] = caseum

        for address in self.new_bacteria_addresses:
            # Check if the event is still in the impacted addresses (will have been removed if something else has
            # affected it)
            if address in self.impacted_addresses:
                bac = Bacterium(address, 'slow')
                automaton.bacteria.append(bac)
                automaton.work_grid[address]['contents'] = bac

