from CAPE.Agent import *


class Bacterium(Agent):
    def __init__(self, address, metabolism):
        Agent.__init__(self, address)
        self.metabolism = metabolism
        self.resting = False
        self.division_neighbourhood = 'mo'

    def output_code(self):
        # 1.0 == Fast, 2.0 == slow
        code = 1.0 + (self.metabolism == 'slow')
        # Add .5 if resting
        if self.resting:
            code += 0.5
        return code


class Caseum(Agent):
    def __init__(self, address):
        Agent.__init__(self, address)

    def output_code(self):
        return 100.0


class TCell(Agent):
    def __init__(self, address):
        Agent.__init__(self, address)

    def output_code(self):
        return 3.0


class Macrophage(Agent):
    def __init__(self, address, state):
        Agent.__init__(self, address)
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
