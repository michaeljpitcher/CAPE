import unittest
from CAPE.CAPEAutomaton import *


class AutomatonTestCase(unittest.TestCase):

    def setUp(self):
        self.shape = (10,10)
        self.attributes=['a','b','c']
        self.formats = ['float','int','object']
        self.time_params={}
        self.time_params['initial_time'] = 0.0
        self.time_params['time_step'] = 0.1
        self.time_params['time_limit'] = 1.0
        self.model_params = {}
        self.model_params['param1'] = 1
        self.model_params['param2'] = 10
        self.model_params['param3'] = 100
        self.initialise = {}
        self.initialise['a'] = {}
        self.initialise['a'][(1,1)] = 8.8
        self.initialise['b'] = {}
        self.initialise['b'][(2, 2)] = 13
        self.initialise['c'] = {}
        self.initialise['c'][(3, 3)] = dict()
        self.initialise['c'][(3, 3)]['test'] = 99
        self.initialise['c'][(4, 4)] = Agent()
        self.automaton = Automaton(self.shape, self.attributes, self.formats, self.time_params, self.model_params,
                                   self.initialise)

    def test_initialise(self):
        self.assertItemsEqual(self.attributes, self.automaton.attributes)
        self.assertItemsEqual(self.model_params, self.automaton.model_parameters)
        self.assertEqual(self.automaton.time, self.time_params['initial_time'])
        self.assertEqual(self.automaton.time_step, self.time_params['time_step'])
        self.assertEqual(self.automaton.time_limit, self.time_params['time_limit'] / self.time_params['time_step'])

        self.assertSequenceEqual(self.automaton.grid.shape, self.shape)
        self.assertSequenceEqual(self.automaton.work_grid.shape, self.shape)

        for x in range(10):
            for y in range(10):
                # A
                if (x == 1 and y == 1):
                    self.assertEqual(self.automaton.grid[(x,y)]['a'], 8.8)
                else:
                    self.assertEqual(self.automaton.grid[(x,y)]['a'], 0.0)
                # B
                if (x == 2 and y == 2):
                    self.assertEqual(self.automaton.grid[(x, y)]['b'], 13)
                else:
                    self.assertEqual(self.automaton.grid[(x, y)]['b'], 0)
                # C
                if (x == 3 and y == 3):
                    self.assertItemsEqual(self.automaton.grid[(x,y)]['c'].keys(), ['test'])
                    self.assertEqual(self.automaton.grid[(x, y)]['c']['test'], 99)
                elif (x == 4 and y == 4):
                    self.assertTrue(isinstance(self.automaton.grid[(x, y)]['c'], Agent))

    def test_run_not_overriden(self):
        with self.assertRaises(NotImplementedError) as context:
            self.automaton.run()

    def test_run_override(self):
        class TestAutomaton(Automaton):
            def __init__(self, shape, attributes, formats, time_parameters, model_parameters, initialisation):
                Automaton.__init__(self, shape, attributes, formats, time_parameters, model_parameters, initialisation)

            def update_agents(self):
                pass

            def update_cells(self):
                pass

        test_automaton = TestAutomaton(self.shape, self.attributes, self.formats, self.time_params, self.model_params,
                                       self.initialise)

        test_automaton.run()
        self.assertEqual(test_automaton.time, self.time_params['initial_time'] +
                         self.time_params['time_limit'] / self.time_params['time_step'])



if __name__ == '__main__':
    unittest.main()
