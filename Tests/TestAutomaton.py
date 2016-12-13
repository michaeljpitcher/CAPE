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

    def test_neighbour_addresses(self):
        pass

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

    def test_neighbour_relatives(self):
        self.assertItemsEqual(self.automaton.moore_relative.keys(), [1, 2, 3])
        self.assertItemsEqual(self.automaton.moore_relative[1],
                              [(-1, -1), (-1, 0), (-1, 1), (0, -1), (0, 1), (1, -1),
                               (1, 0), (1, 1)])
        self.assertItemsEqual(self.automaton.moore_relative[2],
                              [(-2, -2), (-2, -1), (-2, 0), (-2, 1), (-2, 2), (-1, -2), (-1, 2), (0, -2), (0, 2),
                               (1, -2), (1, 2), (2, -2), (2, -1), (2, 0), (2, 1), (2, 2)])
        self.assertItemsEqual(self.automaton.moore_relative[3],
                              [(-3, -3), (-3, -2), (-3, -1), (-3, 0), (-3, 1), (-3, 2), (-3, 3), (-2, -3), (-2, 3),
                               (-1, -3), (-1, 3), (0, -3), (0, 3), (1, -3), (1, 3), (2, -3), (2, 3), (3, -3), (3, -2),
                               (3, -1), (3, 0), (3, 1), (3, 2), (3, 3)])

        self.assertItemsEqual(self.automaton.von_neumann_relative.keys(), [1, 2, 3])
        self.assertItemsEqual(self.automaton.von_neumann_relative[1],
                              [(-1, 0), (0, -1), (0, 1), (1, 0)])
        self.assertItemsEqual(self.automaton.von_neumann_relative[2],
                              [(-2, 0), (-1, -1), (-1, 1), (0, -2), (0, 2), (1, -1), (1, 1), (2, 0)])
        self.assertItemsEqual(self.automaton.von_neumann_relative[3],
                              [(-3, 0), (-2, -1), (-2, 1), (-1, -2), (-1, 2), (0, -3), (0, 3), (1, -2), (1, 2), (2, -1),
                               (2, 1), (3, 0)])

    def test_moore_neighbours(self):
        neighbours_0_0_1 = self.automaton.moore_neighbours((0,0), 1)
        self.assertItemsEqual(neighbours_0_0_1.keys(), [(0,1),(1,0),(1,1)])
        neighbours_9_9_1 = self.automaton.moore_neighbours((9, 9), 1)
        self.assertItemsEqual(neighbours_9_9_1.keys(), [(8, 9), (9, 8), (8, 8)])
        neighbours_0_9_1 = self.automaton.moore_neighbours((0, 9), 1)
        self.assertItemsEqual(neighbours_0_9_1.keys(), [(1, 9), (0, 8), (1, 8)])
        neighbours_9_0_1 = self.automaton.moore_neighbours((9, 0), 1)
        self.assertItemsEqual(neighbours_9_0_1.keys(), [(9, 1), (8, 0), (8, 1)])
        neighbours_0_5_1 = self.automaton.moore_neighbours((0, 5), 1)
        self.assertItemsEqual(neighbours_0_5_1.keys(), [(0,4), (0,6), (1,4), (1,5), (1,6)])
        neighbours_5_0_1 = self.automaton.moore_neighbours((5, 0), 1)
        self.assertItemsEqual(neighbours_5_0_1.keys(), [(4,0), (6,0), (4,1), (5,1), (6,1)])
        neighbours_5_9_1 = self.automaton.moore_neighbours((5, 9), 1)
        self.assertItemsEqual(neighbours_5_9_1.keys(), [(4, 9), (6, 9), (4, 8), (5, 8), (6, 8)])
        neighbours_9_5_1 = self.automaton.moore_neighbours((9, 5), 1)
        self.assertItemsEqual(neighbours_9_5_1.keys(), [(9, 4), (9, 6), (8, 4), (8, 5), (8, 6)])
        neighbours_5_5_1 = self.automaton.moore_neighbours((5, 5), 1)
        self.assertItemsEqual(neighbours_5_5_1.keys(), [(4,4),(4,5),(4,6),(5,4),(5,6),(6,4),(6,5),(6,6)])
        
    def test_von_neumann_neighbours(self):
        neighbours_0_0_1 = self.automaton.von_neumann_neighbours((0,0), 1)
        self.assertItemsEqual(neighbours_0_0_1.keys(), [(0,1),(1,0)])
        neighbours_9_9_1 = self.automaton.von_neumann_neighbours((9, 9), 1)
        self.assertItemsEqual(neighbours_9_9_1.keys(), [(8, 9), (9, 8)])
        neighbours_0_9_1 = self.automaton.von_neumann_neighbours((0, 9), 1)
        self.assertItemsEqual(neighbours_0_9_1.keys(), [(1, 9), (0, 8)])
        neighbours_9_0_1 = self.automaton.von_neumann_neighbours((9, 0), 1)
        self.assertItemsEqual(neighbours_9_0_1.keys(), [(9, 1), (8, 0)])
        neighbours_0_5_1 = self.automaton.von_neumann_neighbours((0, 5), 1)
        self.assertItemsEqual(neighbours_0_5_1.keys(), [(0,4), (0,6), (1,5)])
        neighbours_5_0_1 = self.automaton.von_neumann_neighbours((5, 0), 1)
        self.assertItemsEqual(neighbours_5_0_1.keys(), [(4,0), (6,0), (5,1)])
        neighbours_5_9_1 = self.automaton.von_neumann_neighbours((5, 9), 1)
        self.assertItemsEqual(neighbours_5_9_1.keys(), [(4, 9), (6, 9), (5, 8)])
        neighbours_9_5_1 = self.automaton.von_neumann_neighbours((9, 5), 1)
        self.assertItemsEqual(neighbours_9_5_1.keys(), [(9, 4), (9, 6), (8, 5)])
        neighbours_5_5_1 = self.automaton.von_neumann_neighbours((5, 5), 1)
        self.assertItemsEqual(neighbours_5_5_1.keys(), [(4,5),(5,4),(5,6),(6,5)])

        print self.automaton.moore_neighbours((0,0), 1).keys()
        print self.automaton.moore_neighbours((0, 0), 2).keys()


if __name__ == '__main__':
    unittest.main()
