from TBAutomaton.TBAutomaton import *
import cProfile
import numpy as np

import ConfigParser
import os
import time


def initialise():
    """
    Convert the defined lists/settings into a series of coordinates for initialisation of the automaton
    :param config:
    :param total_shape:
    :return:
    """

    available_addresses = []
    for a in itertools.product(range(total_shape[0]), range(total_shape[1])):
        available_addresses.append(a)

    # BLOOD_VESSELS
    blood_vessels_method = config.get("InitialiseSection", "blood_vessels")
    blood_vessel_addresses = []

    if blood_vessels_method == 'hard_code':
        bv_list = config.get("InitialiseSection", 'blood_vessels_hard_code').split('/')
        for b in bv_list:
            address = tuple(int(c) for c in b.split(","))
            available_addresses.remove(address)
            blood_vessel_addresses.append(address)
    elif blood_vessels_method == 'from_file':
        path = config.get("InitialiseSection", 'blood_vessels_from_file')
        # Add the values to the list
        list_of_vessels = [line.rstrip('\n') for line in open(path)]
        for index in range(len(list_of_vessels)):
            if float(list_of_vessels[index]) > 0.0:
                try:
                    address = np.unravel_index(index, total_shape)
                except ValueError:
                    raise Exception(index, total_shape)
                blood_vessel_addresses.append(address)
                available_addresses.remove(address)
    elif blood_vessels_method == 'random':
        number = config.getint("InitialiseSection", "blood_vessels_random_number")
        assert len(available_addresses) > number
        for i in range(number):
            address = available_addresses.pop(np.random.randint(0, len(available_addresses)))
            blood_vessel_addresses.append(address)

    # FAST BACTERIA
    bacteria_method = config.get("InitialiseSection", "bacteria")
    fast_addresses = []
    slow_addresses = []

    if bacteria_method == 'hard_code':
        fast_list = config.get("InitialiseSection", "bacteria_fast_hard_code").split("/")

        assert len(available_addresses) > len(fast_list)
        for a in fast_list:
            address = tuple(int(c) for c in a.split(","))
            if address in available_addresses:
                fast_addresses.append(address)
                available_addresses.remove(address)
            else:
                # TODO - avoid conflict
                pass

        slow_list = config.get("InitialiseSection", "bacteria_slow_hard_code").split("/")
        assert len(available_addresses) > len(slow_list)
        for a in slow_list:
            address = tuple(int(c) for c in a.split(","))
            if address not in blood_vessel_addresses:
                slow_addresses.append(address)
                available_addresses.remove(address)
            else:
                # TODO - avoid conflict
                pass
    elif bacteria_method == 'random':
        number_fast = config.getint("InitialiseSection", "bacteria_fast_random_number")
        assert len(available_addresses) > number_fast
        for i in range(number_fast):
            address = available_addresses.pop(np.random.randint(0, len(available_addresses)))
            fast_addresses.append(address)

        number_slow = config.getint("InitialiseSection", "bacteria_slow_random_number")
        assert len(available_addresses) > number_slow
        for i in range(number_slow):
            address = available_addresses.pop(np.random.randint(0, len(available_addresses)))
            slow_addresses.append(address)

    # MACROPHAGES
    macrophage_method = config.get("InitialiseSection", "macrophages")
    macrophage_addresses = []
    if macrophage_method == 'random':
        number = config.getint("InitialiseSection", "macrophages_random_number")

        # Make sure there's enough room
        assert len(available_addresses) > number

        for i in range(number):
            address = available_addresses.pop(np.random.randint(0, len(available_addresses)))
            macrophage_addresses.append(address)
    # TODO: hard-code macrophages?

    return blood_vessel_addresses, fast_addresses, slow_addresses, macrophage_addresses

print '------------------------'
print 'TB Simulation Automaton'
print '------------------------'
whole_start_time = time.time()
print "Begin:   {", whole_start_time, "}"

config = ConfigParser.RawConfigParser()
if not config.read('config.properties'):
    raise IOError("Config file (config.properties) not found")

# LOAD PARAMETERS
parameters = {}
# Get all options in parameters section and add to the dictionary
for i in config.options("ModelParametersSection"):
    parameters[i] = config.getfloat("ModelParametersSection", i)

# TIME PARAMETERS
time_parameters = {}
# Get all options in time parameters section
for i in config.options("TimeParametersSection"):
    time_parameters[i] = config.getfloat("TimeParametersSection", i)

# LOAD GRID ATTRIBUTES
total_shape = [int(a) for a in config.get("GridSection", "total_shape").split(",")]

# LOAD RUN PARAMETERS
profile = config.getboolean("RunParametersSection", "profile")
numpy_seed = None
if not config.getboolean("RunParametersSection", "random"):
    numpy_seed = config.getint("RunParametersSection", "non_random_seed")

number_of_runs = config.getint("RunParametersSection", "number_runs")

main_output_location = config.get("RunParametersSection", "output_location")
# Make the main output folder
if not os.path.exists(main_output_location):
    os.makedirs(main_output_location)

random = config.getboolean("RunParametersSection", "random")
debug = config.getboolean("RunParametersSection", "debug")

# LOAD INITIALISATION
blood_vessels, fast_bacteria, slow_bacteria, macrophages = initialise()

for n in range(number_of_runs):
    output_location = main_output_location + "/" + str(n)
    if not os.path.exists(output_location):
        os.makedirs(output_location)

    if random:
        numpy_seed = config.getint("RunParametersSection", "non_random_seed")
        automaton = TBAutomaton(total_shape, time_parameters, parameters, output_location, blood_vessels, macrophages,
                            fast_bacteria, slow_bacteria, numpy_seed=numpy_seed, debug=debug)
    else:
        automaton = TBAutomaton(total_shape, time_parameters, parameters, output_location, blood_vessels, macrophages,
                                fast_bacteria, slow_bacteria, debug=debug)

    if profile:
        pr = cProfile.Profile()
        pr.enable()
        automaton.run()
        pr.disable()
        pr.print_stats(sort='cumtime')
    else:
        automaton.run()

whole_end_time = time.time()
print "End:     {", whole_end_time, "}"
print "Duration: ", whole_end_time - whole_start_time
