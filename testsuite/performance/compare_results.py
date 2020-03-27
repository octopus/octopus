#!/usr/bin/env python3
import sys
import yaml

if len(sys.argv) != 3:
    print("Error: need two arguments, reference and current yaml files.")
    sys.exit(1)

with open(sys.argv[1]) as reference_file:
    reference_data = list(yaml.safe_load_all(reference_file))

with open(sys.argv[2]) as current_file:
    current_data = list(yaml.safe_load_all(current_file))

time_data = {}
for reference_test, current_test in zip(reference_data, current_data):
    index = reference_test['schema'].index('total_time')
    time_data[reference_test['run']] = {}
    for key in reference_test['data'].keys():
        reference_time = reference_test['data'][key][index]
        current_time = current_test['data'][key][index]
        time_data[reference_test['run']][key] = {
                'reference': reference_time, 'current': current_time}

num_slower_tests = 0
num_tests = 0
for run, times in time_data.items():
    for key, time in times.items():
        num_tests += 1
        if time['current'] < time['reference']:
            print("Run {}, key {} is slower than reference "
                  "(by {:.1%}, {} s vs {} s).".format(
                      run, key, 1-time['current']/time['reference'],
                      time['current'], time['reference']))
            num_slower_tests += 1

if num_slower_tests > 0:
    print("Error: {:d} runs of {:d} are slower than the reference.".format(
        num_slower_tests, num_tests))
    sys.exit(1)
