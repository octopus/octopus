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

times = {}
for reference_test, current_test in zip(reference_data, current_data):
    index = reference_test['schema'].index('total_time')
    reference_time = reference_test['data']['COMPLETE_RUN'][index]
    current_time = current_test['data']['COMPLETE_RUN'][index]
    times[reference_test['run']] = {'reference': reference_time,
                                    'current': current_time}

any_test_slower = False
for run, time in times.items():
    if time['current'] < time['reference']:
        print("Run {} is slower than reference ({} vs {}).".format(
            run, time['current'], time['reference']))
        any_test_slower = True

if any_test_slower:
    print("Error: some runs are slower than the reference.")
    sys.exit(1)
