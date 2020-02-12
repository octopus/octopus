#!/usr/bin/env python3
import sys
import os
from string import Template
import yaml
import json
import hashlib
import glob


def get_hash(obj):
    return hashlib.sha1(
        json.dumps(obj, sort_keys=True).encode('ascii')).hexdigest()


def mkdir_p(path):
    os.makedirs(path, exist_ok=True)


def outer_product_from_dict(d):
    '''Generate list of dictionaries as outer product from d'''
    result = [[]]
    for key in d:
        result = [x+[y] for x in result for y in d[key]]
    return [{key: value for key, value in zip(d.keys(), r)}
            for r in result]


def keep_combination(combination, test_name):
    """Decides if we want to keep a combination and also modifies it

    The function gets the combination with the relevant input parameters
    and also the name of the test to be able to adapt to it.
    """
    if combination['spin_components'] == 'unpolarized':
        combination['electrons_per_state'] = 2
    if combination['spin_components'] == 'spinors':
        combination['test_type'] = 'complex'
    if combination['periodic_dimensions'] == 0 and \
            combination['number_kpoints'] > 1:
        return False
    if combination['dimensions'] == 2 and \
            combination['periodic_dimensions'] not in [0, 2]:
        return False
    if combination['dimensions'] == 4 and \
            combination['periodic_dimensions'] != 0:
        return False
    return True


def modify_combinations(combinations, test_name):
    elements_to_keep = []
    for index, combination in enumerate(combinations):
        if keep_combination(combination, test_name):
            elements_to_keep.append(index)
    combinations = [c for i, c in enumerate(combinations)
                    if i in elements_to_keep]
    return combinations


def get_combinations(test_name, path='tests'):
    filename = '{}.combinations.yaml'.format(test_name)
    with open(os.path.join(path, filename), 'r') as f_in:
        data = yaml.safe_load(f_in)
    combinations_dict = data['parameters']
    for parameter in combinations_dict:
        # check environment variable to override
        variable = 'OPRT_' + parameter
        if variable in os.environ:
            print("Information: {} is overridden by environment variable"
                  .format(parameter))
            combinations_dict[parameter] = yaml.safe_load(os.environ[variable])
    combinations_list = modify_combinations(
        outer_product_from_dict(combinations_dict), test_name)
    combinations = {}
    for combination in combinations_list:
        combinations[get_hash(combination)] = combination
    return combinations


def get_timings_tags(test_name, path='tests'):
    filename = '{}.combinations.yaml'.format(test_name)
    with open(os.path.join(path, filename), 'r') as f_in:
        data = yaml.safe_load(f_in)
    return data['timings_tags']


def get_template(test_name, path='tests'):
    filename = '{}.inp'.format(test_name)
    with open(os.path.join(path, filename), 'r') as f_in:
        input_template = Template(f_in.read())
    return input_template


def write_input_files(combinations, test_name, template, test_path):
    timing_tags = get_timings_tags(test_name, test_path)
    test_path = os.path.join('runs', test_name)
    for hash, combination in combinations.items():
        path = os.path.join(test_path, hash)
        mkdir_p(path)
        with open(os.path.join(path, 'inp'), 'w') as inputfile:
            inputfile.write(template.substitute(**combination))
        timing_path = os.path.join(path, 'profiling')
        mkdir_p(timing_path)
        with open(os.path.join(timing_path,
                               'process_timings.sh'), 'w') as inputfile:
            inputfile.write("head -n2 time.000000.yaml > time.yaml\n")
            for tag in timing_tags:
                inputfile.write("grep {} time.000000.yaml >> time.yaml\n"
                                .format(tag))
            inputfile.write("echo > /dev/null\n")
    combination_list = []
    for hash, combination in combinations.items():
        combination_list.append({'hash': hash, **combination})
    with open(os.path.join(test_path, 'combinations.yaml'), 'w') as output:
        yaml.safe_dump(combination_list, output)


def write_make_targets(combinations, test_name, filename='targets.inc'):
    test_path = os.path.join('runs', test_name)
    target_paths = []
    for hash, combination in combinations.items():
        path = os.path.join(test_path, hash)
        target_paths.append(os.path.join(path, 'profiling', 'time.yaml'))
    targets_string = ' '.join(target_paths)
    with open(os.path.join(test_path, filename), 'w') as make_targets:
        make_targets.write('targets += ' + targets_string + '\n')


def create_test(test_name, test_path):
    combinations = get_combinations(test_name, test_path)
    template = get_template(test_name, test_path)
    write_input_files(combinations, test_name, template, test_path)
    write_make_targets(combinations, test_name)
    print("Created {} runs for test {}.".format(
        len(combinations), test_name))


def get_test_names(path='tests'):
    test_files = glob.glob(os.path.join(path, '*.inp'))
    test_names = []
    for full_filename in test_files:
        filename = os.path.split(full_filename)[1]
        test_name = os.path.splitext(filename)[0]
        test_names.append(test_name)
    return test_names


def create_central_targets():
    with open('targets.inc', 'w') as make_targets:
        make_targets.write('include runs/*/*inc\n')


if __name__ == '__main__':
    if "testsuite" in os.environ:
        test_path = os.path.join(os.environ["testsuite"], "tests")
    else:
        test_path = "tests"
    if len(sys.argv) == 1:
        test_names = get_test_names(path=test_path)
    elif sys.argv[1] == 'all':
        test_names = get_test_names(path=test_path)
    else:
        test_names = sys.argv[1:]
    for test_name in test_names:
        create_test(test_name, test_path)
    create_central_targets()
