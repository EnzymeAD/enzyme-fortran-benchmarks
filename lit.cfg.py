# -*- Python -*-
# Configuration file for the 'lit' test runner.
#
# Portions adapted from enzyme/test/lit.cfg.py and
# enzyme-benchmarks/lit.cfg.py in the Enzyme project
# (https://github.com/EnzymeAD/Enzyme), which is licensed under the Apache
# License v2.0 with LLVM Exceptions.

import os
import platform

import lit
import lit.formats

# name: The name of this test suite.
config.name = 'BenchEnzymeFortran'

# testFormat: The test format to use to interpret tests.
#
# For now we require '&&' between commands, until they get globally killed and
# the test runner updated.
LLVM_VERSION_MAJOR = lit.__versioninfo__[0]
execute_external = platform.system() != 'Windows'
if LLVM_VERSION_MAJOR >= 23:
    config.test_format = lit.formats.ShTest(execute_external, force_execute_external=True)
else:
    config.test_format = lit.formats.ShTest(execute_external)

# suffixes: A list of file extensions to treat as test files.
config.suffixes = ['.f90']

# test_source_root: The root path where tests are located.
config.test_source_root = os.path.dirname(__file__)

# test_exec_root: The root path where tests should be run.
config.test_exec_root = os.path.join(config.bench_obj_root, 'benchmarks')

# Tweak the PATH to include the tools dir and the scripts dir.
base_paths = [config.llvm_tools_dir, config.environment['PATH']]
config.environment['PATH'] = os.pathsep.join(base_paths)

config.environment['LD_LIBRARY_PATH'] = os.pathsep.join(
    (config.llvm_libs_dir, config.environment.get('LD_LIBRARY_PATH', '')))
