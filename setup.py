""" instaler """
# -*- coding: utf-8 -*-
# :Project:   elongation_simulators -- Packaging
# :Author:    Fabio Hedayioglu <fheday@gmail.com>
# :License:   MIT License
# :Copyright: © 2020 Fabio Hedayioglu
#

import sys
import sysconfig
from setuptools import setup
from setuptools.command.build_ext import build_ext
from pybind11.setup_helpers import Pybind11Extension,\
                            ParallelCompile  # noqa:E402


WIN = sys.platform.startswith("win32") and "mingw" not in sysconfig.get_platform()

# Avoid a gcc warning below:
# cc1plus: warning: command line option ‘-Wstrict-prototypes’ is valid
# for C/ObjC but not for C++


class BuildExt(build_ext):
    """ Class to remove warning on gcc."""
    def build_extensions(self):
        if '-Wstrict-prototypes' in self.compiler.compiler_so:
            self.compiler.compiler_so.remove('-Wstrict-prototypes')
        super().build_extensions()


EXTRA_COMPILE_ARGS = None
if WIN:
    EXTRA_COMPILE_ARGS = ["/std:c++17", "/O2", "/Ot", "/GL", "/DCOMIPLE_PYTHON_MODULE", "/I./src/eigen-3.3.7/eigen3/"]
else:
    EXTRA_COMPILE_ARGS = ["-O3", "-ffast-math", "-ftree-vectorize", "-Wall",
                          "-g2", "-flto=auto", "-DCOMIPLE_PYTHON_MODULE"]

ext_modules = [
    Pybind11Extension(
        "sequence_simulator",
        ["src/concentrationsreader.cpp", "src/mrna_reader.cpp", "src/elongation_codon.cpp",
         "src/initiationterminationcodon.cpp", "src/mrnaelement.cpp", "src/sequence_simulator.cpp",
         "src/codon_simulator.cpp", "src/elongation_simulation_manager.cpp",
         "src/elongation_simulation_processor.cpp", "./src/jsoncpp/jsoncpp.cpp"],
        include_dirs=["./src/jsoncpp/", "./src/eigen-3.3.7/"],
        extra_compile_args=EXTRA_COMPILE_ARGS
    ),
    Pybind11Extension(
        "codon_simulator",
        ["src/concentrationsreader.cpp", "src/mrna_reader.cpp", "src/codon_simulator.cpp",
         "./src/jsoncpp/jsoncpp.cpp"],
        include_dirs=["./src/jsoncpp/", "./src/eigen-3.3.7/"],
        extra_compile_args=EXTRA_COMPILE_ARGS
    )
]


if sys.version_info < (3,):
    raise NotImplementedError("Only Python 3+ is supported.")


with open('version.txt', encoding='utf-8') as f:
    VERSION = f.read()

with open('README.md', encoding='utf-8') as f:
    LONG_DESCRIPTION = f.read()

with open('CHANGES.rst', encoding='utf-8') as f:
    CHANGES = f.read()

# Optional multithreaded build
ParallelCompile("NPY_NUM_BUILD_JOBS").install()

CMDCLASS = {}
if not WIN:
    CMDCLASS = {"build_ext": BuildExt}


setup(
    name='elongation_simulator',
    version=VERSION,
    description='High-performance Codon simulator"\
     " and elongation simulator for eukaryotic organism',
    long_description=LONG_DESCRIPTION + '\n\n' + CHANGES,
    long_description_content_type='text/markdown',
    license='MIT License',
    keywords='elongation translation',
    author='Fabio Hedayioglu',
    author_email='fheday@gmail.com',
    maintainer='Fabio Hedayioglu',
    maintainer_email='fheday@gmail.com',
    url='https://github.com/fheday/elongation_simulator/',
    classifiers=[
        'Development Status :: 5 - Production/Stable',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: C++',
        'Programming Language :: Python :: 3 :: Only',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.9',
        'Programming Language :: Python :: 3.10',
        'Programming Language :: Python :: 3.11',
        'Programming Language :: Python :: 3.12',
        'Programming Language :: Python :: 3.13',
    ],
    install_requires=['setuptools', 'pybind11', 'pytest', 'numpy', 'pyqt5'],
    packages=["concentrations", "elongation"],
    package_dir={"concentrations": "concentrations", "elongation": "elongation"},
    scripts=['concentrations/basepairingeditor.py', 'elongation/simulationbuilder.py'],
    cmdclass=CMDCLASS,
    ext_modules=ext_modules,
    zip_safe=False,
    include_package_data=True
)
