""" instaler """
# -*- coding: utf-8 -*-
# :Project:   elongation_simulators -- Packaging
# :Author:    Fabio Hedayioglu <fheday@gmail.com>
# :License:   MIT License
# :Copyright: © 2020 Fabio Hedayioglu
#

import os
import sys
import sysconfig
from setuptools import setup
from pybind11.setup_helpers import Pybind11Extension, build_ext, ParallelCompile, naive_recompile # noqa:E402
DIR = os.path.abspath(os.path.dirname(__file__))
sys.path.append(os.path.join(DIR, "pybind11"))
del sys.path[-1]


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
    EXTRA_COMPILE_ARGS = ["/O2", "/Ot", "/GL", "/DCOMIPLE_PYTHON_MODULE", "/I./eigen-3.3.7/eigen3/"]
else:
    EXTRA_COMPILE_ARGS = ["-O3", "-ffast-math", "-ftree-vectorize", "-Wall",\
                          "-g2", "-flto", "-DCOMIPLE_PYTHON_MODULE"]

ext_modules = [
    Pybind11Extension(
        "translation",
        ["concentrationsreader.cpp", "mrna_reader.cpp", "elongation_codon.cpp", "initiationterminationcodon.cpp",\
         "mrnaelement.cpp", "translation.cpp", "ribosomesimulator.cpp", "elongation_simulation_manager.cpp",\
         "elongation_simulation_processor.cpp", "./jsoncpp/jsoncpp.cpp"],
        include_dirs=["./jsoncpp/", "./eigen-3.3.7/", "./pybind11/"],
        extra_compile_args=EXTRA_COMPILE_ARGS
    ),
    Pybind11Extension(
        "ribosomesimulator",
        ["concentrationsreader.cpp", "mrna_reader.cpp", "ribosomesimulator.cpp", "./jsoncpp/jsoncpp.cpp"],
        include_dirs=["./jsoncpp/", "./eigen-3.3.7/", "./pybind11/"],
        extra_compile_args=EXTRA_COMPILE_ARGS
    )
]


if sys.version_info < (3,):
    raise NotImplementedError("Only Python 3+ is supported.")


with open('version.txt', encoding='utf-8') as f:
    VERSION = f.read()

with open('README.rst', encoding='utf-8') as f:
    LONG_DESCRIPTION = f.read()

with open('CHANGES.rst', encoding='utf-8') as f:
    CHANGES = f.read()

# Optional multithreaded build
ParallelCompile(default=0, needs_recompile=naive_recompile).install()

CMDCLASS = {}
if not WIN:
    CMDCLASS = {"build_ext": BuildExt}



setup(
    name='elongation_simulators',
    version=VERSION,
    description='Python wrapper around translation and ribosome simulators',
    long_description=LONG_DESCRIPTION + '\n\n' + CHANGES,
    #long_description_content_type='text/x-rst',
    license='MIT License',
    keywords='elongation translation',
    author='Fabio Hedayioglu',
    author_email='fheday@gmail.com',
    maintainer='Fabio Hedayioglu',
    maintainer_email='fheday@gmail.com',
    url='https://github.com/fheday/elongation_simulator/',
    classifiers=[
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: C++',
        'Programming Language :: Python :: 3 :: Only',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9',
        'Programming Language :: Python :: 3.10',
        'Programming Language :: Python :: 3.11',
        'Programming Language :: Python',
    ],
   install_requires=['pybind11', 'pytest', 'numpy'],
   packages=["concentrations", "elongation"],
   package_dir={"concentrations":"concentrations", "elongation":"elongation"},
   cmdclass=CMDCLASS,
   ext_modules=ext_modules,
   zip_safe=False,
   include_package_data=True
)
