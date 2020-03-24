# -*- coding: utf-8 -*-
# :Project:   elongation_simulators -- Packaging
# :Author:    Fabio Hedayioglu <fheday@gmail.com>
# :License:   MIT License
# :Copyright: © 2020 Fabio Hedayioglu
#

from setuptools import setup, Extension, find_packages
import os
import sys
import subprocess
from pathlib import Path

from setuptools.command.build_ext import build_ext


class CMakeExtension(Extension):
    def __init__(self, name):
        Extension.__init__(self, name, sources=[])


class CMakeBuild(build_ext):
    def run(self):
        env = os.environ.copy()
        env['PATH'] += ":~/.local/bin/"
        try:
            out = subprocess.run(['cmake', '--version'], env=env)
        except OSError:
            raise RuntimeError(
                "CMake must be installed to build the following extensions: " +
                ", ".join(e.name for e in self.extensions))

        build_directory = os.path.abspath(self.build_temp)

        cmake_args = [
            '-DCMAKE_LIBRARY_OUTPUT_DIRECTORY=' + build_directory,
            '-DPYTHON_EXECUTABLE=' + sys.executable,
            '-DCMAKE_PREFIX_PATH='+ '~/.local/lib/python3.6/site-packages/'
        ]

        cfg = 'Debug' if self.debug else 'Release'
        build_args = ['--config', cfg]

        cmake_args += ['-DCMAKE_BUILD_TYPE=' + cfg]

        # Assuming Makefiles
        build_args += ['--', 'ribosomesimulator', 'translation', '-j4']

        self.build_args = build_args

        env['CXXFLAGS'] = '{} -DVERSION_INFO=\\"{}\\"'.format(
            env.get('CXXFLAGS', ''),
            self.distribution.get_version())
        if not os.path.exists(self.build_temp):
            os.makedirs(self.build_temp)

        # CMakeLists.txt is in the same directory as this setup.py file
        cmake_list_dir = os.path.abspath(os.path.dirname(__file__))
        print('-'*10, 'Running CMake prepare', '-'*40)
        subprocess.check_call(['cmake', cmake_list_dir] + cmake_args,
                              cwd=self.build_temp, env=env)

        print('-'*10, 'Building extensions', '-'*40)
        cmake_cmd = ['cmake', '--build', '.'] + self.build_args
        subprocess.check_call(cmake_cmd,
                              cwd=self.build_temp, env=env)

        # Move from build temp to final position
        for ext in self.extensions:
            self.move_output(ext)

    def move_output(self, ext):
        build_temp = Path(self.build_temp).resolve()
        dest_path = Path(self.get_ext_fullpath(ext.name)).resolve()
        source_path = build_temp / self.get_ext_filename(ext.name)
        dest_directory = dest_path.parents[0]
        dest_directory.mkdir(parents=True, exist_ok=True)
        self.copy_file(source_path, dest_path)

if sys.version_info < (3,):
    raise NotImplementedError("Only Python 3+ is supported.")


with open('./version.txt', encoding='utf-8') as f:
    VERSION = f.read()

with open('./README.rst', encoding='utf-8') as f:
    LONG_DESCRIPTION = f.read()

with open('./CHANGES.rst', encoding='utf-8') as f:
    CHANGES = f.read()

setup(
    name='elongation_simulators',
    version=VERSION,
    description='Python wrapper around translation and ribosome simulators',
    long_description=LONG_DESCRIPTION + '\n\n' + CHANGES,
    long_description_content_type='text/x-rst',
    license='MIT License',
    keywords='elongation translation',
    author='Fabio Hedayioglu',
    author_email='fheday@gmail.com',
    maintainer='Fabio Hedayioglu',
    maintainer_email='fheday@gmail.com',
    url='https://github.com/fheday',
    classifiers=[
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: C++',
        'Programming Language :: Python :: 3 :: Only',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.5',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python',
    ],

    packages=find_packages(),
    install_requires=['cmake', 'pybind11', 'pybind11-cmake'],
    ext_modules=[
                CMakeExtension('ribosomesimulator'),
                CMakeExtension('translation')
    ],
    cmdclass=dict(build_ext=CMakeBuild),
#  zip_safe=False,
)
