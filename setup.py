# Copyright 2010 Leonardo Uieda
#
# This file is part of Geothermics.
#
# Fatiando a Terra is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# Geothermics is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with Geothermics.  If not, see <http://www.gnu.org/licenses/>.

import os

from numpy.distutils.extension import Extension
from numpy.distutils.core import setup

# Define the paths
fortran_dir = os.path.join('fortran')

# Define the extention modules
diffusionfd = Extension('geothermics._diffusionfd',
                        sources=[os.path.join(fortran_dir,
                                              'diffusionfd.f95')])

ext_modules = []
ext_modules.append(diffusionfd)


if __name__ == '__main__':

    setup(name='geothermics',
          ext_modules=ext_modules
         )