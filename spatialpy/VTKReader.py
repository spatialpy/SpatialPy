'''
SpatialPy is a Python 3 package for simulation of
spatial deterministic/stochastic reaction-diffusion-advection problems
Copyright (C) 2021 SpatialPy developers.

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU GENERAL PUBLIC LICENSE Version 3 as
published by the Free Software Foundation.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU GENERAL PUBLIC LICENSE Version 3 for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
'''

import numpy
import math


class VTKReader:
    """VTKReader.py: SpatialPy minimal VTK legacy file reader."""
    """Reference: https://vtk.org/wp-content/uploads/2015/04/file-formats.pdf"""

    def __init__(self, filename=None, debug=False):

        self.filename = filename
        self.pointdatatype = None
        self.numpoints = None
        self.points = None
        self.arrays = None
        self.debug = debug
        self.datatypes = {
            "int": "int32",
            "float": "float32",
            "double": "float64",
        }

    def setfilename(self, filename):
        """Set filename.
        Args:
            (str) filename
        """

        self.filename = filename

    def getarrayname(self, i):
        """Get (str) array name.
        Args:
            (int) index"""

        arrayids = list(self.arrays.keys())

        if i <= len(arrayids):
            return arrayids[i]
        else:
            return None

    def getarrays(self):
        """Get (dict) of arrays."""

        return self.arrays

    def getnumpoints(self):
        """Get (int) number of points."""

        return self.numpoints

    def getpoints(self):
        """Get (list) points."""
        return self.points

    def isvalidnum(self, numstr):
        """Test if string is a valid numeric value.
        Args:
            (str) string to test
        Return:
            (boolean) True/False
        """

        try:
            float(numstr)
        except:
            return False

        return True

    def readnumeric(self, fd):
        """Read numeric data to be returned in a list.
        Args:
            file descriptor
        Return:
            (list) numeric data read"""

        numericlist = []

        for line in fd:
            #print('line={0}'.format(line))
            l = line.strip().split()
            if(len(l)==0):
                break
            if self.isvalidnum(l[0]):
                numericlist.extend(l)
            else:
                break

        return numericlist

    def readpoints(self, fd, row, col, datatype):
        """Read numeric data to be returned in a shaped numpy array.
        Args:
            file descriptor
            (int) number of rows
            (int) number of columns
            (str) VTK datatype
        Return:
            (numpy array) of shape row by col of (numpy) dtype
        """

        pointlist = numpy.array(self.readnumeric(fd), dtype=self.datatypes[datatype])
        if col > 1:
            pointlist = pointlist.reshape(row, col)

        return pointlist

    # This should be refactored at some point and probably eliminated
    # In factor of readnumeric()
    def populatearrays(self, vtkdata, arraydata, col, row, name, datatype):
        """Populate arrays with data to be added to dict of shaped arrays.
        Args:
            (dict) array data
            (array) array data to be added to dict
            (int) number of rows
            (int) number of columns
            (str) name of dict key
            (str) VTK datatype
        """

        array = numpy.array(arraydata, dtype=datatype)
        if col > 1:
            array = array.reshape(row, col)
        vtkdata[name] = array

    def readarrays(self, fd):
        """Read all arrays and data in file.
        Args:
            file descriptor
        Return:
            (dict) array data
        """
        vtkdata = {}
        arraydata = []

        for line in fd:
            #if self.debug: print("line={0}".format(line),end='')
            if line.isspace():
                continue
            try:
                name, col, row, datatype = line.strip().split()
                # got another PointData section
            except Exception as e:
                print("Error: {1}".format(e))
                print("on line >>>")
                print(line, end='')
                print("<<< {0}".format(fd.tell()))
                raise e
            col = int(col)
            row = int(row)

            # now read row*col number of values
            arraydata = []
            while len(arraydata) < row*col:
                line = fd.readline()
                arraydata.extend(line.strip().split())

            #if self.debug: print("populatearrays(name={0})".format(name))
            self.populatearrays(vtkdata, arraydata, col, row, name, datatype)

        return vtkdata

    def readfile(self):
        """Read VTK file."""

        with open(self.filename) as fd:
            if self.debug: print("open({0})".format(self.filename))
            tmp =  fd.readline()
            #if self.debug:  print("line={0}".format(tmp),end='')
            tmp =  fd.readline()
            #if self.debug:  print("line={0}".format(tmp),end='')

            # We only output ASCII so we can ignore BINARY
            tmp =  fd.readline()
            #if self.debug: print("line={0}".format(tmp),end='')
            if tmp.strip().upper() != "ASCII":
                raise VTKReaderIOError("{0} doesn't look like a valid ASCII VTK file.".format(self.filename))

            tmp =  fd.readline()
            #if self.debug:  print("line={0}".format(tmp),end='')

            tmp =  fd.readline()
            #if self.debug: print("line={0}".format(tmp),end='')
            _, self.numpoints, self.pointdatatype = tmp.strip().split()
            self.numpoints = int(self.numpoints)

            #if self.debug: print("self.readpoints(numpoints={0})".format(self.numpoints),end='')
            self.points = self.readpoints(fd, self.numpoints, 3, self.pointdatatype)
            #if self.debug: print("self.points.shape = {0}".format(self.points.shape))

            for line in fd:
                if line[:5] != "FIELD":
                    #print("skipping line: {0}".format(line),end='')
                    continue
                else:
                    #if self.debug:  print("break skipping on line: {0}".format(line),end='')
                    break

            self.arrays = self.readarrays(fd)
            #if self.debug: print("self.arrays.keys = {0}".format(self.arrays.keys()))


class VTKReaderError(Exception):
    """Base class for exceptions in VTKReader module."""

    pass


class VTKReaderIOError(VTKReaderError):
    """Exception raised for I/O errors."""

    def __init__(self, message):
        self.message = message
