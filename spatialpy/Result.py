import os
import re
import shutil
import subprocess
import sys
import tempfile
import types
import warnings
import uuid


import numpy
import scipy.io
import scipy.sparse

from model import *

import inspect

try:
    # This is only needed if we are running in an Ipython Notebook
    import IPython.display
except:
    pass

try:
    import h5py
except:
    raise Exception("PyURDME requires h5py.")

try:
    import dolfin
    import mshr
except:
    raise Exception("PyURDME requires FeniCS/Dolfin.")

try:
    #dolfin.parameters["linear_algebra_backend"] = "uBLAS"
except:
    #dolfin.parameters["linear_algebra_backend"] = "Eigen"

import pickle
import json
import functools

# module-level variable to for javascript export in IPython/Jupyter notebooks
__pyurdme_javascript_libraries_loaded = False
def load_pyurdme_javascript_libraries():
    global __pyurdme_javascript_libraries_loaded
    if not __pyurdme_javascript_libraries_loaded:
        __pyurdme_javascript_libraries_loaded = True
        import os.path
        import IPython.display
        with open(os.path.join(os.path.dirname(__file__),'data/three.js_templates/js/three.js')) as fd:
            bufa = fd.read()
        with open(os.path.join(os.path.dirname(__file__),'data/three.js_templates/js/render.js')) as fd:
            bufb = fd.read()
        with open(os.path.join(os.path.dirname(__file__),'data/three.js_templates/js/OrbitControls.js')) as fd:
            bufc = fd.read()
        IPython.display.display(IPython.display.HTML('<script>'+bufa+bufc+bufb+'</script>'))


def deprecated(func):
    '''This is a decorator which can be used to mark functions
     as deprecated. It will result in a warning being emitted
     when the function is used.'''

    @functools.wraps(func)
    def new_func(*args, **kwargs):
        warnings.warn_explicit(
             "Call to deprecated function {}.".format(func.__name__),
             category=DeprecationWarning,
             filename=func.func_code.co_filename,
             lineno=func.func_code.co_firstlineno + 1
         )
        return func(*args, **kwargs)
    return new_func


# Set log level to report only errors or worse
#dolfin.set_log_level(dolfin.ERROR)
import logging
logging.getLogger('FFC').setLevel(logging.ERROR)
logging.getLogger('UFL').setLevel(logging.ERROR)

class URDMEResult(dict):
    """ Result object for a URDME simulation, extends the dict object. """

    def __init__(self, model=None, filename=None, loaddata=False):
        self.model = model
        self.sol = None
        self.U = None
        self.tspan = None
        self.data_is_loaded = False
        self.sol_initialized = False
        self.filename = filename
        if filename is not None and loaddata:
            self.read_solution()
        self.stdout = None
        self.stderr = None

    def __ne__(self, other):
        return not self.__eq__(other)

    def __eq__(self, other, verbose=False):
        try:
            tspan = self.get_timespan()
            if numpy.any(tspan != other.get_timespan()):
                if verbose: print "tspan does not match"
                return False
            for t in tspan:
                for sname in self.model.listOfSpecies:
                    if numpy.any(self.get_species(sname, timepoints=t) != other.get_species(sname, timepoints=t)):
                        if verbose: print "Species {0} does not match at t={1}".format(sname, t)
                        return False
            return True
        except ValueError as e:
            if verbose: print "value error: {0}".format(e)
            return False


    def get_endtime_model(self):
        """ Return a URDME model object with the initial conditions set to the final time point of the
            result object.
        """
        if self.model is None:
            raise Exception("can not continue a result with no model")
        # create a soft copy
        model_str = pickle.dumps(self.model)
        model2 = pickle.loads(model_str)
        # set the initial conditions
        model2.u0 = numpy.zeros(self.model.u0.shape)
        for s, sname in enumerate(self.model.listOfSpecies):
            model2.u0[s,:] = self.get_species(sname, timepoints=-1)
        return model2


    def __getstate__(self):
        """ Used by pickle to get state when pickling. We need to read the contents of the
        output file since we can't pickle file objects. """

        try:
            with open(self.filename,mode='rb') as fh:
                filecontents = fh.read()
        except Exception as e:
            raise Exception(("Error pickling model. Failed to read result file:",str(e)))

        state = self.__dict__
        state["filecontents"] = filecontents

        state["v2d"] = self.get_v2d()
        state["d2v"] = self.get_d2v()

        return state


    def __setstate__(self, state):
        """ Used by pickle to set state when unpickling. """

        # If the object contains filecontents, write those to a new tmp file.
        try:
            filecontents = state.pop("filecontents",None)
            fd = tempfile.NamedTemporaryFile(delete=False, dir=os.environ.get('PYURDME_TMPDIR'))
            with open(fd.name, mode='wb') as fh:
                fh.write(filecontents)
            state["filename"] = fd.name
        except Exception as e:
            print "Error unpickling model, could not recreate the solution file."
            raise e

        for k,v in state.items():
            self.__dict__[k] = v

    def get_v2d(self):
        """ Return the vertex-to-dof mapping. """
        if not hasattr(self, 'v2d'):
            fs = self.model.mesh.get_function_space()
            self.v2d = dolfin.vertex_to_dof_map(fs)

        return self.v2d

    def get_d2v(self):
        """ Return the dof-to-vertex mapping. """
        if not hasattr(self, 'd2v'):
            fs = self.model.mesh.get_function_space()
            self.d2v = dolfin.dof_to_vertex_map(fs)

        return self.d2v

    def _reorder_dof_to_voxel(self, M, num_species=None):
        """ Reorder the colums of M from dof ordering to vertex ordering. """

        v2d = self.get_v2d()
        if len(M.shape) == 1:
            num_timepoints = 1
        else:
            num_timepoints = M.shape[0]
        num_vox = self.model.mesh.get_num_voxels()
        if num_species is None:
            num_species = self.model.get_num_species()
        num_dofs = num_vox*num_species
        C = numpy.zeros((num_timepoints, num_dofs), dtype=numpy.float64)

        for vox_ndx in range(num_vox):
            for cndx in range(num_species):
                try:
                    if len(M.shape) == 1:
                        C[:, vox_ndx*num_species+cndx] = M[v2d[vox_ndx]*num_species+cndx]
                    else:
                        C[:, vox_ndx*num_species+cndx] = M[:, v2d[vox_ndx]*num_species+cndx]
                except IndexError as e:
                    import traceback
                    #traceback.print_stack()
                    print traceback.format_exc()
                    print "C.shape: ", C.shape
                    print "M.shape: ", M.shape
                    print "num_timepoints: ", num_timepoints
                    print "vox_ndx={0},num_species={1},cndx={2}".format(vox_ndx,num_species,cndx)
                    print "v2d[vox_ndx]={0}".format(v2d[vox_ndx])
                    print "vox_ndx*num_species+cndx={0}".format(vox_ndx*num_species+cndx)
                    print "v2d[vox_ndx]*num_species+cndx={0}".format(v2d[vox_ndx]*num_species+cndx)
                    raise e
        return C

    def read_solution(self):
        """ Read the tspan and U matrix into memory. """

        resultfile = h5py.File(self.filename, 'r')
        U = resultfile['U']
        U = numpy.array(U)

        tspan = resultfile['tspan']
        tspan = numpy.array(tspan).flatten()
        resultfile.close()

        # Reorder the dof from dof ordering to voxel ordering
        U = self._reorder_dof_to_voxel(U)

        self.U = U
        self.tspan = tspan
        self.data_is_loaded = True

    def get_timespan(self):
        if self.tspan is not None:
            resultfile = h5py.File(self.filename, 'r')
            tspan = resultfile['tspan']
            tspan = numpy.array(tspan).flatten()
            resultfile.close()
            self.tspan = tspan
        return self.tspan

    def get_species(self, species, timepoints="all", concentration=False):
        """ Returns a slice (view) of the output matrix U that contains one species for the timepoints
            specified by the time index array. The default is to return all timepoints.
            Data is loaded by slicing directly in the hdf5 dataset, i.e. it the entire
            content of the file is not loaded in memory and the U matrix
            is never added to the object.
            if concentration is False (default), the integer, raw, trajectory data is returned,
            if set to True, the concentration (=copy_number/volume) is returned.
        """

        if isinstance(species, Species):
            spec_name = species.name
        else:
            spec_name = species

        species_map = self.model.get_species_map()
        num_species = self.model.get_num_species()

        spec_indx = species_map[spec_name]

        resultfile = h5py.File(self.filename, 'r')
        #Ncells = self.model.mesh.num_vertices()  # Need dof ordering numVoxels
        U = resultfile['U']
        Ncells = U.shape[1]/num_species

        if timepoints  ==  "all":
            Uslice= U[:,(spec_indx*Ncells):(spec_indx*Ncells+Ncells)]
        else:
            Uslice = U[timepoints,(spec_indx*Ncells):(spec_indx*Ncells+Ncells)]

        if concentration:
            Uslice = self._copynumber_to_concentration(Uslice)

        # Reorder the dof from dof ordering to voxel ordering
        Uslice = self._reorder_dof_to_voxel(Uslice, num_species=1)

        # Make sure we return 1D slices as flat arrays
        dims = numpy.shape(Uslice)
        if dims[0] == 1:
            Uslice = Uslice.flatten()

        resultfile.close()
        return Uslice


    def __setattr__(self, k, v):
        if k in self.keys():
            self[k] = v
        elif not hasattr(self, k):
            self[k] = v
        else:
            raise AttributeError, "Cannot set '%s', cls attribute already exists" % ( k, )

    def __setupitems__(self, k):
        if k == 'sol' and not self.sol_initialized:
            self._initialize_sol()
        elif (k == 'U' or k == 'tspan') and not self.data_is_loaded:
            if self.filename is None:
                raise AttributeError("This result object has no data file.")
            self.read_solution()

    def __getitem__(self, k):
        self.__setupitems__(k)
        if k in self.keys():
            return self.get(k)
        raise KeyError("Object has no attribute {0}".format(k))

    def __getattr__(self, k):
        self.__setupitems__(k)
        if k in self.keys():
            return self.get(k)
        raise AttributeError("Object has no attribute {0}".format(k))

    def __del__(self):
        """ Deconstructor. """
            #   if not self.data_is_loaded:
        try:
            # Clean up data file
            os.remove(self.filename)
        except OSError as e:
            #print "URDMEResult.__del__: Could not delete result file'{0}': {1}".format(self.filename, e)
            pass

    @deprecated
    def _initialize_sol(self):
        """ Initialize the sol variable. This is a helper function for export_to_vtk(). """

        # Create Dolfin Functions for all the species
        sol = {}

        if self.model is None:
            raise URDMEError("URDMEResult.model must be set before the sol attribute can be accessed.")
        numvox = self.model.mesh.num_vertices()
        fs = self.model.mesh.get_function_space()
        vertex_to_dof_map = self.get_v2d()
        dof_to_vertex_map = self.get_d2v()

        # The result is loaded in dolfin Functions, one for each species and time point
        for i, spec in enumerate(self.model.listOfSpecies):

            species = self.model.listOfSpecies[spec]
            spec_name = species.name

            spec_sol = {}
            for j, time in enumerate(self.tspan):

                func = dolfin.Function(fs)
                func_vector = func.vector()

                S = self.get_species(spec, [j])

                for voxel in range(numvox):
                    ix  = vertex_to_dof_map[voxel]
                    try:
                        func_vector[ix] = float(S[voxel]/self.model.dofvol[ix])

                    except IndexError as e:
                        print "func_vector.size(): ", func_vector.size()
                        print "dolfvox: ",dolfvox
                        print "S.shape: ",S.shape
                        print "voxel: ",voxel
                        print "vertex_to_dof_map[voxel]", vertex_to_dof_map[voxel]
                        print "self.model.dofvol.shape: ", self.model.dofvol.shape
                        raise e

                spec_sol[time] = func

            sol[spec] = spec_sol
        self.sol = sol
        self.sol_initialized = True
        return sol

    def export_to_csv(self, folder_name):
        """ Dump trajectory to a set CSV files, the first specifies the mesh (mesh.csv) and the rest specify trajectory data for each species (species_S.csv for species named 'S').
            The columns of mesh.csv are: 'Voxel ID', 'X', 'Y', 'Z', 'Volume', 'Subdomain'.
            The columns of species_S.csv are: 'Time', 'Voxel 0', Voxel 1', ... 'Voxel N'.
        """
        import csv
        subprocess.call(["mkdir", "-p", folder_name])
        #['Voxel ID', 'X', 'Y', 'Z', 'Volume', 'Subdomain']
        with open(os.path.join(folder_name,'mesh.csv'), 'w+') as csvfile:
            writer = csv.writer(csvfile, delimiter=',')
            writer.writerow(['Voxel ID', 'X', 'Y', 'Z', 'Volume', 'Subdomain'])
            vol = self.model.get_solver_datastructure()['vol']
            for ndx in range(self.model.mesh.get_num_voxels()):
                row = [ndx]+self.model.mesh.coordinates()[ndx,:].tolist()+[vol[ndx]]+[self.model.sd[ndx]]
                writer.writerow(row)

        for spec in self.model.listOfSpecies:
            #['Time', 'Voxel 0', Voxel 1', ... 'Voxel N']
            with open(os.path.join(folder_name,'species_{0}.csv'.format(spec)), 'w+') as csvfile:
                data = self.get_species(spec)
                (num_t,num_vox) = data.shape
                writer = csv.writer(csvfile, delimiter=',')
                row = ['Time']
                for v in range(num_vox):
                    row.append('Voxel {0}'.format(v))
                writer.writerow(row)
                timespan = self.get_timespan()
                for t in range(num_t):
                    writer.writerow([timespan[t].tolist()] + data[t,:].tolist())

    def export_to_vtk(self, species, folder_name):
        """ Dump the trajectory to a collection of vtk files in the folder folder_name (created if non-existant).
            The exported data is #molecules/volume, where the volume unit is implicit from the mesh dimension. """

        #self._initialize_sol()
        subprocess.call(["mkdir", "-p", folder_name])
        fd = dolfin.File(os.path.join(folder_name, "trajectory.xdmf").encode('ascii', 'ignore'))
        func = dolfin.Function(self.model.mesh.get_function_space())
        func_vector = func.vector()
        vertex_to_dof_map = self.get_v2d()

        for i, time in enumerate(self.tspan):
            solvector = self.get_species(species,i,concentration=True)
            for j, val in enumerate(solvector):
                func_vector[vertex_to_dof_map[j]] = val
            fd << func

    def export_to_xyx(self, filename, species=None, file_format="VMD"):
        """ Dump the solution attached to a model as a xyz file. This format can be
            read by e.g. VMD, Jmol and Paraview. """

        if self.U is None:
            raise URDMEError("No solution found in the model.")

        #outfile = open(filename,"w")
        dims = numpy.shape(self.U)
        Ndofs = dims[0]
        Mspecies = len(self.model.listOfSpecies)
        Ncells = Ndofs / Mspecies

        coordinates = self.model.mesh.get_voxels()
        coordinatestr = coordinates.astype(str)

        if species == None:
            species = list(self.model.listOfSpecies.keys())

        if file_format == "VMD":
            outfile = open(filename, "w")
            filestr = ""
            for i, time in enumerate(self.tspan):
                number_of_atoms = numpy.sum(self.U[:, i])
                filestr += (str(number_of_atoms) + "\n" + "timestep " + str(i) + " time " + str(time) + "\n")
                for j, spec in enumerate(species):
                    for k in range(Ncells):
                        for mol in range(self.U[k * Mspecies + j, i]):
                            # Sample a random position in a sphere of radius computed from the voxel volume
                            # TODO: Sample volume
                            linestr = spec + "\t" + '\t'.join(coordinatestr[k, :]) + "\n"
                            filestr += linestr

            outfile.write(filestr)
            outfile.close()

        elif file_format == "ParaView":
            foldername = filename
            os.mkdir(foldername)
            for i, time in enumerate(self.tspan):
                outfile = open(foldername + "/" + filename + "." + str(i), "w")
                number_of_atoms = numpy.sum(self.U[:, i])
                filestr = ""
                filestr += (str(number_of_atoms) + "\n" + "timestep " + str(i) + " time " + str(time) + "\n")
                for j, spec in enumerate(self.model.listOfSpecies):
                    for k in range(Ncells):
                        for mol in range(model.U[k * Mspecies + j, i]):
                            linestr = spec + "\t" + '\t'.join(coordinatestr[k, :]) + "\n"
                            filestr += linestr
                outfile.write(filestr)
                outfile.close()



    def _export_to_particle_js(self,species,time_index, colors=None):
        """ Create a html string for displaying the particles as small spheres. """
        import random
        with open(os.path.dirname(os.path.abspath(__file__))+"/data/three.js_templates/particles.html",'r') as fd:
            template = fd.read()

        factor, coordinates = self.model.mesh.get_scaled_normalized_coordinates()
        dims = numpy.shape(coordinates)
        if dims[1]==2:
            is3d = 0
            vtxx = numpy.zeros((dims[0],3))
            for i, v in enumerate(coordinates):
                vtxx[i,:]=(list(v)+[0])
            coordinates = vtxx
        else:
            is3d = 1

        h = self.model.mesh.get_mesh_size()

        x=[]
        y=[]
        z=[]
        c=[]
        radius = []

        total_num_particles = 0
        #colors = ["blue","red","yellow", "green"]
        if colors == None:
            colors =  get_N_HexCol(len(species))

        if not isinstance(species, list):
           species = [species]
        for j,spec in enumerate(species):
            timeslice = self.get_species(spec, time_index)
            ns = numpy.sum(timeslice)
            total_num_particles += ns

            for i, particles in enumerate(timeslice):
                # "Radius" of voxel
                hix = h[i]*factor
                hiy = hix
                hiz = hix*is3d

                for particle in range(int(particles)):
                    x.append((coordinates[i,0]+random.uniform(-1,1)*hix))
                    y.append((coordinates[i,1]+random.uniform(-1,1)*hiy))
                    z.append((coordinates[i,2]+random.uniform(-1,1)*hiz))
                    if self.model.listOfSpecies[spec].reaction_radius:
                        radius.append(factor*self.model.listOfSpecies[spec].reaction_radius)
                    else:
                        radius.append(0.01)

                    c.append(colors[j])

        template = template.replace("__X__",str(x))
        template = template.replace("__Y__",str(y))
        template = template.replace("__Z__",str(z))
        template = template.replace("__COLOR__",str(c))
        template = template.replace("__RADIUS__",str(radius))

        return template


    def export_to_three_js(self, species, time_index):
        """ Return a json serialized document that can
            be read and visualized by three.js.
        """

        colors = self._compute_solution_colors(species,time_index)
        return self.model.mesh.export_to_three_js(colors=colors)

    def _copynumber_to_concentration(self,copy_number_data):
        """ Scale comy numbers to concentrations (in unit mol/volume),
            where the volume unit is defined by the user input.
            Dof-ordering is assumed in both solution and volumes.
        """

        shape = numpy.shape(copy_number_data)
        if len(shape) == 1:
            shape = (1,shape[0])

        scaled_sol = numpy.zeros(shape)
        scaled_sol[:,:] = copy_number_data
        dims = numpy.shape(scaled_sol)

        for t in range(dims[0]):
            timeslice = scaled_sol[t,:]
            for i,cn in enumerate(timeslice):
                scaled_sol[t, i] = float(cn)/(6.022e23*self.model.dofvol[i])

        return scaled_sol


    def _compute_solution_colors(self,species, time_index):
        """ Create a color list for species at time. """

        timeslice = self.get_species(species,time_index, concentration = True)
        colors = _compute_colors(timeslice)
        return colors

    def display_particles(self,species, time_index, width=500):
        load_pyurdme_javascript_libraries()
        hstr = self._export_to_particle_js(species, time_index)
        displayareaid=str(uuid.uuid4())
        hstr = hstr.replace('###DISPLAYAREAID###',displayareaid)
        hstr = hstr.replace('###WIDTH###',str(width))
        height = int(width*0.75)

        html = '<div style="width: {0}px; height: {1}px;" id="{2}" ></div>'.format(width, height, displayareaid)
        IPython.display.display(IPython.display.HTML(html+hstr))


    def display(self, species, time_index, opacity=1.0, wireframe=True, width=500, camera=[0,0,1]):
        """ Plot the trajectory as a PDE style plot. """
        load_pyurdme_javascript_libraries()
        data = self.get_species(species,time_index,concentration=True)
        fun = DolfinFunctionWrapper(self.model.mesh.get_function_space())
        vec = fun.vector()
        (nd,) = numpy.shape(data)
        if nd == len(vec):
            for i in range(nd):
                vec[i]=data[i]
        else:
            #v2d= self.get_v2d()
            for i in range(len(vec)):
                vec[i] = data[i] # shouldn't we use v2d or d2v here?  But it doesn't work if I do.
        fun.display(opacity=opacity, wireframe=wireframe, width=width, camera=camera)


