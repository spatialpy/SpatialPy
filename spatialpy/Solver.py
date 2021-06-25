import os
import shutil
import signal
import subprocess
import sys
import tempfile
import threading
import time
import tempfile
import getpass
import re


from spatialpy.Model import *
from spatialpy.Result import *

def read_from_stdout(stdout,verbose=True):
    ''' Used with subprocess.Popen and threading to capture all output and print
        to the screen notebook without waiting or storing the output ina buffer.'''
    try:
        while True:
            line = stdout.readline()
            if line != b'':
                if verbose:
                    print(line.decode(),end='')
            else:
                #got empty line, ending
                return
    except Exception as e:
        print("read_from_stdout(): {0}".format(e))


class Solver:
    """ SpatialPy solver object."""

    def __init__(self, model, debug_level=0):
        """ Constructor. """
        # TODO: fix class checking
        # if not isinstance(model, Model):
        #    raise SimulationError("Solver constructors must take a Model as an argument.")
        # if not issubclass(self.__class__, Solver):
        #    raise SimulationError("Solver classes must be a subclass of SpatialPy.Solver.")

        self.model = model
        self.is_compiled = False
        self.debug_level = debug_level
        self.model_name = self.model.name
        self.build_dir = None
        self.executable_name = 'ssa_sdpd'
        self.h = None  # basis function width

        self.SpatialPy_ROOT = os.path.dirname(
            os.path.abspath(__file__))+"/ssa_sdpd-c-simulation-engine"
        self.SpatialPy_ROOTDIR =  self.SpatialPy_ROOT.replace(" ","\\ ");
        self.SpatialPy_ROOTINC =  self.SpatialPy_ROOT.replace(" ","\\\\ ");
        self.SpatialPy_ROOTPARAM =  self.SpatialPy_ROOT.replace(" ","?");
        #print("SpatialPy_ROOTDIR = "+self.SpatialPy_ROOTDIR)
        #print("SpatialPy_ROOTPARAM = "+self.SpatialPy_ROOTPARAM)

        tmpdir = tempfile.gettempdir()
        self.core_dir = os.path.join(os.path.join(tmpdir, 'spatialpy_core'), getpass.getuser())

        if not os.path.isdir(os.path.join(tmpdir, 'spatialpy_core')):
            os.mkdir(os.path.join(tmpdir, 'spatialpy_core'))
        if not os.path.isdir(self.core_dir):
            os.mkdir(self.core_dir)

    def __del__(self):
        """ Deconstructor.  Removes the compiled solver."""
        try:
            if self.build_dir is not None:
                try:
                    shutil.rmtree(self.build_dir)
                except OSError as e:
                    print("Could not delete '{0}'".format(
                        self.solver_base_dir))
        except Exception as e:
            pass

    def run_debugger(self):
        """ Start a gdbgui debugger at port 5000. """
        self.debugger_url = 'http://127.0.0.1:5000'
        if not hasattr(self, 'debugger_process'):
            self.debugger_process = subprocess.Popen('gdbgui -r', shell=True)
        print(f'Your debugger is running at {self.debugger_url}')

    def compile(self, debug=False, profile=False):
        """ Compile the model."""

        # Create a unique directory each time call to compile.
        self.build_dir = tempfile.mkdtemp(
            prefix='spatialpy_build_', dir=os.environ.get('SPATIALPY_TMPDIR'))

        if self.debug_level >= 1:
            print("Compiling Solver.  Build dir: {0}".format(self.build_dir))

        # Write the propensity file
        self.propfilename = re.sub('[^\w\_]', '', self.model_name) # Match except word characters \w = ([a-zA-Z0-9_]) and \_ = _ replace with ''
        self.propfilename = self.propfilename + '_generated_model'
        self.prop_file_name = self.build_dir + '/' + self.propfilename + '.c'
        if self.debug_level > 1:
            print("Creating propensity file {0}".format(self.prop_file_name))
        self.create_propensity_file(file_name=self.prop_file_name)

        # Build the solver
        makefile = self.SpatialPy_ROOTDIR+'/build/Makefile'
        cmd_list = ['cd', self.core_dir, '&&', 'make', 'CORE', '-f',  makefile,
            'ROOT="' + self.SpatialPy_ROOTPARAM+'"',
            'ROOTINC="' + self.SpatialPy_ROOTINC+'"',
            'BUILD='+self.core_dir, '&&'
            'cd', self.build_dir, '&&', 'make', '-I', self.core_dir, '-f', makefile,
            'ROOT="' + self.SpatialPy_ROOTPARAM+'"',
            'ROOTINC="' + self.SpatialPy_ROOTINC+'"',
            'COREDIR="' + self.core_dir + '"',
            'MODEL=' + self.prop_file_name, 'BUILD='+self.build_dir]
        if profile:
          cmd_list.append('GPROFFLAG=-pg')
        if profile or debug:
          cmd_list.append('GDB_FLAG=-g')
        cmd = " ".join(cmd_list)
        if self.debug_level > 1:
            print("cmd: {0}\n".format(cmd))
        try:
            handle = subprocess.Popen(cmd, stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT, shell=True)
            stdout, _ = handle.communicate()
            return_code = handle.wait()
            if return_code != 0:
                try:
                    print(stdout.decode("utf-8"))
                except Exception as e:
                    pass
                raise SimulationError(
                    "Compilation of solver failed, return_code={0}".format(
                        return_code))

            if self.debug_level > 1:
                print(stdout.decode("utf-8"))

        except OSError as e:
            print("Error, execution of compilation raised an exception: {0}".format(
                e))
            print("cmd = {0}".format(cmd))
            raise SimulationError("Compilation of solver failed")

        self.is_compiled = True


    def run(self, number_of_trajectories=1, seed=None, timeout=None,
                number_of_threads=None, debug=False, profile=False, verbose=True):
        """ Run one simulation of the model.
        Args:
            number_of_trajectories: (int) How many trajectories should be simulated.
            seed: (int) the random number seed (incremented by one for multiple runs).
            timeout: (int) maximum number of seconds the solver can run.
            number_of_threads: (int) the number threads the solver will use.
            debug: (bool) start a gdbgui debugger (also compiles with debug symbols
                if compilation hasn't happened)
            profile: (bool) output gprof profiling data if available
        Returns:
            Result object.
                or, if number_of_trajectories > 1
            a list of Result objects
        """
        if number_of_trajectories > 1:
            result_list = []
        # Check if compiled, call compile() if not.
        if not self.is_compiled:
            self.compile(debug=debug, profile=profile)

        # Execute the solver
        for run_ndx in range(number_of_trajectories):
            outfile = tempfile.mkdtemp(
                prefix='spatialpy_result_', dir=os.environ.get('SPATIALPY_TMPDIR'))
            result = Result(self.model, outfile)
            solver_cmd = 'cd {0}'.format(
                outfile) + ";" + os.path.join(self.build_dir, self.executable_name)

            if number_of_threads is not None:
                solver_cmd += " -t " + str(number_of_threads)

            if seed is not None:
                solver_cmd += " -s "+str(seed+run_ndx)

            if self.debug_level > 1:
                print('cmd: {0}\n'.format(solver_cmd))

            start = time.monotonic()
            return_code = None
            try:
                with subprocess.Popen(solver_cmd, shell=True,
                        stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                        start_new_session=True) as process:
                    try:
                        # start thread to read process stdout to stdout
                        t = threading.Thread(target=read_from_stdout,
                            args=(process.stdout,verbose,))
                        t.start()
                        if timeout is not None:
                            return_code = process.wait(timeout=timeout)
                        else:
                            return_code = process.wait()
                        t.join()
                    except KeyboardInterrupt:
                        # send signal to the process group
                        os.killpg(process.pid, signal.SIGINT)
                        print('Terminated by user after seconds: {:.2f}'.format(
                            time.monotonic() - start))
                    except subprocess.TimeoutExpired:
                        result.timeout = True
                        # send signal to the process group
                        os.killpg(process.pid, signal.SIGINT)
                        raise SimulationTimeout("SpatialPy solver timeout exceded.")
            except OSError as e:
                print("Error, execution of solver raised an exception: {0}".format(
                    e))
                print("cmd = {0}".format(solver_cmd))

            if self.debug_level >= 1:  # output time
                print('Elapsed seconds: {:.2f}'.format(
                    time.monotonic() - start))

            if return_code is not None and return_code != 0:
                print("solver_cmd = {0}".format(solver_cmd))
                raise SimulationError(
                    "Solver execution failed, return code = {0}".format(return_code))

            result.success = True
            if profile:
                self.read_profile_info(result)
            if number_of_trajectories > 1:
                result_list.append(result)
            else:
                return result

        return result_list

    def read_profile_info(self, result):
        profile_data_path = os.path.join(result.result_dir, 'gmon.out')
        exe_path = os.path.join(self.build_dir, self.executable_name)
        cmd = f'gprof {exe_path} {profile_data_path}'
        print(cmd)
        proc = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE)
        stdout, stderr = proc.communicate()
        proc.wait()
        print(f'Gprof report for {result.result_dir}')
        print(stdout.decode('utf-8'))

    def create_propensity_file(self, file_name=None):
        """ Generate the C propensity file that is used to compile the solvers.
        """

        num_types = len(self.model.listOfTypeIDs)
        if self.model.enable_pde:
            num_chem_species = len(self.model.listOfSpecies)
            num_chem_rxns = len(self.model.listOfReactions)
        else:
            num_chem_species = 0
            num_chem_rxns = 0
        if self.model.enable_rdme:
            num_stoch_species = len(self.model.listOfSpecies)
            num_stoch_rxns = len(self.model.listOfReactions)
        else:
            num_stoch_species = 0
            num_stoch_rxns = 0
        num_data_fn = len(self.model.listOfDataFunctions)


        template = open(os.path.abspath(os.path.dirname(
            __file__)) + '/ssa_sdpd-c-simulation-engine/propensity_file_template.cpp', 'r')
        propfile = open(file_name, "w")
        propfilestr = template.read()

        speciesdef = ""
        i = 0
        for S in self.model.listOfSpecies:
            speciesdef += "#define " + S + " " + "x[" + str(i) + "]" + "\n"
            i += 1

        propfilestr = propfilestr.replace("__DEFINE_SPECIES__", speciesdef)

        propfilestr = propfilestr.replace(
            "__NUMBER_OF_REACTIONS__", str(self.model.get_num_reactions()))
        propfilestr = propfilestr.replace(
            "__NUMBER_OF_SPECIES__", str(self.model.get_num_species()))
        propfilestr = propfilestr.replace(
            "__NUMBER_OF_VOXELS__", str(self.model.domain.get_num_voxels()))

        # Make sure all paramters are evaluated to scalars before we write them to the file.
        self.model.resolve_parameters()
        parameters = ""
        for p in self.model.listOfParameters:
            parameters += "const double " + p + " = " + \
                str(self.model.listOfParameters[p].value) + ";\n"
        propfilestr = propfilestr.replace(
            "__DEFINE_PARAMETERS__", str(parameters))

        # Reactions
        funheader = "double __NAME__(const int *x, double t, const double vol, const double *data_fn, int sd)"

        funcs = ""
        funcinits = ""
        i = 0
        for R in self.model.listOfReactions:
            func = ""
            rname = self.model.listOfReactions[R].name
            func += funheader.replace("__NAME__", rname) + "\n{\n"
            if self.model.listOfReactions[R].restrict_to == None or (isinstance(self.model.listOfReactions[R].restrict_to, list) and len(self.model.listOfReactions[R].restrict_to) == 0):
                func += "return "
                func += self.model.listOfReactions[R].propensity_function
                func += ";"
            else:
                func += "if("
                if isinstance(self.model.listOfReactions[R].restrict_to, list) and len(self.model.listOfReactions[R].restrict_to) > 0:
                    for sd in self.model.listOfReactions[R].restrict_to:
                        func += "sd == " + str(sd) + "||"
                    func = func[:-2]
                elif isinstance(self.model.listOfReactions[R].restrict_to, int):
                    func += "sd == " + \
                        str(self.model.listOfReactions[R].restrict_to)
                else:
                    raise SimulationError(
                        "When restricting reaction to types, you must specify either a list or an int")
                func += "){\n"
                func += "return "
                func += self.model.listOfReactions[R].propensity_function
                func += ";"
                func += "\n}else{"
                func += "\n\treturn 0.0;}"

            func += "\n}"
            funcs += func + "\n\n"
            funcinits += "    ptr[" + \
                str(i) + "] = (PropensityFun) " + rname + ";\n"
            i += 1

        propfilestr = propfilestr.replace("__DEFINE_REACTIONS__", funcs)
        propfilestr = propfilestr.replace("__DEFINE_PROPFUNS__", funcinits)

        # TODO make deterministic chemical reaction functions work
        funheader = "double det__NAME__(const double *x, double t, const double vol, const double *data_fn, int sd)"
        deterministic_chem_rxn_functions = ""
        deterministic_chem_rxn_function_init = ""
        ############################################
        i = 0
        for R in self.model.listOfReactions:
            func = ""
            rname = self.model.listOfReactions[R].name
            func += funheader.replace("__NAME__", rname) + "\n{\n"
            if self.model.listOfReactions[R].restrict_to == None or (isinstance(self.model.listOfReactions[R].restrict_to, list) and len(self.model.listOfReactions[R].restrict_to) == 0):
                func += "return "
                func += self.model.listOfReactions[R].ode_propensity_function
                func += ";"
            else:
                func += "if("
                if isinstance(self.model.listOfReactions[R].restrict_to, list) and len(self.model.listOfReactions[R].restrict_to) > 0:
                    for sd in self.model.listOfReactions[R].restrict_to:
                        func += "sd == " + str(sd) + "||"
                    func = func[:-2]
                elif isinstance(self.model.listOfReactions[R].restrict_to, int):
                    func += "sd == " + \
                        str(self.model.listOfReactions[R].restrict_to)
                else:
                    raise SimulationError(
                        "When restricting reaction to types, you must specify either a list or an int")
                func += "){\n"
                func += "return "
                func += self.model.listOfReactions[R].ode_propensity_function
                func += ";"
                func += "\n}else{"
                func += "\n\treturn 0.0;}"

            func += "\n}"
            deterministic_chem_rxn_functions += func + "\n\n"
            deterministic_chem_rxn_function_init += "    ptr[" + \
                str(i) + "] = (ChemRxnFun) det" + rname + ";\n"
            i += 1
        ############################################
        propfilestr = propfilestr.replace("__DEFINE_CHEM_FUNS__", deterministic_chem_rxn_functions)
        propfilestr = propfilestr.replace("__DEFINE_CHEM_FUN_INITS__", deterministic_chem_rxn_function_init)

        # End of pyurdme replacements
        # SSA-SDPD values here
        init_particles = ""
        if self.model.domain.type is None:
            self.model.domain.type = numpy.ones(self.model.domain.get_num_voxels())
        for i in range(len(self.model.domain.type)):
            if self.model.domain.type[i] == 0:
                raise SimulationError(
                    "Not all particles have been defined in a type. Mass and other properties must be defined")
            init_particles += "init_create_particle(sys,id++,{0},{1},{2},{3},{4},{5},{6},{7},{8});".format(
                self.model.domain.coordinates()[i,0],self.model.domain.coordinates()[i,1],self.model.domain.coordinates()[i,2],
                self.model.domain.type[i],self.model.domain.nu[i],self.model.domain.mass[i],
                (self.model.domain.mass[i] / self.model.domain.vol[i]),int(self.model.domain.fixed[i]),num_chem_species )+"\n"
        propfilestr = propfilestr.replace("__INIT_PARTICLES__", init_particles)

        # process initial conditions here
        self.model.apply_initial_conditions()
        nspecies = self.model.u0.shape[0]
        ncells = self.model.u0.shape[1]

        input_constants = ""

        outstr = "unsigned int input_u0[{0}] = ".format(nspecies*ncells)
        outstr += "{"
        if len(self.model.listOfSpecies) > 0:
            for i in range(ncells):
                for s in range(nspecies):
                    if i+s > 0:
                        outstr += ','
                    outstr += str(int(self.model.u0[s, i]))
        outstr += "};"
        input_constants += outstr + "\n"

        data_fn_defs = ""
        if len(self.model.listOfSpecies) > 0:
            for ndf in range(len(self.model.listOfDataFunctions)):
                data_fn_defs += "#define {0} data_fn[{1}]\n".format(
                    self.model.listOfDataFunctions[ndf].name, ndf)
        propfilestr = propfilestr.replace(
            "__DATA_FUNCTION_DEFINITIONS__", data_fn_defs)

        if len(self.model.listOfSpecies) > 0:
            N = self.model.create_stoichiometric_matrix()
            if(min(N.shape) > 0):
                Nd = N.todense() # this will not work if Nrxn or Nspecies is zero
                outstr = "static int input_N_dense[{0}] = ".format(
                    Nd.shape[0] * Nd.shape[1])
                outstr += "{"
                for i in range(Nd.shape[0]):
                    for j in range(Nd.shape[1]):
                        if j+i > 0:
                            outstr += ','
                        outstr += "{0}".format(int(Nd[i, j]))
                outstr += "};\n"
                outstr += "static size_t input_irN[{0}] = ".format(len(N.indices))
                outstr += "{"
                for i in range(len(N.indices)):
                    if i > 0:
                        outstr += ','
                    outstr += str(N.indices[i])
                outstr += "};"
                input_constants += outstr + "\n"
                outstr = "static size_t input_jcN[{0}] = ".format(len(N.indptr))
                outstr += "{"
                for i in range(len(N.indptr)):
                    if i > 0:
                        outstr += ','
                    outstr += str(N.indptr[i])
                outstr += "};"
                input_constants += outstr + "\n"
                outstr = "static int input_prN[{0}] = ".format(len(N.data))
                outstr += "{"
                for i in range(len((N.data))):
                    if i > 0:
                        outstr += ','
                    outstr += str(int(N.data[i]))
                outstr += "};"
                input_constants += outstr + "\n"
            else:
                input_constants += "static int input_N_dense[0] = {};\n"
                input_constants += "static size_t input_irN[0] = {};\n"
                input_constants += "static size_t input_jcN[0] = {};\n"
                input_constants += "static int input_prN[0] = {};\n"

            G = self.model.create_dependency_graph()
            outstr = "static size_t input_irG[{0}] = ".format(len(G.indices))
            outstr += "{"
            for i in range(len(G.indices)):
                if i > 0:
                    outstr += ','
                outstr += str(G.indices[i])
            outstr += "};"
            input_constants += outstr + "\n"

            outstr = "static size_t input_jcG[{0}] = ".format(len(G.indptr))
            outstr += "{"
            for i in range(len(G.indptr)):
                if i > 0:
                    outstr += ','
                outstr += str(G.indptr[i])
            outstr += "};"
            input_constants += outstr + "\n"
        if(len(self.model.listOfSpecies)>0):
            outstr = "const char* const input_species_names[] = {"
            for i, s in enumerate(self.model.listOfSpecies.keys()):
                if i > 0:
                    outstr += ","
                outstr += '"'+s+'"'
            outstr += ", 0};"
            input_constants += outstr + "\n"
        num_types = len(self.model.listOfTypeIDs)
        outstr = "const int input_num_subdomain = {0};".format(num_types)
        input_constants += outstr + "\n"
        outstr = "const double input_subdomain_diffusion_matrix[{0}] = ".format(
            len(self.model.listOfSpecies)*num_types)
        outstr += "{"
        for i, sname in enumerate(self.model.listOfSpecies.keys()):
            s = self.model.listOfSpecies[sname]
            for j, sd_id in enumerate(self.model.listOfTypeIDs):
                if i+j > 0:
                    outstr += ','
                try:
                    if s not in self.model.listOfDiffusionRestrictions or \
                       sd_id in self.model.listOfDiffusionRestrictions[s]:
                        outstr += "{0}".format(s.diffusion_constant)
                    else:
                        outstr += "0.0"
                except KeyError as e:
                    print("error: {0}".format(e))
                    print(self.model.listOfDiffusionRestrictions)
                    raise Exception()

        outstr += "};"
        input_constants += outstr + "\n"
        propfilestr = propfilestr.replace(
            "__INPUT_CONSTANTS__", input_constants)

        system_config = "debug_flag = {0};\n".format(self.debug_level)
        system_config += "ParticleSystem *system = new ParticleSystem({0},{1},{2},{3},{4},{5});\n".format(
            num_types, num_chem_species, num_chem_rxns, 
            num_stoch_species, num_stoch_rxns, num_data_fn
            )
        system_config += "system->static_domain = {0};\n".format(int(self.model.staticDomain))
        if(len(self.model.listOfSpecies) > 0):
            system_config += "system->subdomain_diffusion_matrix = input_subdomain_diffusion_matrix;\n"
            system_config += "system->stoichiometric_matrix = input_N_dense;\n"
            system_config += "system->chem_rxn_rhs_functions = ALLOC_ChemRxnFun();\n"
            system_config += "system->stoch_rxn_propensity_functions = ALLOC_propensities();\n"
            system_config += "system->species_names = input_species_names;\n";

        system_config += "system->dt = {0};\n".format(self.model.timestep_size)
        system_config += "system->nt = {0};\n".format(self.model.num_timesteps)
        system_config += "system->output_freq = {0};\n".format(self.model.output_freq)
        if self.h is None:
            self.h = self.model.domain.find_h()
        if self.h == 0.0:
            raise ModelError('h (basis function width) can not be zero.')
        system_config +="system->h = {0};\n".format(self.h)
        system_config +="system->rho0 = {0};\n".format(self.model.domain.rho0)
        system_config +="system->c0 = {0};\n".format(self.model.domain.c0)
        system_config +="system->P0 = {0};\n".format(self.model.domain.P0)
        #// bounding box
        system_config += "system->xlo = {0};\n".format(self.model.domain.xlim[0])
        system_config += "system->xhi = {0};\n".format(self.model.domain.xlim[1])
        system_config += "system->ylo = {0};\n".format(self.model.domain.ylim[0])
        system_config += "system->yhi = {0};\n".format(self.model.domain.ylim[1])
        system_config += "system->zlo = {0};\n".format(self.model.domain.zlim[0])
        system_config += "system->zhi = {0};\n".format(self.model.domain.zlim[1])
        #
        if self.model.domain.gravity is not None:
            for i in range(3):
                system_config += "system->gravity[{0}] = {1};\n".format(i,self.model.domain.gravity[i])


        propfilestr = propfilestr.replace("__SYSTEM_CONFIG__", system_config)

        init_rdme=''
        if(self.model.enable_rdme and len(self.model.listOfSpecies) > 0):
            init_rdme = '''
            initialize_rdme(system, input_irN, input_jcN, input_prN, input_irG,
                            input_jcG, input_u0 );'''
        propfilestr = propfilestr.replace("__INIT_RDME__", init_rdme)


        init_bc = ""
        for bc in self.model.listOfBoundaryConditions:
            init_bc += bc.expression()
        for ndf, df in enumerate(self.model.listOfDataFunctions):
            init_bc += "me->data_fn[{0}] = {1};\n".format(ndf,df.expression())
        propfilestr = propfilestr.replace("__BOUNDARY_CONDITIONS__", init_bc)

        #### Write the data to the file ####
        propfile.write(propfilestr)
        propfile.close()


class SimulationError(Exception):
    pass


class SimulationTimeout(SimulationError):
    pass
