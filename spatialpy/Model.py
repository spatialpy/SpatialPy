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

#This module defines a model that simulates a discrete, stoachastic, mixed biochemical reaction network in python.

import uuid
from collections import OrderedDict
from spatialpy.Solver import Solver
import numpy
import scipy
import warnings
import math


def export_StochSS(spatialpy_model, filename=None, return_stochss_model=False):
    """
    SpatialPy model to StochSS converter

    Args:
        spatialpy_model : spatialpy.Model
            SpatialPy model to be converted to StochSS
        filename : str
            Path to the exported stochss model
        return_stochss_model : bool
            Whether or not to return the model
    """
    try:
        from spatialpy.stochss.StochSSexport import export
    except ImportError:
        raise ImportError('StochSS export conversion not imported successfully')

    return export(spatialpy_model, path=filename, return_stochss_model=return_stochss_model)


class Model():
    """ Representation of a spatial biochemical model. """
    reserved_names = ['vol']
    special_characters = ['[', ']', '+', '-', '*', '/', '.', '^']


    def __init__(self, name=""):
        """ Create an empty SpatialPy model. """

        # The name that the model is referenced by (should be a String)
        self.name = name

        ######################
        # Dictionaries with Species, Reactions and Parameter objects.
        # Species, Reaction, and Parameter names are used as keys.
        self.listOfParameters = OrderedDict()
        self.listOfSpecies    = OrderedDict()
        self.listOfReactions  = OrderedDict()

        ######################
        # Dict that holds flattended parameters and species for
        # evaluation of expressions in the scope of the model.
        self.namespace = OrderedDict([])
        self.species_map = {}

        ######################
        self.domain = None
        self.listOfTypeIDs = [1] # starts with type '1'
        self.listOfDiffusionRestrictions = {}
        self.listOfDataFunctions = []
        self.listOfInitialConditions = []
        self.listOfBoundaryConditions = []

        ######################
        self.staticDomain = True
        self.enable_rdme = True
        self.enable_pde = True #TODO

        ######################
        self.tspan = None
        self.timestep_size = 1e-5
        self.num_timesteps = None
        self.output_freq = None


    def __str__(self):
        divider = f"\n{'*'*10}\n"

        def decorate(header):
            return f"\n{divider}{header}{divider}"

        print_string = f"{self.name}"
        if len(self.listOfSpecies):
            print_string += decorate("Species")
            for _, species in self.listOfSpecies.items():
                print_string += f"\n{str(species)}"
        if len(self.listOfInitialConditions):
            print_string += decorate("Initial Conditions")
            for initial_condition in self.listOfInitialConditions:
                print_string += f"\n{str(initial_condition)}"
        if len(self.listOfDiffusionRestrictions):
            print_string += decorate("Diffusion Restrictions")
            for species, types in self.listOfDiffusionRestrictions.items():
                print_string += f"\n{species.name} is restricted to: {str(types)}"
        if len(self.listOfParameters):
            print_string += decorate("Parameters")
            for _, parameter in self.listOfParameters.items():
                print_string += f"\n{str(parameter)}"
        if len(self.listOfReactions):
            print_string += decorate("Reactions")
            for _, reaction in self.listOfReactions.items():
                print_string += f"\n{str(reaction)}"
        if self.domain is not None:
            print_string += decorate("Domain")
            print_string += f"\n{str(self.domain)}"

        return print_string


    def run(self, number_of_trajectories=1, seed=None, timeout=None, number_of_threads=None, debug_level=0, debug=False, profile=False):
        """ Simulate the model.
        Args:
            number_of_trajectories: How many trajectories should be run.
            seed: (int) The random seed given to the solver.
            number_of_threads: (int) The number threads the solver will use.
            debug_level: (int) Level of output from the solver: 0, 1, or 2. Default: 0.
        Returns:
            A SpatialPy.Result object with the results of the simulation.
        """

        sol = Solver(self, debug_level=debug_level)

        return sol.run(number_of_trajectories=number_of_trajectories, seed=seed, timeout=timeout,
                       number_of_threads=number_of_threads, debug=debug, profile=profile)


    def set_timesteps(self, step_size, num_steps):
        """" Set the simlation time span parameters
        Args:
            step_size: float, size of each timestep in seconds
            num_steps: int, total number of steps to take

        Note: the number of output times will be num_steps+1 as the first
              output will be at time zero.
        """
        if self.timestep_size is None:
            raise ModelError("timestep_size is not set")

        steps_per_output = math.ceil(step_size/self.timestep_size)

        self.num_timesteps = math.ceil(num_steps *  steps_per_output)
        self.output_freq = steps_per_output

    def timespan(self, time_span, timestep_size=None):
        """
        Set the time span of simulation. The SSA-SDPD engine does not support
        non-uniform timespans.

        tspan : numpy ndarray
            Evenly-spaced list of times at which to sample the species
            populations during the simulation.
        """

        self.tspan = time_span
        if timestep_size is not None:
            self.timestep_size = timestep_size

        items_diff = numpy.diff(time_span)
        items = map(lambda x: round(x, 10), items_diff)
        isuniform = (len(set(items)) == 1)

        if isuniform:
            self.set_timesteps( items_diff[0], len(items_diff) )
        else:
            raise ModelError("Only uniform timespans are supported")



    def set_type(self, geometry_ivar, type_id, mass=None, nu=None, fixed=False):
        """ Add a type definition to the model.  By default, all regions are set to
        type 0.
        Args:
            geometry_ivar: an instance of a 'spatialpy.Geometry' subclass.  The 'inside()' method
                       of this object will be used to assign type_id to points.
            type_id: (usually an int) the identifier for this type
            mass: (float) the mass of each particle in the type
            nu: (float) the viscosity of each particle in the type
            fixed: (bool) are the particles in this type immobile
        Return values:
            Number of domain points that were tagged with this type_id
        """

        if self.domain is None:
            raise Exception("SpatialPy models must have a domain before types can be attached");
        if type_id not in self.listOfTypeIDs:
            # index is the "particle type", value is the "type ID"
            self.listOfTypeIDs.append(type_id)
        # apply the type to all points, set type for any points that match
        count = 0
        on_boundary = self.domain.find_boundary_points()
        for v_ndx in range(self.domain.get_num_voxels()):
            if geometry_ivar.inside( self.domain.coordinates()[v_ndx,:], on_boundary[v_ndx]):
                self.domain.type[v_ndx] = type_id
                if (mass is not None):
                    self.domain.mass[v_ndx] = mass
                if (nu is not None):
                    self.domain.nu[v_ndx] = nu
                self.domain.fixed[v_ndx] = fixed
                count +=1
        if count == 0:
            warnings.warn("Type with type_id={0} has zero particles in it".format(type_id))
        return count

    def restrict(self, species, listOfTypes):
        """ Set the diffusion coefficient to zero for 'species' in all types not in
            'listOfTypes'. This effectively restricts the movement of 'species' to
            the types specified in 'listOfTypes'.
        Args:
            species: an instance of a 'spatialpy.Species'.
            listOfTypes: a list, each object in the list should be a 'type_id'
        """
        #x = Species()
        #if not isinstance(species, Species):
        #if str(type(species)) != 'Species':
        #    raise ModelError("First argument to restrict() must be a Species object, not {0}".format(str(type(species))))
        if not isinstance(listOfTypes,list):
            self.listOfDiffusionRestrictions[species] = [listOfTypes]
        else:
            self.listOfDiffusionRestrictions[species] = listOfTypes

    def add_domain(self, domain):
        '''
        Add a domain to the model

        domain : Domain
            The Domain object to be added to the model
        '''
        if type(domain).__name__ != 'Domain':
            raise ModelError("Unexpected parameter for add_domain. Parameter must be a Domain.")

        self.domain = domain
        self.listOfTypeIDs = list(set(domain.type))

    def add_data_function(self, data_function):
        """ Add a scalar spatial function to the simulation. This is useful if you have a
            spatially varying in put to your model. Argument is a instances of subclass of the
            spatialpy.DataFunction class. It must implement a function 'map(x)' which takes a
            the spatial positon 'x' as an array, and it returns a float value.
        """
        self.listOfDataFunctions.append(data_function)

    def add_initial_condition(self, ic):
        """ Add an initial condition object to the initialization of the model."""
        self.listOfInitialConditions.append(ic)

    def add_boundary_condition(self, bc):
        """ Add an BoundaryCondition object to the model."""
        bc.model = self
        self.listOfBoundaryConditions.append(bc)

    def update_namespace(self):
        """ Create a dict with flattened parameter and species objects. """

        for param in self.listOfParameters:
            self.namespace[param]=self.listOfParameters[param].value
        # Dictionary of expressions that can be evaluated in the scope of this model.
        self.expressions = {}

    def get_species(self, sname):
        return self.listOfSpecies[sname]

    def get_num_species(self):
        return len(self.listOfSpecies)

    def get_all_species(self):
        return self.listOfSpecies

    def add_species(self, obj):
        """
        Adds a species, or list of species to the model.

        Attributes
        ----------
        obj : Species, or list of Species
            The species or list of species to be added to the model object.
        """


        if isinstance(obj, list):
            for S in obj:
                self.add_species(S)
        elif type(obj).__name__ == 'Species':
            problem = self.problem_with_name(obj.name)
            if problem is not None:
                raise problem
            self.species_map[obj] = len(self.listOfSpecies)
            self.listOfSpecies[obj.name] = obj
        else:
            raise ModelError("Unexpected parameter for add_species. Parameter must be Species or list of Species.")
        return obj


    def delete_species(self, obj):
        """ Remove a Species from model.listOfSpecies. """
        self.listOfSpecies.pop(obj)

    def delete_all_species(self):
        self.listOfSpecies.clear()

    def get_parameter(self,pname):
        try:
            return self.listOfParameters[pname]
        except:
            raise ModelError("No parameter named "+pname)

    def get_all_parameters(self):
        return self.listOfParameters

    def problem_with_name(self, name):
        if name in Model.reserved_names:
            return ModelError('Name "{}" is unavailable. It is reserved for internal SpatialPy use. Reserved Names: ({}).'.format(name, Model.reserved_names))
        if name in self.listOfSpecies:
            return ModelError('Name "{}" is unavailable. A species with that name exists.'.format(name))
        if name in self.listOfParameters:
            return ModelError('Name "{}" is unavailable. A parameter with that name exists.'.format(name))
        if name.isdigit():
            return ModelError('Name "{}" is unavailable. Names must not be numeric strings.'.format(name))
        for special_character in Model.special_characters:
            if special_character in name:
                return ModelError('Name "{}" is unavailable. Names must not contain special characters: {}.'.format(name, Model.special_characters))



    def add_parameter(self,params):
        """ Add Parameter(s) to model.listOfParameters. Input can be either a single
            Parameter object or a list of Parameter objects.
        """
        if isinstance(params,list):
            for p in params:
                self.add_parameter(p)
        else:
            #if isinstance(params, type(Parameter())):
            x = Parameter()
            if str(type(params)) == str(type(x)):
                problem = self.problem_with_name(params.name)
                if problem is not None:
                    raise problem
                # make sure that you don't overwrite an existing parameter??
                if params.name in self.listOfParameters.keys():
                    raise ParameterError("Parameter '{0}' has already been added to the model.".format(params.name))
                self.listOfParameters[params.name] = params
            else:
                #raise ParameterError("Could not resolve Parameter expression {} to a scalar value.".format(params))
                raise ParameterError("Parameter '{0}' needs to be of type '{2}', it is of type '{1}'".format(params.name,str(type(params)),str(type(x))))
        return params

    def delete_parameter(self, obj):
        self.listOfParameters.pop(obj)

    def set_parameter(self,pname,expression):
        """ Set the expression of an existing paramter. """
        p = self.listOfParameters[pname]
        p.expression = expression
        p.evaluate()

    def resolve_parameters(self):
        """ Attempt to resolve all parameter expressions to scalar floating point values.
            Must be called prior to exporting the model.  """
        self.update_namespace()
        for param in self.listOfParameters:
            try:
                self.listOfParameters[param].evaluate(self.namespace)
            except:
                raise ParameterError("Could not resolve Parameter expression "+param + "to a scalar value.")

    def delete_all_parameters(self):
        self.listOfParameters.clear()

    def add_reaction(self,reacs):
        """ Add Reaction(s) to the model. Input can be single instance, a list of instances
            or a dict with name, instance pairs. """
        if isinstance(reacs, list):
            for r in reacs:
                if r.name is None or r.name == "":
                    r.name = 'rxn' + str(uuid.uuid4()).replace('-', '_')
                self.listOfReactions[r.name] = r
        elif type(reacs).__name__ == "Reaction":
                if reacs.name is None or reacs.name == "":
                    reacs.name = 'rxn' + str(uuid.uuid4()).replace('-', '_')
                self.listOfReactions[reacs.name] = reacs
        else:
            raise ModelError("add_reaction() takes a spatialpy.Reaction object or list of objects")

    def get_reaction(self, rname):
        return self.listOfReactions[rname]

    def get_num_reactions(self):
        return len(self.listOfReactions)

    def get_all_reactions(self):
        return self.listOfReactions

    def delete_reaction(self, obj):
        self.listOfReactions.pop(obj)

    def delete_all_reactions(self):
        self.listOfReactions.clear()

    def __ne__(self, other):
        return not self.__eq__(other)

    def __eq__(self, other):
        return (self.listOfParameters == other.listOfParameters and \
            self.listOfSpecies == other.listOfSpecies and \
            self.listOfReactions == other.listOfReactions and \
            self.name == other.name)

    def create_stoichiometric_matrix(self):
        """ Generate a stoichiometric matrix in sparse CSC format. """

        if self.get_num_reactions() > 0:
            ND = numpy.zeros((self.get_num_species(), self.get_num_reactions()))
            for i, r in enumerate(self.listOfReactions):
                R = self.listOfReactions[r]
                reactants = R.reactants
                products  = R.products

                for s in reactants:
                    ND[self.species_map[s], i] -= reactants[s]
                for s in products:
                    ND[self.species_map[s], i] += products[s]

            N = scipy.sparse.csc_matrix(ND)
        else:
            N = numpy.zeros((self.get_num_species(), self.get_num_reactions()))

        return N


    def create_dependency_graph(self):
        """ Construct the sparse dependency graph. """
        # We cannot safely generate a dependency graph (without attempting to analyze the
        # propensity string itself) if the model contains custom propensities.
        mass_action_model = True
        for name, reaction in self.listOfReactions.items():
            if not reaction.massaction:
                GF = numpy.ones((self.get_num_reactions(),
                    self.get_num_reactions() + self.get_num_species()))
                mass_action_model = False

        if mass_action_model:
            GF = numpy.zeros((self.get_num_reactions(),
                self.get_num_reactions() + self.get_num_species()))
            species_map = self.species_map

            involved_species = []
            reactants = []
            for name, reaction in self.listOfReactions.items():
                temp = []
                temp2 = []
                for s in reaction.reactants:
                    temp.append(species_map[s])
                    temp2.append(species_map[s])
                for s in reaction.products:
                    temp.append(species_map[s])
                involved_species.append(temp)
                reactants.append(temp2)

            species_to_reactions = []
            for sname,species in self.listOfSpecies.items():
                temp = []
                for j, x in enumerate(reactants):
                    if species_map[species] in x:
                        temp.append(j)
                species_to_reactions.append(temp)

            reaction_to_reaction = []
            for name, reaction in self.listOfReactions.items():
                temp = []
                for s in reaction.reactants:
                    if species_to_reactions[species_map[s]] not in temp:
                        temp = temp+species_to_reactions[species_map[s]]

                for s in reaction.products:
                    if species_to_reactions[species_map[s]] not in temp:
                        temp = temp+ species_to_reactions[species_map[s]]

                temp = list(set(temp))
                reaction_to_reaction.append(temp)

            # Populate G
            for j, spec in enumerate(species_to_reactions):
                for s in spec:
                    GF[s, j] = 1

            for i,reac in enumerate(reaction_to_reaction):
                for r in reac:
                    GF[r, self.get_num_species()+i] = 1


        try:
            G = scipy.sparse.csc_matrix(GF)
        except Exception as e:
            G = GF

        return G

    def apply_initial_conditions(self):
        """ Initalize the u0 matrix (zeros) and then apply each initial condition"""
        # initalize
        ns = self.get_num_species()
        nv = self.domain.get_num_voxels()
        self.u0 = numpy.zeros((ns, nv))
        # apply initial condition functions
        for ic in self.listOfInitialConditions:
            ic.apply(self)





class Species():
    """ Model of a biochemical species. """

    def __init__(self,name=None,diffusion_constant=None,diffusion_coefficient=None,D=None):
        # A species has a name (string) and an initial value (positive integer)
        if name is None:
            raise ModelError("Species must have a name")
        else:
            self.name = name
        if diffusion_constant is not None:
            self.diffusion_constant=diffusion_constant
        elif  diffusion_coefficient is not None:
            self.diffusion_constant=diffusion_coefficient
        elif D is not None:
            self.diffusion_constant=D
        else:
            raise ModelError("Species must have a diffusion_constant")


    def __str__(self):
        print_string = f"{self.name}: {str(self.diffusion_constant)}"
        return print_string

class Parameter():
    """
        Model of a rate paramter.
        A parameter can be given as a String expression (function) or directly as a scalar value.
        If given a String expression, it should be evaluable in the namespace of a parent Model.

    """

    def __init__(self,name="",expression=None,value=None):

        self.name = name
        # We allow expression to be passed in as a non-string type. Invalid strings
        # will be caught below. It is perfectly fine to give a scalar value as the expression.
        # This can then be evaluated in an empty namespace to the scalar value.
        self.expression = expression
        if expression != None:
            self.expression = str(expression)

        self.value = value

        # self.value is allowed to be None, but not self.expression. self.value
        # might not be evaluable in the namespace of this parameter, but defined
        # in the context of a model or reaction.
        if self.expression == None:
            #raise TypeError
            self.value = 0

        if self.value is None:
            self.evaluate()

    def evaluate(self,namespace={}):
        """ Evaluate the expression and return the (scalar) value """
        try:
            self.value = (float(eval(self.expression, namespace)))
        except:
            self.value = None

    def set_expression(self,expression):
        self.expression = expression
        # We allow expression to be passed in as a non-string type. Invalid strings
        # will be caught below. It is perfectly fine to give a scalar value as the expression.
        # This can then be evaluated in an empty namespace to the scalar value.
        if expression is not None:
            self.expression = str(expression)

        if self.expression is None:
            raise TypeError

        self.evaluate()

    def __str__(self):
        print_string = f"{self.name}: {str(self.expression)}"
        return print_string


class Reaction():
    """
        Models a biochemical reaction. A reaction conatains dictinaries of species (reactants and products) and parameters.
        The reaction's propensity function needs to be evaluable and result in a non-negative scalar value
        in the namespace defined by the union of its Reactant, Product and Parameter dictionaries.

    """

    def __init__(self, name = "", reactants = {}, products = {}, propensity_function=None, massaction=None, rate=None, annotation=None,restrict_to=None):
        """
            Initializes the reaction using short-hand notation.

            Input:
                name:                       string that the model is referenced by
                parameters:                 a list of parameter instances
                propensity_function:        String with the expression for the reaction's propensity
                reactants:                  List of (species,stoichiometry) tuples
                product:                    List of (species,stoichiometry) tuples
                annotation:                 Description of the reaction (meta)

                massaction True,{False}     is the reaction of mass action type or not?
                rate                        if mass action, rate is a reference to a parameter instance.

            If massaction is set to true, propensity_function is not a valid argument. Instead, the
            propensity function is constructed automatically. For mass-action, zeroth, first and second
            order reactions are supported, attempting to used higher orders will result in an error.

            Raises: ReactionError

        """

        # Metadata
        self.name = name
        self.annotation = ""

        self.massaction = massaction

        self.propensity_function = propensity_function
        self.ode_propensity_function = propensity_function

        if self.propensity_function is None:
            if rate is None:
                errmsg = "Reaction "+self.name +": You must either set the reaction to be mass-action or specifiy a propensity function."
                raise ReactionError(errmsg)
            else:
                # If they don't give us a propensity function and do give a rate, assume mass-action.
                self.massaction = True
        else:
            if rate is not None:
                errmsg = "Reaction "+self.name +": You cannot set the propensity type to mass-action and simultaneously set a propensity function."
                raise ReactionError(errmsg)
            else:
                self.massaction = False

        self.reactants = {}
        if reactants is not None:
            for r in reactants:
                rtype = type(r).__name__
                if rtype=='instance':
                    self.reactants[r.name] = reactants[r]
                else:
                    self.reactants[r]=reactants[r]

        self.products = {}
        if products is not None:
            for p in products:
                rtype = type(p).__name__
                if rtype=='instance':
                    self.products[p.name] = products[p]
                else:
                    self.products[p]=products[p]

        if self.massaction:
            self.type = "mass-action"
            if rate is None:
                raise ReactionError("Reaction "+self.name +": A mass-action propensity has to have a rate.")
            self.marate = rate
            self.create_mass_action()
        else:
            self.type = "customized"

        self.restrict_to = restrict_to

    def __str__(self):
        print_string = f"{self.name}, Active in: {str(self.restrict_to)}"
        if len(self.reactants):
            print_string += f"\n\tReactants"
            for species, stoichiometry in self.reactants.items():
                name = species if isinstance(species, str) else species.name
                print_string += f"\n\t\t{name}: {stoichiometry}"
        if len(self.products):
            print_string += f"\n\tProducts"
            for species, stoichiometry in self.products.items():
                name = species if isinstance(species, str) else species.name
                print_string += f"\n\t\t{name}: {stoichiometry}"
        print_string += f"\n\tPropensity Function: {self.propensity_function}"
        return print_string


    def create_mass_action(self):
        """ Create a mass action propensity function given self.reactants and a single parameter value.
        """
        # We support zeroth, first and second order propensities only.
        # There is no theoretical justification for higher order propensities.
        # Users can still create such propensities if they really want to,
        # but should then use a custom propensity.
        total_stoch = 0
        for r in self.reactants:
            total_stoch += self.reactants[r]
        if total_stoch > 2:
            raise ReactionError("Reaction: A mass-action reaction cannot involve more than two of one species or one "
                                "of two species.")
        # Case EmptySet -> Y

        propensity_function = self.marate.name
        self.ode_propensity_function = self.marate.name

        # There are only three ways to get 'total_stoch==2':
        for r in self.reactants:
            # Case 1: 2X -> Y
            if self.reactants[r] == 2:
                propensity_function = ("0.5*" + propensity_function +
                                       "*" + r.name + "*(" + r.name + "-1)/vol")
            else:
                # Case 3: X1, X2 -> Y;
                propensity_function += "*" + r.name
            self.ode_propensity_function += "*" + r.name

        # Set the volume dependency based on order.
        order = len(self.reactants)
        if order == 2:
            propensity_function += "/vol"
        elif order == 0:
            propensity_function += "*vol"

        self.propensity_function = propensity_function

    def set_type(self,type):
        if type not in {'mass-action','customized'}:
            raise ReactionError("Invalid reaction type.")
        self.type = type

    def add_reactant(self,S,stoichiometry):
        if stoichiometry <= 0:
            raise ReactionError("Reaction "+self.name+"Stoichiometry must be a positive integer.")
        self.reactants[S.name]=stoichiometry

    def add_product(self,S,stoichiometry):
        self.products[S.name]=stoichiometry

    def annotate(self,annotation):
        self.annotation = annotation


# Module exceptions
class ModelError(Exception):
    pass

class SpeciesError(ModelError):
    pass

class ReactionError(ModelError):
    pass

class ParameterError(ModelError):
    pass
