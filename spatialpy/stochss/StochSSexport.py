import ast
import json
import copy

def __add_boundary_conditions(model, boundary_conditions):
    for boundary_condition in boundary_conditions:
        s_bound_cond = {"compID":model['defaultID'],
                        "name": type(boundary_condition).__name__,
                        "expression": boundary_condition.expression(),
                        "annotation": ""}
        model['boundaryConditions'].append(s_bound_cond)
        model['defaultID'] += 1


def __add_domain(model, domain):
    boundary_condition = {"reflext_x":True, "reflect_y":True, "reflect_z":True}
    particles = __get_particles(domain=domain)
    domain = {
        "size": None,
        "rho_0": domain.rho0,
        "c_0": domain.c0,
        "p_0": domain.P0,
        "gravity": [0] * 3 if domain.gravity is None else domain.gravity,
        "x_lim": list(domain.xlim),
        "y_lim": list(domain.ylim),
        "z_lim": list(domain.zlim),
        "boundary_condition": boundary_condition,
        "types": [],
        "particles": particles
    }

    model['domain'] = domain


def __add_initial_conditions(model, initial_conditions, types):
    for initial_condition in initial_conditions:
        species = __get_species(species=model['species'], name=initial_condition.species.name)
        s_initial_condition = {"specie":species,
                               "count":initial_condition.count}
        if "Place" in str(type(initial_condition)):
            initial_condition['icType'] = "Place"
            initial_condition['x'] = initial_condition.location[0]
            initial_condition['y'] = initial_condition.location[1]
            initial_condition['z'] = initial_condition.location[2]
            initial_condition['types'] = types
        else:
            if initial_condition.types is None:
                initial_condition['types'] = types
            else:
                initial_condition['types'] = initial_condition.types
            if "Scatter" in str(type(initial_condition)):
                initial_condition['icType'] = "Scatter"
            else:
                initial_condition['icType'] = "Distribute Uniformly per Voxel"
            initial_condition['x'] = 0
            initial_condition['y'] = 0
            initial_condition['z'] = 0

        model['initialConditions'].append(s_initial_condition)


def __add_parameters(model, parameters):
    for name, parameter in parameters.items():
        try:
            expression = ast.literal_eval(parameter.expression)
        except:
            expression = parameter.expression
        s_parameter = {"compID":model['defaultID'],
                       "name":name,
                       "expression":str(expression),
                       "annotation": ""}
        model['parameters'].append(s_parameter)
        model['defaultID'] += 1


def __add_reactions(model, reactions, types):
    for name, reaction in reactions.items():
        if reaction.restrict_to is not None:
            types = reaction.restrict_to
        s_reaction = {"compID":model['defaultID'],
                      "name":name,
                      "reactionType": "custom-propensity",
                      "massaction": False,
                      "propensity": reaction.propensity_function,
                      "annotation": "",
                      "rate": {},
                      "types":types if isinstance(types, list) else [types],
                      "reactants": [],
                      "products": []}

        for key in ['reactants', 'products']:
            __add_stoich_species(s_reaction=s_reaction, reaction=reaction,
                                                key=key, species=model['species'])
        __add_summary(reaction=s_reaction)

        model['reactions'].append(s_reaction)
        model['defaultID'] += 1


def __add_species(model, species, diffusion_restrictions, types):
    for name, specie in species.items():
        if specie in diffusion_restrictions.keys():
            types = diffusion_restrictions[specie]
        s_species = {"compID":model['defaultID'],
                     "name":name,
                     "value":0,
                     "mode":None,
                     "switchTol": 0.03,
                     "switchMin": 100,
                     "isSwitchTol": True,
                     "annotation": "",
                     "diffusionConst":specie.diffusion_constant,
                     "types": types}
        model['species'].append(s_species)
        model['defaultID'] += 1


def __add_stoich_species(s_reaction, reaction, key, species):
    source = reaction.reactants if key == "reactants" else reaction.products
    for specie, ratio in source.items():
        stoich_species = {"ratio":ratio,
                          "specie":__get_species(species=species, name=specie.name)}
        s_reaction[key].append(stoich_species)


def __add_summary(reaction):
    r_summary = __build_element(reaction['reactants'])
    p_summary = __build_element(reaction['products'])
    reaction['summary'] = f"{r_summary} \\rightarrow {p_summary}"


def __add_types(model, types):
    default_type = {"fixed":False, "mass":1, "name":"Un-Assigned",
                    "nu":0, "typeID":0, "volume":1}
    model['domain']['types'].append(default_type)
    for sp_type in types:
        s_type = copy.deepcopy(default_type)
        s_type['typeID'] = sp_type
        s_type['name'] = str(sp_type)

        model['domain']['types'].append(s_type)


def __build_element(stoich_species):
    if not stoich_species:
        return "\\emptyset"

    elements = []
    for species in stoich_species:
        name = species['specie']['name']
        ratio = species['ratio']
        element = f"{ratio}{name}" if ratio > 1 else name
        elements.append(element)
    return '+'.join(elements)


def __get_particles(domain):
    s_particles = []
    for i, point in enumerate(domain.vertices):
        s_particle = {"fixed":bool(domain.fixed[i]),
                      "mass":domain.mass[i],
                      "nu":domain.nu[i],
                      "particle_id":i,
                      "point":list(point),
                      "type":int(domain.type[i]),
                      "volume":domain.vol[i]}

        s_particles.append(s_particle)
    return s_particles


def __get_species(species, name):
    return list(filter(lambda specie: specie['name'] == name, species))[0]


def __write_to_file(model, path):
    with open(path, "w") as model_file:
        json.dump(model, model_file, indent=4, sort_keys=True)


def export(model, path=None, return_stochss_model=False):
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

    if path is None:
        path = f"{model.name}.smdl"

    s_model = {"is_spatial": True,
               "defaultID": 1,
               "defaultMode": "",
               "annotation": "",
               "volume": 1,
               "modelSettings": {
                    "endSim": model.tspan[-1],
                    "timeStep": model.timestep_size
                },
                "species": [],
                "initialConditions": [],
                "parameters": [],
                "reactions": [],
                "rules": [],
                "eventsCollection": [],
                "functionDefinitions": [],
                "boundaryConditions": []
              }

    __add_domain(model=s_model, domain=model.domain)
    s_model['domain']['static'] = model.staticDomain
    __add_types(model=s_model, types=model.listOfTypeIDs)
    __add_boundary_conditions(model=s_model, boundary_conditions=model.listOfBoundaryConditions)
    __add_species(model=s_model, species=model.get_all_species(), types=model.listOfTypeIDs,
                  diffusion_restrictions=model.listOfDiffusionRestrictions)
    __add_initial_conditions(model=s_model, types=model.listOfTypeIDs,
                             initial_conditions=model.listOfInitialConditions)
    __add_parameters(model=s_model, parameters=model.get_all_parameters())
    __add_reactions(model=s_model, reactions=model.get_all_reactions(), types=model.listOfTypeIDs)

    if return_stochss_model:
        return s_model

    __write_to_file(model=s_model, path=path)
    return path
