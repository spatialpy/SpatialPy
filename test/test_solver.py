#!/usr/bin/env python3

import pickle
import unittest

import spatialpy


class diffusion_debug(spatialpy.Model):

    def __init__(self, model_name="diffusion_debug_test"):
        spatialpy.Model.__init__(self, model_name)

        D_const = 0.01

        A = spatialpy.Species(name="A", diffusion_constant=D_const)
        self.add_species([A])

        self.domain = spatialpy.Domain.create_2D_domain(
            xlim=[-1, 1], ylim=[-1, 1], nx=50, ny=50, type_id=1.0,
            mass=1.0, nu=1.0, fixed=True,  rho0=1.0, c0=1.0, P0=1.0
        )

        self.add_initial_condition(spatialpy.PlaceInitialCondition(A, 100000, [0,0,0]))

        self.timestep_size=.1
        self.num_timesteps=10
        self.output_freq=1


# class testPeriodicDiffusion(spatialpy.Model):
#     def __init__(self, model_name="test1D"):
#         spatialpy.Model.__init__(self, model_name)
#         X = self.add_species(spatialpy.Species(name="X",  diffusion_constant=0.001))
#         self.domain = spatialpy.Domain.generate_unit_interval_mesh(nx=100, periodic=True)
#         self.add_initial_condition(spatialpy.PlaceInitialCondition(X, 1000))
#         #self.set_initial_condition_place_near({X:1000}, 0.1)
#         self.timespan(range(10))


class TestSolverFunctionality(unittest.TestCase):

    def setUp(self):
        self.model = diffusion_debug()
        #self.periodic_model = testPeriodicDiffusion()

    def test_solver_io(self):
        """ Test that the initial value in the solver output file is the same as the input initial value. """
        model = diffusion_debug()
        result = model.run()
        A = result.get_species("A", 0)
        self.assertFalse((A-model.u0).any())

    # def test_zero_diffusion(self):
    #     """ Test that nothing happens if the diffusion is set to zero. """
    #     model = diffusion_debug(diffusion_constant=0.0)
    #     result = model.run()
    #     A = result.get_species("A", -1)
    #     self.assertFalse((A - model.u0).any())

    def test_same_seed(self):
        """ Test that the output is the same if the same seed is used, edxplicit solver creation  """
        solver = spatialpy.Solver(self.model)
        result1 = solver.run(seed=1)
        result2 = solver.run(seed=1)
        self.assertTrue(result1 == result2)

    def test_same_seed2(self):
        """ Test that the output is the same if the same seed is used, use model.run()  """
        result1 = self.model.run(seed=1)
        result2 = self.model.run(seed=1)
        self.assertTrue(result1 == result2)

    def test_different_seeds(self):
        """ Test that the output is different if different seeds are used. """
        solver = spatialpy.Solver(self.model)
        result1 = solver.run(seed=1)
        result2 = solver.run(seed=100)
        self.assertFalse(result1 == result2)

    def test_default_seed(self):
        """ Test that the output is different if no seed is given (default set on C level). """
        solver = spatialpy.Solver(self.model)
        result1 = solver.run()
        result2 = solver.run()
        self.assertFalse(result1 == result2)

    def test_mesh_pickle(self):
        meshstr = pickle.dumps(self.model.domain)
        domain = pickle.loads(meshstr)

    def test_model_pickle(self):
        """ Test that the model is picklable. We do not compare models directly, but rather the results after simulation. """
        model = self.model
        model_str = pickle.dumps(model)
        model2 = pickle.loads(model_str)
        result1 = model.run(seed=1)
        result2 = model2.run(seed=1)
        self.assertTrue(result1 == result2)

    def test_solver_pickle(self):
        """ Test that the model, solver and result objects are pickleable. """
        sol = spatialpy.Solver(self.model)
        sol_str = pickle.dumps(sol)
        sol2 = pickle.loads(sol_str)
        result1 = sol.run(seed=1)
        result2 = sol2.run(seed=1)
        self.assertTrue(result1 == result2)

    def test_result_pickle(self):
        """ Test that the result object is picklable. """
        sol = spatialpy.Solver(self.model)
        result = sol.run(seed=1)
        result_str = pickle.dumps(result)
        result2 = pickle.loads(result_str)
        self.assertTrue(result == result2)

    def test_run_ensemble(self):
        """ Test the running of ensembles of runs """
        result_list = self.model.run(3)
        self.assertEqual(len(result_list), 3)

    # def test_1D_periodic_boundary(self):
    #     """ Test if periodic boundary conditions are working. """
    #     result = self.periodic_model.run()
    #     self.assertTrue(result.get_species("X", timepoints=1)[-1] > 0)

    # def test_1D_periodic_boundary_pickle(self):
    #     """ Test if periodic boundary conditions are working. """
    #     model_str = pickle.dumps(self.periodic_model)
    #     model2 = pickle.loads(model_str)
    #     self.periodic_model.assemble()
    #     model2.assemble()
    #     self.assertEqual(self.periodic_model.domain.get_num_dof_voxels(), model2.domain.get_num_dof_voxels())


if __name__ == '__main__':
    unittest.main()
