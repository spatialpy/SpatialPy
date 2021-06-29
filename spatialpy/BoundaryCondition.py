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
from spatialpy.Model import ModelError


class BoundaryCondition():
    """ Set spatial regions of the domain where a property of
        particles are held constant (updated each simulation step)
    """
    def __init__(self,
                 xmin=None, xmax=None,
                 ymin=None, ymax=None,
                 zmin=None, zmax=None,
                 type_id=None,
                 species=None,
                 deterministic=True,
                 property=None,
                 value=None,
                 model=None):
        """
            Conditions (one or more of the following must be set):
                 xmin, xmax: (float) min or max value in the x dimension
                 ymin, ymax: (float) min or max value in the y dimension
                 zmin, zmax: (float) min or max value in the z dimension
                 type_id: type (subdomain) of the partciles
            Targets (one of the following must be set):
                property: (str), 'nu', 'rho','v'
                spesicies: (str) name of a chemical species.
                           Must also set deterministic=True/False flag.
            Assignment:
                value: (float or float[3]), value property will take in region
                       defined by the conditions

        """
        self.xmin = xmin
        self.xmax = xmax
        self.ymin = ymin
        self.ymax = ymax
        self.zmin = zmin
        self.zmax = zmax
        self.type_id = type_id
        self.species = species
        self.property = property
        self.deterministic = deterministic
        self.value = value
        self.model = model


    def expression(self):
        if( self.species is not None and self.property is not None):
            raise ModelError("Can not set both species and property")
        if self.value is None:
            raise ModelError("Must set value")
        cond=[]
        if(self.xmin is not None): cond.append("(me->x[0] >= {0})".format(self.xmin))
        if(self.xmax is not None): cond.append("(me->x[0] <= {0})".format(self.xmax))
        if(self.ymin is not None): cond.append("(me->x[1] >= {0})".format(self.ymin))
        if(self.ymax is not None): cond.append("(me->x[1] <= {0})".format(self.ymax))
        if(self.zmin is not None): cond.append("(me->x[2] >= {0})".format(self.zmin))
        if(self.zmax is not None): cond.append("(me->x[2] <= {0})".format(self.zmax))
        if(self.type_id is not None): cond.append("(me->type == {0})".format(int(self.type_id)))
        if(len(cond)==0): raise ModelError('need at least one condition on the BoundaryCondition')
        bcstr = "if(" + '&&'.join(cond) + "){"
        if self.species is not None:
            if self.deterministic:
                s_ndx = self.model.species_map[self.model.listOfSpecies[self.species]]
                bcstr += "me->C[{0}] = {1};".format(s_ndx,self.value)
            else:
                raise Exception("BoundaryConditions don't work for stochastic species yet")
        elif self.property is not None:
            if(self.property == 'v'):
                for i in range(3):
                    bcstr+= "me->v[{0}]={1};".format(i,self.value[i])
            elif(self.property == 'nu'):
                bcstr+= "me->nu={0};".format(self.value)
            elif(self.property == 'rho'):
                bcstr+= "me->rho={0};".format(self.value)
            else:
                raise Exception("Unable handle boundary condition for property '{0}'".format(self.property))
        bcstr+= "}"
        return bcstr
