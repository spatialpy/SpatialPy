/**
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
**/

/* (C) Copyright 2019 Robert Sauter
 * SPDX-License-Identifier: MIT
 */

/** Pthread-barrier implementation for macOS */

#ifndef PTHREAD_BARRIER_H
#define PTHREAD_BARRIER_H

#include <pthread.h>

#ifdef __APPLE__

#ifdef __cplusplus
extern "C" {
#endif

#ifndef PTHREAD_BARRIER_SERIAL_THREAD
# define PTHREAD_BARRIER_SERIAL_THREAD -1
#endif

typedef pthread_mutexattr_t pthread_barrierattr_t;

/* structure for internal use that should be considered opaque */
typedef struct {
	pthread_mutex_t mutex;
	pthread_cond_t cond;
	unsigned count;
	unsigned left;
	unsigned round;
} pthread_barrier_t;

int pthread_barrier_init(pthread_barrier_t *__restrict barrier,
                         const pthread_barrierattr_t * __restrict attr,
                         unsigned count);
int pthread_barrier_destroy(pthread_barrier_t *barrier);

int pthread_barrier_wait(pthread_barrier_t *barrier);

int pthread_barrierattr_init(pthread_barrierattr_t *attr);
int pthread_barrierattr_destroy(pthread_barrierattr_t *attr);
int pthread_barrierattr_getpshared(const pthread_barrierattr_t *__restrict attr,
                                   int *__restrict pshared);
int pthread_barrierattr_setpshared(pthread_barrierattr_t *attr,
                                   int pshared);


#ifdef  __cplusplus
}
#endif

#endif /* __APPLE__ */


#endif /* PTHREAD_BARRIER_H */
