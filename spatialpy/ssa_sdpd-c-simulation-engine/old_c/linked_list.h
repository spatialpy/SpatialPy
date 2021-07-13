/* *****************************************************************************************
SSA-SDPD simulation engine
Copyright 2018 Brian Drawert (UNCA)

This program is distributed under the terms of the GNU GENERAL PUBLIC LICENSE Version 3.
See the file LICENSE.txt for details.
***************************************************************************************** */
#ifndef linked_list_h
#define linked_list_h
//include <pthread.h>
#include "pthread_barrier.h"

typedef struct __linked_list_t linked_list_t;
typedef struct __node_t node_t;
typedef struct __neighbor_list_t neighbor_list_t;
typedef struct __ordered_list_t ordered_list_t;
typedef struct __neighbor_node_t neighbor_node_t;
typedef struct __ordered_node_t ordered_node_t;
#include "particle.h"

// Create data structure for a node of the list
struct __node_t {
    particle_t* data;
    node_t* next;
    node_t* prev;
};
// Data structre for the linked list type
struct __linked_list_t {
    node_t* head;
    node_t* tail;
    size_t count;
};


// Linked list structure for each particle to hold
// the list of neighbors (including distance and
// diffusion parameters: dWdr, D_i_j)
struct __neighbor_node_t {
    particle_t* data;
    double dist;
    double dWdr;
    double D_i_j;
    neighbor_node_t* next;
    neighbor_node_t* prev;
};
struct __neighbor_list_t {
    neighbor_node_t* head;
    neighbor_node_t* tail;
    size_t count;
};


//Linked list structure to hold an ordered
// list (based on 'tt', smallest values at the head)
struct __ordered_node_t {
    particle_t* data;
    double tt;
    ordered_node_t* prev;
    ordered_node_t* next;
};
struct __ordered_list_t {
    ordered_node_t* head;
    ordered_node_t* tail;
    size_t count;
};

// Functions to manipulate the linked list

//constructor
linked_list_t* create_linked_list();
neighbor_list_t* create_neighbor_list();
ordered_list_t* create_ordered_list();
// remove all elements
void empty_linked_list( linked_list_t*ll);
void empty_neighbor_list( neighbor_list_t*ll);
void empty_ordered_list( ordered_list_t*ll);
// destructor
void destroy_linked_list( linked_list_t* ll );
void destroy_neighbor_list( neighbor_list_t* ll );
void destroy_ordered_list( ordered_list_t* ll );
// add a new node to the end of the linked list
node_t* linked_list_add( linked_list_t* ll, particle_t* data_in);
neighbor_node_t* neighbor_list_add( neighbor_list_t* ll, particle_t* data_in);
ordered_node_t* ordered_list_add( ordered_list_t* ll, particle_t* data_in);
// Delete a node from the linked list
void linked_list_delete( linked_list_t* ll, node_t* to_delete);
void neighbor_list_delete( neighbor_list_t* ll, neighbor_node_t* to_delete);
void ordered_list_delete( ordered_list_t* ll, ordered_node_t* to_delete);

// search for a node by it's data field
//node* linked_list_search( linked_list* ll, char* search_string );

// get node by index
//node* linked_list_get( linked_list* ll, int index);

// get first node on list and remove it from list
//node * linked_list_pop( linked_list * ll);

// in-place bubble sort
void linked_list_sort(linked_list_t*ll, int sort_ndx);
void neighbor_list_sort(neighbor_list_t*ll);
void ordered_list_sort(ordered_list_t*ll);
// move a single element in an otherwise sorted list
void ordered_list_bubble_up_down(ordered_list_t*ll, ordered_node_t*n);



#endif /* linked_list_h*/
