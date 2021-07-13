/* *****************************************************************************
SSA-SDPD simulation engine
Copyright 2018 Brian Drawert (UNCA)

This program is distributed under the terms of the GNU General Public License.
See the file LICENSE.txt for details.
***************************************************************************** */
#include <string.h>  // for strcmp and strcpy
#include <stdlib.h>  // for malloc and free
#include <stdio.h>   // for printf
#include "linked_list.h"
#include "particle.h"
#include <math.h>
#include "dSFMT/dSFMT.h"


//#define DEBUG_UPDATE


//constructor
linked_list_t* create_linked_list(){
    linked_list_t* ll = (linked_list_t*) malloc( sizeof(linked_list_t));
    ll->count = 0;
    ll->head = NULL;
    ll->tail = NULL;
    return ll;
}
neighbor_list_t* create_neighbor_list(){
    neighbor_list_t* ll = (neighbor_list_t*) malloc( sizeof(neighbor_list_t));
    ll->count = 0;
    ll->head = NULL;
    ll->tail = NULL;
    return ll;
}
ordered_list_t* create_ordered_list(){
    ordered_list_t* ll = (ordered_list_t*) malloc( sizeof(ordered_list_t));
    ll->count = 0;
    ll->head = NULL;
    ll->tail = NULL;
    return ll;
}

// destructor
void destroy_linked_list( linked_list_t* ll ){
    // empty the linked list
    empty_linked_list(ll);
    // un-allocate the memory
    free(ll);
}
void destroy_neighbor_list( neighbor_list_t* ll ){
    // empty the linked list
    empty_neighbor_list(ll);
    // un-allocate the memory
    free(ll);
}
void destroy_ordered_list( ordered_list_t* ll ){
    // empty the linked list
    empty_ordered_list(ll);
    // un-allocate the memory
    free(ll);
}
// empty
void empty_linked_list( linked_list_t*ll){
    while( ll->count > 0){
        linked_list_delete( ll, ll->head );
    }
}
void empty_neighbor_list( neighbor_list_t*ll){
    while( ll->count > 0){
        neighbor_list_delete( ll, ll->head );
    }
}
void empty_ordered_list( ordered_list_t*ll){
    while( ll->count > 0){
        ordered_list_delete( ll, ll->head );
    }
}

// add a new node to the end of the linked list
node_t* linked_list_add( linked_list_t* ll, particle_t* data_in){
    node_t* n = (node_t *) malloc( sizeof(node_t) );
    n->data = data_in;
    n->next = NULL;
    n->prev = NULL;
    // Traverse the list to find the end node
    if(ll->head == NULL){
        ll->head = n;
        ll->tail = n;
    }else{
        ll->tail->next = n;
        n->prev = ll->tail;
        ll->tail = n;
    }
    // increase the size of the list
    ll->count++;
    return n;
}
neighbor_node_t* neighbor_list_add( neighbor_list_t* ll, particle_t* data_in){
    neighbor_node_t* n = (neighbor_node_t *) malloc( sizeof(neighbor_node_t) );
    n->data = data_in;
    n->next = NULL;
    n->prev = NULL;
    // Traverse the list to find the end node
    if(ll->head == NULL){
        ll->head = n;
        ll->tail = n;
    }else{
        ll->tail->next = n;
        n->prev = ll->tail;
        ll->tail = n;
    }
    // increase the size of the list
    ll->count++;
    return n;
}
ordered_node_t* ordered_list_add( ordered_list_t* ll, particle_t* data_in){
    ordered_node_t* n = (ordered_node_t *) malloc( sizeof(ordered_node_t) );
    n->data = data_in;
    n->next = NULL;
    n->prev = NULL;
    // Traverse the list to find the end node
    if(ll->head == NULL){
        ll->head = n;
        ll->tail = n;
    }else{
        ll->tail->next = n;
        n->prev = ll->tail;
        ll->tail = n;
    }
    // increase the size of the list
    ll->count++;
    return n;
}

// Delete a node from the linked list
void linked_list_delete( linked_list_t* ll, node_t* to_delete){
    node_t* prev_node;
    if( ll->head == NULL){
        printf("Error, linked_list_delete() empty list\n");
        return;
    }else if( to_delete == ll->head ){
        ll->head = ll->head->next;
    }else{
        for( prev_node=ll->head; prev_node->next!=NULL; prev_node=prev_node->next ){
            if(prev_node->next == to_delete){
                break;
            }
        }
        if( prev_node->next == NULL){
            printf("Error, linked_list_delete(), could not find item in list\n");
            return;
        }
        prev_node->next = to_delete->next;  // connect the list
        if(prev_node->next != NULL){ prev_node->next->prev = prev_node; }
    }

    //free and reduce size
    ll->count--;
    free(to_delete);
}
void neighbor_list_delete( neighbor_list_t* ll, neighbor_node_t* to_delete){
    neighbor_node_t* prev_node;
    if( ll->head == NULL){
        printf("Error, neighbor_list_t_delete() empty list\n");
        return;
    }else if( to_delete == ll->head ){
        ll->head = ll->head->next;
    }else{
        for( prev_node=ll->head; prev_node->next!=NULL; prev_node=prev_node->next ){
            if(prev_node->next == to_delete){
                break;
            }
        }
        if( prev_node->next == NULL){
            printf("Error, neighbor_list_t_delete(), could not find item in list\n");
            return;
        }
        prev_node->next = to_delete->next;  // connect the list
        if(prev_node->next != NULL){ prev_node->next->prev = prev_node; }
    }

    //free and reduce size
    ll->count--;
    free(to_delete);
}
void ordered_list_delete( ordered_list_t* ll, ordered_node_t* to_delete){
    ordered_node_t* prev_node;
    if( ll->head == NULL){
        printf("Error, ordered_list_t_delete() empty list\n");
        return;
    }else if( to_delete == ll->head ){
        ll->head = ll->head->next;
    }else{
        for( prev_node=ll->head; prev_node->next!=NULL; prev_node=prev_node->next ){
            if(prev_node->next == to_delete){
                break;
            }
        }
        if( prev_node->next == NULL){
            printf("Error, ordered_list_t_delete(), could not find item in list\n");
            return;
        }
        prev_node->next = to_delete->next;  // connect the list
        if(prev_node->next != NULL){ prev_node->next->prev = prev_node; }
    }

    //free and reduce size
    ll->count--;
    free(to_delete);
}

node_t *split(node_t *head){
    node_t *fast = head, *slow = head ;
    while(fast->next && fast->next->next){
	fast = fast->next->next ;
	slow = slow->next ;
    }
    node_t *temp = slow->next ;
    slow->next = NULL ;
    return temp ;
}

neighbor_node_t *neighbor_split(neighbor_node_t *head){
    neighbor_node_t *fast = head, *slow = head ;
    while(fast->next && fast->next->next){
	fast = fast->next->next ;
	slow = slow->next ;
    }
    neighbor_node_t *temp = slow->next ;
    slow->next = NULL ;
    return temp ;
}

ordered_node_t *ordered_split(ordered_node_t *head){
    ordered_node_t *fast = head, *slow = head ;
    while(fast->next && fast->next->next){
	fast = fast->next->next ;
	slow = slow->next ;
    }
    ordered_node_t *temp = slow->next ;
    slow->next = NULL ;
    return temp ;
}

node_t *merge(node_t *first, node_t *second){
    if (!first){return second ;}
    if (!second){return first ;}
    if (first->data->x[0] < second->data->x[0]){
	first->next = merge(first->next, second) ;
	first->next->prev = first ;
	return first ;
    }else{
	second->next = merge(first, second->next) ;
	second->next->prev = second ;
	second->prev = NULL ;
	return second ;
    }
}

neighbor_node_t *neighbor_merge(neighbor_node_t *first, neighbor_node_t *second){
    if (!first){return second ;}
    if (!second){return first ;}
    if (first->dist < second->dist){
	first->next = neighbor_merge(first->next, second) ;
	first->next->prev = first ;
	return first ;
    }else{
	second->next = neighbor_merge(first, second->next) ;
	second->next->prev = second ;
	second->prev = NULL ;
	return second ;
    }
}

ordered_node_t *ordered_merge(ordered_node_t *first, ordered_node_t *second){
    if (!first){return second ;}
    if (!second){return first ;}
    if (first->tt < second->tt){
	first->next = ordered_merge(first->next, second) ;
	first->next->prev = first ;
	return first ;
    }else{
	second->next = ordered_merge(first, second->next) ;
	second->next->prev = second ;
	second->prev = NULL ;
	return second ;
    }
}

node_t* linked_list_sort__mergesort(node_t *head){
    if(!head || !head->next){return head ;}
    node_t *second = split(head) ;
    
    head = linked_list_sort__mergesort(head) ;
    second = linked_list_sort__mergesort(second) ;

    return merge(head, second) ;
}

neighbor_node_t *neighbor_list_sort__mergesort(neighbor_node_t *head){
    if(!head || !head->next){return head ;}
    neighbor_node_t *second = neighbor_split(head) ;
    
    head = neighbor_list_sort__mergesort(head) ;
    second = neighbor_list_sort__mergesort(second) ;

    return neighbor_merge(head, second) ;
}

ordered_node_t *ordered_list_sort__mergesort(ordered_node_t *head){
    if(!head || !head->next){return head ;}
    ordered_node_t *second = ordered_split(head) ;
    
    head = ordered_list_sort__mergesort(head) ;
    second = ordered_list_sort__mergesort(second) ;

    return ordered_merge(head, second) ;
}

void linked_list_sort(linked_list_t*ll){
    ll->head = linked_list_sort__mergesort(ll->head) ;
}

void neighbor_list_sort(neighbor_list_t*ll){
    ll->head = neighbor_list_sort__mergesort(ll->head) ;
}

void ordered_list_sort(ordered_list_t*ll){
    ll->head = ordered_list_sort__mergesort(ll->head) ;
}


// move a single element in an otherwise sorted list
void ordered_list_bubble_up_down(ordered_list_t*ll, ordered_node_t*n){
    ordered_list_sort(ll) ;
/*
    ordered_node_t*n1 = n->next;

    // Remove node from current position
    //printf("BUBBLE...") ;
   
    ordered_node_t*before = n->prev;
    ordered_node_t*after  = n->next;
    if(before != NULL){		// If node before, connect that to one after
        before->next = after;
	printf("1") ;
    }else{
        ll->head = ll->head->next;	//else node is head
	ll->head->prev = NULL ;
	printf("2") ;
    }
    if(after != NULL){		// If node after, connect that to one before
        after->prev = before;		// else node is tail
	printf("3") ;
    }else{			// if nothing after (is tail), set tail to node before
        ll->tail = ll->tail->prev;
	ll->tail->next = NULL ;
	printf("4") ;
    }
    // Find new position
    //      if tt==inf, move to end
    if(isinf(n->tt) || n->tt >= ll->tail->tt){
        // move to end
	if(n == ll->head){
		ll->head = n->next ;
	}
        ll->tail->next = n;
        n->prev = ll->tail;
        n->next = NULL;
        ll->tail = n;
    }
    //  check to move to beginning
    else if(n->tt <= ll->head->tt || isinf(ll->head->tt)){
        // move to beginning
        ll->head->prev = n;
        n->next = ll->head;
        n->prev = NULL;
        ll->head = n;
    }
    //      Check if we move down, move down linearly
    else if(n->next != NULL && n->next->tt < n->tt){
	while(n1!=NULL && n1->tt < n->tt){
	    n1=n1->next ;
	}
                n->next = n1;
                n->prev = n1->prev;
                n1->prev->next = n;
                n1->prev = n;
    }
    else {
	n1 = n->prev ;
	while(n1!=NULL && n1->tt > n->tt){
		n1=n1->prev ;
	}
                n->prev = n1;
                n->next = n1->next;
                if(n1->next != NULL){n1->next->prev = n;}
                n1->next = n;
    }

*/
}

