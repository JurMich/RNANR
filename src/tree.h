#ifndef TREE_H_
#define TREE_H_

#define TYPE double

/* explanation: each node (tr_node) holds a pointer to a start of a linked
 * list to the next node, link to previous node, Boltzmann's factor
 * for given flat structure and itentifier, which corresponds to an id of
 * chosen flat structure
 */

typedef struct tr_node tr_node;
typedef struct ll_node ll_node;

struct tr_node {int id; tr_node * parent; ll_node * ll_start; TYPE weight;};
struct ll_node {int id; ll_node * ll_next; tr_node * next_node;};

tr_node * create_root();
TYPE tr_node_weight(int id, tr_node * par_node);
TYPE total_weight_par(tr_node * par_node);
tr_node * add_if_nexists(int id, tr_node*par_node);
int child_exists(int id, tr_node*par_node);

tr_node * traceback_to_root(tr_node * leaf, TYPE struct_weight);

#endif