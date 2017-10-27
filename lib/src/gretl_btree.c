/*
 *  gretl -- Gnu Regression, Econometrics and Time-series Library
 *  Copyright (C) 2001 Allin Cottrell and Riccardo "Jack" Lucchetti
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 */

/* The code here is a simplified version of GLib's GTree, modified
   such that ...

   
*/

#include <glib.h>
#include <glib/gslice.h>
#include "gretl_btree.h"

#define MAX_BTREE_HEIGHT 40

typedef struct _BTreeNode  BTreeNode;

struct _BTree {
    BTreeNode *root;
    GCompareFunc key_compare;
    guint nnodes;
    gint ref_count;
};

struct _BTreeNode {
    gpointer key;       /* key for this node */
    gpointer value;     /* value stored at this node */
    BTreeNode *left;    /* left subtree */
    BTreeNode *right;   /* right subtree */
    gint8 balance;      /* height (right) - height (left) */
    guint8 left_child;
    guint8 right_child;
};

static BTreeNode *b_tree_node_balance (BTreeNode *node);
static BTreeNode *b_tree_find_node (BTree *tree,
				    gconstpointer key);
static BTreeNode *b_tree_node_rotate_left (BTreeNode *node);
static BTreeNode *b_tree_node_rotate_right (BTreeNode *node);

static BTreeNode *b_tree_node_new (gpointer key,
				   gpointer value)
{
    BTreeNode *node = g_slice_new(BTreeNode);

    node->balance = 0;
    node->left = NULL;
    node->right = NULL;
    node->left_child = FALSE;
    node->right_child = FALSE;
    node->key = key;
    node->value = value;

    return node;
}

static gint tree_key_comp (gconstpointer a, gconstpointer b)
{
    return a - b > 0 ? 1 : a == b ? 0 : -1;
}

BTree *gretl_btree_new (void)
{
    BTree *tree = g_slice_new(BTree);
  
    tree->root               = NULL;
    tree->key_compare        = (GCompareFunc) tree_key_comp;
    tree->nnodes             = 0;
    tree->ref_count          = 1;
  
    return tree;
}

static inline BTreeNode *b_tree_first_node (BTree *tree)
{
    BTreeNode *tmp;

    if (!tree->root) {
	return NULL;
    }

    tmp = tree->root;

    while (tmp->left_child) {
	tmp = tmp->left;
    }

    return tmp;
}

static inline BTreeNode *b_tree_node_previous (BTreeNode *node)
{
    BTreeNode *tmp = node->left;

    if (node->left_child) {
	while (tmp->right_child) {
	    tmp = tmp->right;
	}
    }

    return tmp;
}

static inline BTreeNode *b_tree_node_next (BTreeNode *node)
{
    BTreeNode *tmp = node->right;

    if (node->right_child) {
	while (tmp->left_child) {
	    tmp = tmp->left;
	}
    }

    return tmp;
}

static void b_tree_remove_all (BTree *tree)
{
    BTreeNode *node, *next;

    g_return_if_fail(tree != NULL);

    node = b_tree_first_node(tree);

    while (node) {
	next = b_tree_node_next(node);
	g_slice_free(BTreeNode, node);
	node = next;
    }

    tree->root = NULL;
    tree->nnodes = 0;
}

void gretl_btree_destroy (BTree *tree)
{
    g_return_if_fail(tree != NULL);

    b_tree_remove_all(tree);
    g_slice_free(BTree, tree);
}

void gretl_btree_insert (BTree *tree,
			 gpointer key,
			 gpointer value)
{
    BTreeNode *node;
    BTreeNode *path[MAX_BTREE_HEIGHT];
    int idx;

    g_return_if_fail(tree != NULL);

    if (!tree->root) {
	tree->root = b_tree_node_new(key, value);
	tree->nnodes++;
	return;
    }

    idx = 0;
    path[idx++] = NULL;
    node = tree->root;

    while (1) {
	int cmp = tree->key_compare(key, node->key);
      
	if (cmp == 0) {
	    node->value = value;
	    return;
	} else if (cmp < 0) {
	    if (node->left_child) {
		path[idx++] = node;
		node = node->left;
	    } else {
		BTreeNode *child = b_tree_node_new(key, value);

		child->left = node->left;
		child->right = node;
		node->left = child;
		node->left_child = TRUE;
		node->balance -= 1;

		tree->nnodes++;

		break;
	    }
	} else {
	    if (node->right_child) {
		path[idx++] = node;
		node = node->right;
	    } else {
		BTreeNode *child = b_tree_node_new(key, value);

		child->right = node->right;
		child->left = node;
		node->right = child;
		node->right_child = TRUE;
		node->balance += 1;

		tree->nnodes++;

		break;
	    }
	}
    }

    /* Restore balance. This is the goodness of a non-recursive
     * implementation, when we are done with balancing we 'break'
     * the loop and we are done.
     */
    while (1) {
	BTreeNode *bparent = path[--idx];
	gboolean left_node = (bparent && node == bparent->left);
      
	g_assert (!bparent || bparent->left == node || bparent->right == node);

	if (node->balance < -1 || node->balance > 1) {
	    node = b_tree_node_balance (node);
	    if (bparent == NULL) {
		tree->root = node;
	    } else if (left_node) {
		bparent->left = node;
	    } else {
		bparent->right = node;
	    }
	}

	if (node->balance == 0 || bparent == NULL) {
	    break;
	}
      
	if (left_node) {
	    bparent->balance -= 1;
	} else {
	    bparent->balance += 1;
	}

	node = bparent;
    }
}

gpointer gretl_btree_lookup (BTree *tree,
			     gconstpointer key)
{
    BTreeNode *node;

    g_return_val_if_fail(tree != NULL, NULL);

    node = b_tree_find_node(tree, key);
  
    return node ? node->value : NULL;
}

static BTreeNode *b_tree_node_balance (BTreeNode *node)
{
    if (node->balance < -1) {
	if (node->left->balance > 0) {
	    node->left = b_tree_node_rotate_left(node->left);
	}
	node = b_tree_node_rotate_right (node);
    } else if (node->balance > 1) {
	if (node->right->balance < 0) {
	    node->right = b_tree_node_rotate_right(node->right);
	}
	node = b_tree_node_rotate_left(node);
    }

    return node;
}

static BTreeNode *b_tree_find_node (BTree *tree,
				    gconstpointer key)
{
    BTreeNode *node = tree->root;
    gint cmp;

    if (!node) {
	return NULL;
    }

    while (1) {
	cmp = tree->key_compare(key, node->key);
	if (cmp == 0) {
	    return node;
	} else if (cmp < 0) {
	    if (!node->left_child) {
		return NULL;
	    }
	    node = node->left;
        } else {
	    if (!node->right_child) {
		return NULL;
	    }
	    node = node->right;
	}
    }
}

static BTreeNode *b_tree_node_rotate_left (BTreeNode *node)
{
    BTreeNode *right;
    gint a_bal;
    gint b_bal;

    right = node->right;

    if (right->left_child) {
	node->right = right->left;
    } else {
	node->right_child = FALSE;
	right->left_child = TRUE;
    }
    right->left = node;

    a_bal = node->balance;
    b_bal = right->balance;

    if (b_bal <= 0) {
	if (a_bal >= 1) {
	    right->balance = b_bal - 1;
	} else {
	    right->balance = a_bal + b_bal - 2;
	}
	node->balance = a_bal - 1;
    } else {
	if (a_bal <= b_bal) {
	    right->balance = a_bal - 2;
	} else {
	    right->balance = b_bal - 1;
	}
	node->balance = a_bal - b_bal - 1;
    }

    return right;
}

static BTreeNode *b_tree_node_rotate_right (BTreeNode *node)
{
    BTreeNode *left;
    gint a_bal, b_bal;

    left = node->left;

    if (left->right_child) {
	node->left = left->right;
    } else {
	node->left_child = FALSE;
	left->right_child = TRUE;
    }
    left->right = node;

    a_bal = node->balance;
    b_bal = left->balance;

    if (b_bal <= 0) {
	if (b_bal > a_bal) {
	    left->balance = b_bal + 1;
	} else {
	    left->balance = a_bal + 2;
	}
	node->balance = a_bal - b_bal + 1;
    } else {
	if (a_bal <= -1) {
	    left->balance = b_bal + 1;
	} else {
	    left->balance = a_bal + b_bal + 2;
	}
	node->balance = a_bal + 1;
    }

    return left;
}




