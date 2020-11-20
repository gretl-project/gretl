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

/* MPI-related functions for "genr" */

static int node_replace_data (NODE *n, void *data, GretlType type)
{
    int err;

    if (n->uv != NULL) {
	err = user_var_replace_value(n->uv, data, type);
    } else {
	fprintf(stderr, "*** replace data: node uv is NULL!\n");
	if (type == GRETL_TYPE_MATRIX) {
	    err = user_matrix_replace_matrix_by_name(n->vname, data);
	} else {
	    err = E_DATA;
	}
    }

    return err;
}

static Gretl_MPI_Op real_get_reduce_op (const char *s)
{
    if (!strcmp(s, "sum")) {
	return GRETL_MPI_SUM;
    } else if (!strcmp(s, "prod")) {
	return GRETL_MPI_PROD;
    } else if (!strcmp(s, "max")) {
	return GRETL_MPI_MAX;
    } else if (!strcmp(s, "min")) {
	return GRETL_MPI_MIN;
    } else if (!strcmp(s, "hcat")) {
	return GRETL_MPI_HCAT;
    } else if (!strcmp(s, "vcat")) {
	return GRETL_MPI_VCAT;
    } else if (!strcmp(s, "acat")) {
	return GRETL_MPI_ACAT;
    } else {
	return 0;
    }
}

static Gretl_MPI_Op real_get_scatter_op (const char *s)
{
    if (!strcmp(s, "bycols")) {
	return GRETL_MPI_HSPLIT;
    } else if (!strcmp(s, "byrows")) {
	return GRETL_MPI_VSPLIT;
    } else {
	return 0;
    }
}

static Gretl_MPI_Op reduce_op_from_string (const char *s)
{
    Gretl_MPI_Op op = real_get_reduce_op(s);

    if (op == 0) {
	s = get_string_by_name(s);
	if (s != NULL) {
	    op = real_get_reduce_op(s);
	}
    }

    return op;
}

static Gretl_MPI_Op scatter_op_from_string (const char *s)
{
    Gretl_MPI_Op op = real_get_scatter_op(s);

    if (op == 0) {
	s = get_string_by_name(s);
	if (s != NULL) {
	    op = real_get_scatter_op(s);
	}
    }

    return op;
}

static gretl_matrix *get_transfer_matrix (NODE *t, int f,
					  int cmplx_ok,
					  parser *p)
{
    gretl_matrix *m = t->v.m;

    if (m == NULL) {
	p->err = E_DATA;
    } else if (m->is_complex && !cmplx_ok) {
	gretl_errmsg_sprintf("%s: %s", getsymb(f),
			     _("complex arguments/operands not supported"));
	p->err = E_CMPLX;
    }

    return m;
}

static NODE *mpi_transfer_node (NODE *l, NODE *r, NODE *r2,
				int f, parser *p)
{
    NODE *ret = NULL;
    GretlType type = 0;
    int root = 0;
    int id = 0;

    if (!gretl_mpi_initialized()) {
	gretl_errmsg_set(_("The MPI library is not loaded"));
	p->err = 1;
	return NULL;
    }

    if (f == F_MPI_SEND) {
	/* we need a destination id as second argument */
	if (l->t == MAT) {
	    type = GRETL_TYPE_MATRIX;
	} else if (l->t == NUM) {
	    type = GRETL_TYPE_DOUBLE;
	} else if (l->t == BUNDLE) {
	    type = GRETL_TYPE_BUNDLE;
	} else if (l->t == ARRAY) {
	    type = GRETL_TYPE_ARRAY;
	} else if (l->t == STR) {
	    type = GRETL_TYPE_STRING;
	} else if (l->t == LIST) {
	    type = GRETL_TYPE_LIST;
	} else {
	    p->err = E_TYPES;
	}
	if (!p->err) {
	    /* destination id */
	    id = node_get_int(r, p);
	}
    } else if (f == F_MPI_RECV) {
	/* the single argument is the source id */
	id = node_get_int(l, p);
    } else if (f == F_BCAST || f == F_REDUCE ||
	       f == F_ALLREDUCE || f == F_SCATTER) {
	/* we need the address of a variable on the left */
	if (l->t != U_ADDR) {
	    p->err = E_TYPES;
	} else {
	    /* switch to 'content' sub-node */
	    l = l->L;
	    if (umatrix_node(l)) {
		/* matrix: all operations OK */
		type = GRETL_TYPE_MATRIX;
	    } else if (ubundle_node(l) && f == F_BCAST) {
		/* bundle: only broadcast OK */
		type = GRETL_TYPE_BUNDLE;
	    } else if (uarray_node(l) && (f == F_REDUCE || f == F_BCAST)) {
		/* array: reduce and broadcast OK */
		type = GRETL_TYPE_ARRAY;
	    } else if (ustring_node(l) && f == F_BCAST) {
		/* string: only broadcast OK */
		type = GRETL_TYPE_STRING;
	    } else if (ulist_node(l) && f == F_BCAST) {
		/* list: only broadcast OK */
		type = GRETL_TYPE_LIST;
	    } else if (uscalar_node(l) && f != F_SCATTER) {
		/* scalar: all ops OK apart from scatter */
		type = GRETL_TYPE_DOUBLE;
	    } else {
		p->err = E_TYPES;
	    }
	}
	if (!p->err && f != F_ALLREDUCE) {
	    /* optional root specification */
	    NODE *rootspec = (f == F_BCAST)? r : r2;

	    if (!null_node(rootspec)) {
		root = node_get_int(rootspec, p);
	    }
	}
	if (!p->err) {
	    /* "self" id */
	    id = gretl_mpi_rank();
	}
    }

    if (p->err) {
	return NULL;
    } else if (f == F_MPI_SEND) {
	void *sendp = NULL;

	if (type == GRETL_TYPE_MATRIX) {
	    sendp = l->v.m;
	} else if (type == GRETL_TYPE_BUNDLE) {
	    sendp = l->v.b;
	} else if (type == GRETL_TYPE_ARRAY) {
	    sendp = l->v.a;
	} else if (type == GRETL_TYPE_STRING) {
	    sendp = l->v.str;
	} else if (type == GRETL_TYPE_LIST) {
	    sendp = l->v.ivec;
	} else if (type == GRETL_TYPE_DOUBLE) {
	    sendp = &l->v.xval;
	} else {
	    p->err = E_TYPES;
	}
	if (!p->err) {
	    ret = aux_scalar_node(p);
	}
	if (!p->err) {
	    p->err = ret->v.xval = gretl_mpi_send(sendp, type, id);
	}
    } else if (f == F_MPI_RECV) {
	void *recvp;

	recvp = gretl_mpi_receive(id, &type, &p->err);

	if (!p->err) {
	    if (type == GRETL_TYPE_MATRIX) {
		ret = aux_matrix_node(p);
		if (!p->err) {
		    ret->v.m = recvp;
		}
	    } else if (type == GRETL_TYPE_BUNDLE) {
		ret = aux_bundle_node(p);
		if (!p->err) {
		    ret->v.b = recvp;
		}
	    } else if (type == GRETL_TYPE_ARRAY) {
		ret = aux_array_node(p);
		if (!p->err) {
		    ret->v.a = recvp;
		}
	    } else if (type == GRETL_TYPE_DOUBLE) {
		ret = aux_scalar_node(p);
		if (!p->err) {
		    ret->v.xval = *(double *) recvp;
		}
	    } else if (type == GRETL_TYPE_STRING) {
		ret = aux_string_node(p);
		if (!p->err) {
		    ret->v.str = recvp;
		}
	    } else if (type == GRETL_TYPE_LIST) {
		ret = aux_list_node(p);
		if (!p->err) {
		    ret->v.ivec = recvp;
		}
	    } else if (type == GRETL_TYPE_INT) {
		ret = aux_scalar_node(p);
		if (!p->err) {
		    ret->v.xval = *(int *) recvp;
		}
	    }
	}
    } else if (f == F_BCAST) {
	gretl_matrix *m = NULL;
	gretl_bundle *b = NULL;
	gretl_array *a = NULL;
	char *s = NULL;
	int *list = NULL;
	double x = NADBL;
	void *bcastp = NULL;

	if (type == GRETL_TYPE_MATRIX) {
	    if (id == root) {
		m = l->v.m;
	    }
	    bcastp = &m;
	} else if (type == GRETL_TYPE_BUNDLE) {
	    if (id == root) {
		b = l->v.b;
	    }
	    bcastp = &b;
	} else if (type == GRETL_TYPE_ARRAY) {
	    if (id == root) {
		a = l->v.a;
	    }
	    bcastp = &a;
	} else if (type == GRETL_TYPE_STRING) {
	    if (id == root) {
		s = l->v.str;
	    }
	    bcastp = &s;
	} else if (type == GRETL_TYPE_LIST) {
	    if (id == root) {
		list = l->v.ivec;
	    }
	    bcastp = &list;
	} else if (type == GRETL_TYPE_DOUBLE) {
	    x = l->v.xval;
	    bcastp = &x;
	} else {
	    p->err = E_TYPES;
	}

	if (!p->err) {
	    ret = aux_scalar_node(p);
	}
	if (!p->err) {
	    p->err = gretl_mpi_bcast(bcastp, type, root);
	    if (!p->err && id != root) {
		if (type == GRETL_TYPE_DOUBLE) {
		    p->err = node_replace_scalar(l, x);
		} else if (type == GRETL_TYPE_MATRIX) {
		    p->err = node_replace_data(l, m, type);
		} else if (type == GRETL_TYPE_BUNDLE) {
		    p->err = node_replace_data(l, b, type);
		} else if (type == GRETL_TYPE_ARRAY) {
		    p->err = node_replace_data(l, a, type);
		} else if (type == GRETL_TYPE_STRING) {
		    p->err = node_replace_data(l, s, type);
		} else if (type == GRETL_TYPE_LIST) {
		    p->err = node_replace_data(l, list, type);
		}
	    }
	    ret->v.xval = p->err;
	}
    } else if (f == F_REDUCE || f == F_ALLREDUCE) {
	ret = aux_scalar_node(p);
	if (!p->err) {
	    Gretl_MPI_Op op = reduce_op_from_string(r->v.str);
	    gretlopt opt = (f == F_REDUCE)? OPT_NONE : OPT_A;
	    gretl_matrix *lm = NULL;
	    gretl_matrix *m = NULL;
	    gretl_array *a = NULL;
	    double x = NADBL;

	    if (type == GRETL_TYPE_ARRAY) {
		p->err = gretl_array_mpi_reduce(l->v.a, &a, op, root);
	    } else if (type == GRETL_TYPE_MATRIX) {
		lm = get_transfer_matrix(l, f, 0, p);
		if (!p->err) {
		    p->err = gretl_matrix_mpi_reduce(lm, &m, op, root, opt);
		}
	    } else {
		p->err = gretl_scalar_mpi_reduce(l->v.xval, &x, op, root, opt);
	    }
	    if (!p->err && (id == root || f == F_ALLREDUCE)) {
		if (type == GRETL_TYPE_ARRAY) {
		    p->err = node_replace_data(l, a, type);
		} else if (type == GRETL_TYPE_MATRIX) {
		    p->err = node_replace_data(l, m, type);
		} else {
		    p->err = node_replace_scalar(l, x);
		}
	    }
	    ret->v.xval = p->err;
	}
    } else if (f == F_SCATTER) {
	gretl_matrix *lm = get_transfer_matrix(l, f, 1, p);

	if (!p->err) {
	    ret = aux_scalar_node(p);
	}
	if (!p->err) {
	    Gretl_MPI_Op op = scatter_op_from_string(r->v.str);
	    gretl_matrix *m = NULL;

	    p->err = ret->v.xval = gretl_matrix_mpi_scatter(lm, &m, op, root);
	    if (!p->err) {
		p->err = node_replace_data(l, m, GRETL_TYPE_MATRIX);
	    }
	}
    } else {
	gretl_errmsg_set("MPI function not yet supported");
	p->err = 1;
    }

    return ret;
}

static NODE *mpi_barrier_node (parser *p)
{
    NODE *ret = NULL;

    if (!gretl_mpi_initialized()) {
	gretl_errmsg_set(_("The MPI library is not loaded"));
	p->err = 1;
	return NULL;
    } else {
	ret = aux_scalar_node(p);
	if (ret != NULL) {
	    ret->v.xval = gretl_mpi_barrier();
	}
    }

    return ret;
}
