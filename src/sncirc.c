//#include "cost_general_functions.c"
#include <R.h>
#include <Rmath.h>
#include <Rinternals.h> // RK addition
#include <R_ext/RS.h>  // RK addition
#include <R_ext/Lapack.h> // RK addition
#include <R_ext/BLAS.h> // RK addition
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h> // log, M_PI
#include <limits.h> // INT_MIN, INT_MAX


double get_from_dmat(double* mat, int i, int j, int w) {
	return mat[i + j * w];
}

void set_dmat(double* mat, int i, int j, int w, double val) {
	mat[i + j * w] = val;
}

int get_from_imat(int* mat, int i, int j, int w) {
	return mat[i + j * w];
}

void set_imat(int* mat, int i, int j, int w, int val) {
	mat[i + j * w] = val;
}

double get_from_3d_dmat(double* mat, int i, int j, int k, int w, int h) {
	return mat[i + (j * w) + (k * (w * h))];
}

void set_3d_dmat(double* mat, int i, int j, int k, int w, int h, double val) {
	mat[i + (j * w) + (k * (w * h))] = val;
}

int get_from_3d_imat(int* mat, int i, int j, int k, int w, int h) {
	return mat[i + (j * w) + (k * (w * h))];
}

void set_3d_imat(int* mat, int i, int j, int k, int w, int h, int val) {
	mat[i + (j * w) + (k * (w * h))] = val;
}

enum subset_index {
	subset_len = 0,
	subset_sum = 1,
	subset_ssq = 2,
};

enum allseg_type {
	allseg_mean = 0,
	allseg_meanvar = 1,
	allseg_bern = 2,
};


void circ_fill_allseg_norm(double* all_seg, int period_len, int* subset_mat, int meanvar) {
	int jend;
	int j_circ;
	int ssq, sum, len;
	double LOG_2_PI = log(2*M_PI);
	for (int i = 0; i < period_len; i++) {
		int ssq = 0, sum = 0, jend = period_len + i - 1, len = 0;
		for (int j = 0; j < jend; j++) {
			j_circ = j % period_len;
			if (j_circ == 0) { j_circ = period_len; }
			j_circ--;
			len += get_from_imat(subset_mat, j_circ, subset_len, 3);
			sum += get_from_imat(subset_mat, j_circ, subset_sum, 3);
			ssq += get_from_imat(subset_mat, j_circ, subset_ssq, 3);
			if (meanvar) {
				double sigsq = ssq - (sum * sum) / len;
				if (sigsq < 0) { sigsq = 0.0000000000000001; }

				set_dmat(all_seg, i, j_circ, period_len, (-0.5 * len * (LOG_2_PI + log(sigsq) - log(len) + 1)));
			}
			else {
				set_dmat(all_seg, i, j_circ, period_len, (-0.5 * (ssq - (sum * sum) / len)));
			}
		}
	}
}


void circ_fill_allseg_bern(double* all_seg, int period_len, int* subset_mat) {
	int jend, j_circ;
	int ssq, sum, len;
	double sl, sl_inv;
	double LOG_2_PI = log(2*M_PI);
	for (int i = 0; i < period_len; i++) {
		ssq = 0, sum = 0, jend = period_len + i - 1, len = 0;
		for (int j = 0; j < jend; j++) {
			j_circ = j % period_len;
			if (j_circ == 0) { j_circ = period_len; }
			j_circ--;
			len += get_from_imat(subset_mat, j_circ, subset_len, 2);
			sum += get_from_imat(subset_mat, j_circ, subset_sum, 2);
			if (sum == len) {sl = 1; sl_inv = 1E-323;}
			else if (sum == 0) {sl = 1E-323; sl_inv = 1;}
			else {sl = sum/len; sl_inv = 1 - sl;}
			set_dmat(all_seg, i, j_circ, period_len, (sum * log(sl)) + (len - sum) * log(sl_inv));
		}
	}
}

double* circ_allseg(double* data, double* all_seg, int data_width, int period_len, enum allseg_type type) {
	int* subset_mat = malloc(3 * period_len * sizeof(int));
	int sum, ssq, len;
	double e, e2;
	for (int i = 0; i < period_len; i++) {
		sum = 0, ssq = 0, len = 0;
		for (int j = 0; j < data_width; j++) {
			if ((e = get_from_dmat(data, j, 1, data_width)) == i) {
				len++;
				e2 = get_from_dmat(data, j, 2, data_width);
				sum += e2;
				ssq += (e2 * e2);
			}
		}
		set_imat(subset_mat, i, subset_len, 3, len);
		set_imat(subset_mat, i, subset_sum, 3, sum);
		set_imat(subset_mat, i, subset_ssq, 3, ssq);
	}

	// double* all_seg = malloc(period_len * period_len * sizeof(double));
	switch (type) {
		case allseg_mean:
			circ_fill_allseg_norm(all_seg, period_len, subset_mat, 0);
			break;
		case allseg_meanvar:
			circ_fill_allseg_norm(all_seg, period_len, subset_mat, 1);
			break;
		case allseg_bern:
			circ_fill_allseg_bern(all_seg, period_len, subset_mat);
			break;
	}
	return all_seg;
}

int mod_index(int index, int mod) {
	int out = index % mod;
	if (out == 0) {out = mod;}
	out--;
	return out;
}

void sncirc(int* dist,
			   double* data,
			   int* data_len,
			   int* period_len,
			   int* max_cpts,
			   int* minseglen,
			   double* pen,
			   int* cptsout,
			   int* error,
			   int* cp,
			   int* cps_m,
			   int* op_cps,
			   int* op_k,
			   int* f_cpts,
			   int* criterion,
			   double* like_m,
			   double* like_m_coll,
			   int* lv,
			   double* all_seg) {

	double* all_seg = circ_allseg(data, all_seg, data_len, *period_len, *dist);

	// double* like_m = malloc((*max_cpts) * (*period_len) * (*period_len) * sizeof(double));
	for (int k = 0; k < *period_len; k++) {
		for (int j = 0; j < *period_len; j++) {
			set_3d_dmat(like_m, 1, j, k, (*max_cpts), (*period_len), get_from_dmat(all_seg, k, j, (*period_len)));
		}
	}

	double like = 0;
	int like_index;
	double tmp = 0;
	// int* cp = malloc((*max_cpts) * (*period_len) * (*period_len) * sizeof(int));
	for (int k = 0; k < (*period_len); k++) {
		for (int m = 1; m < (*max_cpts); m++) {
			for (int j = k + m * (*minseglen); j < (k - 1 + (*period_len)); j++) {
				like = -INFINITY;
				like_index = -1;
				int j_circ = mod_index(j, *period_len);
				for (int v = (k + (m - 1) * (*minseglen)); v < (j - (*minseglen)); v++) {
					int v_circ = mod_index(v, *period_len);
					int v2 = mod_index((v + 1), *period_len);
					tmp = get_from_3d_dmat(like_m, m - 1, v_circ, k, (*max_cpts), (*period_len)) + get_from_dmat(all_seg, v2, j_circ, *period_len);
					if (tmp > like) {
						like = tmp;
						like_index = v;
					}
				}
				set_3d_dmat(like_m, m, j_circ, k, (*max_cpts), (*period_len), like);
				set_3d_imat(cp, m, j_circ, k, (*max_cpts), (*period_len), like_index + (k + (m - 1) * (*minseglen) - 1));
			}
		}
	}

	// double* like_m_coll = malloc((*max_cpts) * (*period_len) * sizeof(double));
	// int* op_k = malloc((*max_cpts) * sizeof(int));
	double max_like_k;
	int max_like_k_index;
	for (int m = 0; m < (*max_cpts); m++) {
		max_like_k = -INFINITY;
		max_like_k_index = -1;
		for (int k = 0; k < (*period_len); k++) {
			int wrap = (k - 1 + (*period_len)) % (*period_len);
			if (wrap == 0) { wrap = ((*period_len) - 1); }
			double val = get_from_3d_dmat(like_m, m, wrap, k, (*max_cpts), (*period_len));
			if (val > max_like_k) {
				max_like_k = val;
				max_like_k_index = k;
			}
		}
		int k_opt = max_like_k_index;
		for (int i = 0; i < (*period_len); i++) {
			set_dmat(like_m_coll, m, i, (*max_cpts), get_from_3d_dmat(like_m, m, i, k_opt, (*max_cpts), (*period_len)));
		}
		op_k[m] = k_opt;
	}




	/// cps_m = m * m
	/// f_cpts = m
	f_cpts[0] = mod_index((op_k[0] - 1), (*period_len));
	for (int m = 1; m < (*max_cpts); m++) {
		f_cpts[m] = mod_index(op_k[m] - 1, (*period_len));
		set_imat(cps_m, m, 1, (*max_cpts), get_from_3d_imat(cp, m, f_cpts[m], op_k[m], (*max_cpts), (*period_len)));
		for (int i = 0; i < m-1; i++) {
			set_imat(cps_m, m, i+1, (*max_cpts), get_from_3d_imat(cp, m, get_from_imat(cps_m, m, i, (*max_cpts)), op_k[m], (*max_cpts), (*period_len)));
		}
	}

	*min_criterion = INT_MAX;
	*op_cps = -1;
	int tmp_crit = 0;
	int k_end;
	for (int i = 0; i < (*max_cpts); i++) {
		k_end = mod_index((op_k[i] - 1) + (*period_len), (*period_len));
		tmp_crit = -2 * get_from_dmat(like_m_coll, i, k_end, *(max_cpts)) + i * (*pen);
		if (tmp_crit < *min_criterion) {
			*min_criterion = tmp_crit;
			*op_cps = i;
		}
	}
	if (*op_cps == -1) {
		Rprintf("Error: op.cps not set\n");
		*error = 3;
		return;
	}

	if (*op_cps == (*max_cpts)) {
		Rprintf("The number of segments indentified is Q, it is advised to increase Q to make sure changepoints are not missed.\n");
	}

	if (*op_cps == 0) {
		cptsout[0] = (*period_len);
	}
	else {
		int index = 0;
		for (int i = 0; i < *(max_cpts); i++) {
			if (get_from_imat(cps_m, *op_cps, i, *(max_cpts)) > 0) {
				cptsout[index++] = get_from_imat(cps_m, *op_cps, i, *(max_cpts));
			}
		}
		cptsout[index++] = *(period_len);
		/// cptsout is sorted in R
	}


	free(all_seg);
	free(cp);
}
