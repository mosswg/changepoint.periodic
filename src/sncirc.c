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


void circ_fill_allseg_norm(double* all_seg, int period_len, double* subset_mat, int meanvar) {
	int jend;
	int j_circ;
	double ssq, sum, len;
	double ssq_t, sum_t, len_t;
	double LOG_2_PI = log(2*M_PI);
	for (int i = 0; i < period_len; i++) {
		ssq = 0, sum = 0, jend = period_len + i - 1, len = 0;
		for (int j = i; j <= jend; j++) {
			j_circ = j % period_len;
			len_t = get_from_dmat(subset_mat, subset_len, j_circ + 1, 3);
			sum_t = get_from_dmat(subset_mat, subset_sum, j_circ + 1, 3);
			ssq_t = get_from_dmat(subset_mat, subset_ssq, j_circ + 1, 3);
			len += len_t;
			sum += sum_t;
			ssq += ssq_t;
			if (len == 0) {
				Rprintf("Error: len zero\n");
				continue;
			}
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


void circ_fill_allseg_bern(double* all_seg, int period_len, double* subset_mat) {
	int jend, j_circ;
	double sum, len;
	double sl, sl_inv;
	for (int i = 0; i < period_len; i++) {
		sum = 0, jend = period_len + i - 1, len = 0;
		for (int j = i; j <= jend; j++) {
			j_circ = j % period_len;
			len += get_from_dmat(subset_mat, subset_len, j_circ + 1, 3);
			sum += get_from_dmat(subset_mat, subset_sum, j_circ + 1, 3);
			if (sum == len) {sl = 1; sl_inv = 1E-323;}
			else if (sum == 0) {sl = 1E-323; sl_inv = 1;}
			else {sl = sum/len; sl_inv = 1 - sl;}
			set_dmat(all_seg, i, j_circ, period_len, (sum * log(sl)) + ((len - sum) * log(sl_inv)));
		}
	}
}

double* circ_allseg(double* data, double* all_seg, int data_len, int period_len, enum allseg_type type) {
	double* subset_mat = malloc(3 * period_len * sizeof(double));
	for (int i = 0; i < (3 * period_len); i++) {
		subset_mat[i] = 0.0f;
	}
	double sum, ssq, len;
	double e, e2;
	for (int i = 1; i <= period_len; i++) {
		sum = 0, ssq = 0, len = 0;
		for (int j = 0; j < data_len; j++) {
			if ((e = get_from_dmat(data, j, 0, data_len)) == i) {
				len++;
				e2 = get_from_dmat(data, j, 1, data_len);
				sum += e2;
				ssq += (e2 * e2);
			}
		}
		set_dmat(subset_mat, subset_len, i, 3, len);
		set_dmat(subset_mat, subset_sum, i, 3, sum);
		set_dmat(subset_mat, subset_ssq, i, 3, ssq);
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
			   double* criterion,
			   double* like_m,
			   double* like_m_coll,
			   int* lv,
			   double* all_seg,
			   int* op_like) {

	circ_allseg(data, all_seg, *data_len, *period_len, *dist);

	// double* like_m = malloc((*max_cpts) * (*period_len) * (*period_len) * sizeof(double));
	for (int k = 0; k < *period_len; k++) {
		for (int j = 0; j < *period_len; j++) {
			set_3d_dmat(like_m, 0, j, k, (*max_cpts), (*period_len), get_from_dmat(all_seg, k, j, (*period_len)));
		}
	}

	double like = 0;
	int like_index, count;
	double tmp = 0;
	// int* cp = malloc((*max_cpts) * (*period_len) * (*period_len) * sizeof(int));
	for (int k = 1; k <= (*period_len); k++) {
		for (int m = 2; m <= (*max_cpts); m++) {
			for (int j = k + m * (*minseglen); j <= (k - 1 + (*period_len)); j++) {
				like = -INFINITY;
				like_index = -1;
				count = 1;
				int j_circ = j % *period_len;
				if (j_circ == 0) { j_circ = *period_len; }
				for (int v = (k + (m - 1) * (*minseglen)); v <= (j - (*minseglen)); v++) {
					int v_circ = v % *period_len;
					if (v_circ == 0) {v_circ = *period_len;}
					int v2 = (v + 1) % *period_len;
					if (v2 == 0) {v2 = *period_len;}
					tmp = get_from_3d_dmat(like_m, m - 2, v_circ - 1, k - 1, (*max_cpts), (*period_len)) + get_from_dmat(all_seg, v2 - 1, j_circ - 1, *period_len);
					if (tmp > like) {
						like = tmp;
						like_index = count;
					}
					count++;
				}
				set_3d_dmat(like_m, m-1, j_circ-1, k-1, (*max_cpts), (*period_len), like);
				int set_val = (like_index + (k + (m - 1) * *minseglen - 1)) % (*period_len);
				if (set_val == 0) { set_val = (*period_len); }
				set_3d_imat(cp, m-1, j_circ-1, k-1, (*max_cpts), (*period_len), set_val);
			}
		}
	}

	// double* like_m_coll = malloc((*max_cpts) * (*period_len) * sizeof(double));
	// int* op_k = malloc((*max_cpts) * sizeof(int));
	double max_like_k;
	int max_like_k_index;
	for (int m = 1; m <= (*max_cpts); m++) {
		max_like_k = -INFINITY;
		max_like_k_index = -1;
		for (int k = 1; k <= (*period_len); k++) {
			int wrap = (k - 1 + (*period_len)) % (*period_len);
			if (wrap == 0) {wrap = (*period_len);}
			double val = get_from_3d_dmat(like_m, m-1, wrap-1, k-1, (*max_cpts), (*period_len));
			if (val > max_like_k) {
				max_like_k = val;
				max_like_k_index = k;
			}
		}
		for (int i = 0; i < (*period_len); i++) {
			set_dmat(like_m_coll, m-1, i, (*max_cpts), get_from_3d_dmat(like_m, m-1, i, max_like_k_index-1, (*max_cpts), (*period_len)));
		}
		op_k[m-1] = max_like_k_index;
	}

	/// cps_m = m * m
	/// f_cpts = m
	f_cpts[0] = (op_k[0] - 1) % (*period_len);
	if (f_cpts[0] == 0) { f_cpts[0] = (*period_len); }
	for (int m = 2; m <= (*max_cpts); m++) {
		int f = (op_k[m - 1] - 1) % (*period_len);
		if (f == 0) { f = (*period_len); }
		f_cpts[m-1] = f;
		set_imat(cps_m, m-1, 0, (*max_cpts), get_from_3d_imat(cp, m-1, f-1, op_k[m-1]-1, (*max_cpts), (*period_len)));
		for (int i = 1; i <= m-1; i++) {
			set_imat(cps_m, m-1, i, (*max_cpts), get_from_3d_imat(cp, m-i-1, get_from_imat(cps_m, m-1, i-1, (*max_cpts))-1, op_k[m-1]-1, (*max_cpts), (*period_len)));
		}
	}

	*op_cps = -1;

	double min_criterion = INFINITY;
	int min_criterion_index = 0;
	int k_end;
	for (int i = 0; i < (*max_cpts); i++) {
		k_end = mod_index((op_k[i] - 1) + (*period_len), (*period_len));
		lv[i] = get_from_dmat(like_m_coll, i, k_end, *(max_cpts));
		criterion[i] = -2 * lv[i]  + i * (*pen);
		if (criterion[i] < min_criterion) {
			min_criterion = criterion[i];
			min_criterion_index = i;
		}
	}
	if (min_criterion_index == 0) {
		*op_cps = 0;
	}
	else {
		*op_cps = min_criterion_index + 1;
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

	if (*op_cps == 0) {
		*op_like = criterion[0];
	}
	else {
		*op_like = criterion[*op_cps - 1];
	}
}
