#define _SAIS_CPP_

#include <cstdlib>

#include "sais.hpp"

#define is_lms(i) ((i) > 0 && t[(i)-1] && !t[(i)])

void sa_is(int *SA, const int *s, int n, int k) {
	int i, j;

	char *t = (char*) malloc(n * sizeof(char));
	init_type(t, s, n);

	int *count = (int*) malloc(k * sizeof(int));
	init_count(count, k, s, n);

	int *bucket = (int*) malloc(k * sizeof(int));

	int lms_end = sort_lms_substrings(SA, bucket, count, k, s, t, n);
	int lms_end_name = name_lms_substrings(SA, lms_end, s, t, n);

	int *SA1 = SA;
	int *s1 = SA + n - lms_end;

	if (lms_end_name < lms_end - 1) {
		sa_is(SA1, s1, lms_end, lms_end_name + 1);
	} else {
		for (i = 0; i < lms_end; i++) SA1[s1[i]] = i;
	}

	for (i = 1, j = 0; i < n; i++) {
		if (is_lms(i)) s1[j++] = i;
	}

	for (i = 0; i < lms_end; i++) SA[i] = s1[SA1[i]];
	for (i = lms_end; i < n; i++) SA[i] = -1;

	init_bucket_end(bucket, count, k);

	for (i = lms_end - 1; i >= 0; i--) {
		j = SA[i];
		SA[i] = -1;
		SA[--bucket[s[j]]] = j;
	}

	induce_sa_l(SA, bucket, count, k, s, t, n);
	induce_sa_s(SA, bucket, count, k, s, t, n);

	free(t);
	free(count);
	free(bucket);
}

static void init_type(char *t, const int *s, int n) {
	t[n-1] = 0;
	t[n-2] = 1;

	for (int i = n - 3; i >= 0; i--) {
		t[i] = s[i] == s[i+1] ? t[i+1] : s[i] > s[i+1];
	}
}

static void init_count(int *count, int k, const int *s, int n) {
	for (int i = 0; i < k; i++) count[i] = 0;
	for (int i = 0; i < n; i++) count[s[i]]++;
}

static void init_bucket_start(int *bucket, const int *count, int k) {
	int sum = 0;

	for (int i = 0; i < k; i++) {
		bucket[i] = sum;
		sum += count[i];
	}
}

static void init_bucket_end(int *bucket, const int *count, int k) {
	int sum = 0;

	for (int i = 0; i < k; i++) {
		sum += count[i];
		bucket[i] = sum;
	}
}

static void induce_sa_l(int *SA, int *bucket, const int *count, int k, const int *s, const char *t, int n) {
	init_bucket_start(bucket, count, k);

	for (int i = 0; i < n; i++) {
		int j = SA[i];

		if (j > 0) {
			j--;
			if (t[j]) SA[bucket[s[j]]++] = j;
		}
	}
}

static void induce_sa_s(int *SA, int *bucket, const int *count, int k, const int *s, const char *t, int n) {
	init_bucket_end(bucket, count, k);

	for (int i = n - 1; i >= 0; i--) {
		int j = SA[i];

		if (j > 0) {
			j--;
			if (!t[j]) SA[--bucket[s[j]]] = j;
		}
	}
}

static int sort_lms_substrings(int *SA, int *bucket, const int *count, int k, const int *s, const char *t, int n) {
	int i;

	init_bucket_end(bucket, count, k);

	for (i = 0; i < n; i++) SA[i] = -1;

	for (i = 1; i < n; i++) {
		if (is_lms(i)) SA[--bucket[s[i]]] = i;
	}

	induce_sa_l(SA, bucket, count, k, s, t, n);
	induce_sa_s(SA, bucket, count, k, s, t, n);

	int lms_end = 0;

	for (i = 0; i < n; i++) {
		if (is_lms(SA[i])) SA[lms_end++] = SA[i];
	}

	for (i = lms_end; i < n; i++) SA[i] = -1;

	return lms_end;
}

static int name_lms_substrings(int *SA, int lms_end, const int *s, const char *t, int n) {
	int i, j;

	int prev = SA[0];
	SA[lms_end+prev/2] = 0;

	int lms_end_name = 0;

	for (i = 1; i < lms_end; i++) {
		int curr = SA[i];

		for (j = 0; ; j++) {
			if (s[curr+j] != s[prev+j] || t[curr+j] != t[prev+j]) {
				lms_end_name++; prev = curr;
				break;
			}

			if (j > 0 && (is_lms(curr+j) || is_lms(prev+j))) {
				break;
			}
		}

		SA[lms_end+curr/2] = lms_end_name;
	}

	for (i = n - 1, j = n - 1; i >= lms_end; i--) {
		if (SA[i] >= 0) SA[j--] = SA[i];
	}

	return lms_end_name;
}
