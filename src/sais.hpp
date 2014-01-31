#ifndef _SAIS_HPP_
#define _SAIS_HPP_

void sa_is(int *SA, const int *s, int n, int k);

#ifdef _SAIS_CPP_

static void init_type(char *t, const int *s, int n);

static void init_count(int *count, int k, const int *s, int n);
static void init_bucket_start(int *bucket, const int *count, int k);
static void init_bucket_end(int *bucket, const int *count, int k);

static void induce_sa_l(int *SA, int *bucket, const int *count, int k, const int *s, const char *t, int n);
static void induce_sa_s(int *SA, int *bucket, const int *count, int k, const int *s, const char *t, int n);

static int sort_lms_substrings(int *SA, int *bucket, const int *count, int k, const int *s, const char *t, int n);
static int name_lms_substrings(int *SA, int lms_end, const int *s, const char *t, int n);

#endif

#endif
