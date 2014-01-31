#ifndef _OVERLAPPER_HPP_
#define _OVERLAPPER_HPP_

#include <map>
#include <unordered_map>
#include <boost/functional/hash.hpp>

#include "io.hpp"

typedef struct {
	int s1;
	int s2;
	int len;
} fragment;

typedef struct {
	int len1;
	int len2;
	int edist;
	const char *type;
} overlap;

typedef std::pair<int, int> int_pair;
typedef boost::hash<int_pair> int_pair_hash;

typedef std::unordered_map<int_pair, fragment, int_pair_hash> fragments_map;

typedef std::map<int_pair, overlap> overlaps_map;

void preprocess_overlaps(fragments_map *fragments, reads_sequence reads_seq, int reads_seq_len, const int *SA, const int *LCP, int threshold);
void get_final_overlaps(overlaps_map *overlaps, fragments_map fragments, const read_data *reads, float threshold);
void save_final_overlaps(overlaps_map overlaps, const char *filename);

void get_true_overlaps(overlaps_map *overlaps, const read_data *reads, int reads_cnt);

void evaluate(overlaps_map final_overlaps, overlaps_map true_overlaps, int reads_cnt);

#ifdef _OVERLAPPER_CPP_

#define OVERLAP_LEN_THR 15
#define K_LEN_THR 20

#endif

#endif
