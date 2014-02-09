#define _OVERLAPPER_CPP_

#include <cstdio>
#include <cstdlib>

#include <pthread.h>
#include <semaphore.h>
#include <fcntl.h>

#include "myers/myers.h"

#include "overlapper.hpp"

#define min2(x,y) ((x) < (y) ? (x) : (y))
#define max2(x,y) ((x) > (y) ? (x) : (y))

#define edist_k(len) ((len) < K_LEN_THR ? (len) : ((len) / 2))

void preprocess_overlaps(fragments_map *fragments, reads_sequence reads_seq, int reads_seq_len, const int *SA, const int *LCP, int threshold) {
		int max_lcp = LCP[0];

		for (int i = 1; i < reads_seq_len; i++) {
			if (LCP[i] > max_lcp) max_lcp = LCP[i];
		}

		for (int len = max_lcp; len >= threshold; len--) {
			int lbound = 1, ubound = 0;

			do {
				while (lbound < reads_seq_len && LCP[lbound] < len) lbound++;

				ubound = lbound;

				while (ubound < reads_seq_len && LCP[ubound] >= len) ubound++;

				for (int i = lbound - 1; i < ubound; i++) {
					for (int j = i + 1; j < ubound; j++) {
						int sa1 = SA[i];
						int sa2 = SA[j];

						if (reads_seq.read_index[sa1] == reads_seq.read_index[sa2])
							continue;

						if (sa1 > sa2) {
							int tmp = sa1;
							sa1 = sa2;
							sa2 = tmp;
						}

						int_pair indices = std::make_pair(reads_seq.read_index[sa1], reads_seq.read_index[sa2]);

						if (fragments->find(indices) == fragments->end()) {
							fragment f;

							f.s1 = reads_seq.read_pos[sa1];
							f.s2 = reads_seq.read_pos[sa2];
							f.len = len;

							fragments->emplace(indices, f);
						}
					}
				}

				lbound = ubound;
			} while (lbound < reads_seq_len);
		}
}

void get_final_overlaps(overlaps_map *overlaps, fragments_map *fragments, const read_data *reads, float threshold, int threads_no) {
	thread_data data;

	data.overlaps = overlaps;
	data.fragments = fragments;
	data.reads = reads;
	data.threshold = threshold;

	pthread_t threads[threads_no];

	for (int i = 0; i < threads_no; i++) {
		pthread_create(&threads[i], NULL, process_overlap, (void*) (&data));
	}

	for (int i = 0; i < threads_no; i++) {
		pthread_join(threads[i], NULL);
	}
}

void save_final_overlaps(overlaps_map overlaps, const char *filename) {
	FILE *fp = fopen(filename, "w");

	if (fp == NULL) {
		perror("fopen");
		exit(-1);
	}

	for (auto it = overlaps.begin(); it != overlaps.end(); it++) {
		int_pair ids = it->first;
		overlap o = it->second;

		fprintf(fp, "%d %d %d %d %s %d\n", ids.first, ids.second, o.len1, o.len2, o.type, o.edist);
	}

	fclose(fp);
}

void get_true_overlaps(overlaps_map *overlaps, const read_data *reads, int reads_cnt) {
	for (int i = 0; i < reads_cnt; i++) {
		for (int j = i + 1; j < reads_cnt; j++) {
			int s1 = reads[i].start;
			int e1 = reads[i].end;
			int s2 = reads[j].start;
			int e2 = reads[j].end;

			int len = min2(e1, e2) - max2(s1, s2) + 1;

			if (len >= OVERLAP_LEN_THR) {
				int_pair ids = std::make_pair(reads[i].id, reads[j].id);

				overlap o;

				o.len1 = len;
				o.len2 = len;
				o.edist = 0;

				if (s1 < s2) {
					o.type = "EB";
				} else {
					o.type = "BE";
				}

				overlaps->insert(std::make_pair(ids, o));
			}
		}
	}
}

void evaluate(overlaps_map final_overlaps, overlaps_map true_overlaps, int reads_cnt) {
	int tp = 0, fp = 0, fn = 0, tn = 0;

	for (auto it = final_overlaps.begin(); it != final_overlaps.end(); it++) {
		if (true_overlaps.count(it->first) == 0) {
			fp++;
		} else {
			tp++;
		}
	}

	for (auto it = true_overlaps.begin(); it != true_overlaps.end(); it++) {
		if (final_overlaps.count(it->first) == 0) {
			fn++;
		}
	}

	tn = (reads_cnt * (reads_cnt - 1)) / 2 - (fp + tp + fn);

	printf("Confusion matrix:\n%d %d\n%d %d\n", tp, fn, fp, tn);

	float accuracy = (tp + tn) / (float) (tp + fp + fn + tn);
	float precision = tp / (float) (tp + fp);
	float recall = tp / (float) (tp + fn);
	float f1 = 2 * precision * recall / (precision + recall);

	printf("Accuracy: %.2f%%\n", accuracy * 100);
	printf("Precision: %.2f%%\n", precision * 100);
	printf("Recall: %.2f%%\n", recall * 100);
	printf("F1 score: %.2f%%\n", f1 * 100);
}

static void* process_overlap(void *arg) {
	thread_data* data = (thread_data*) arg;

	overlaps_map *overlaps = data->overlaps;
	fragments_map *fragments = data->fragments;
	const read_data *reads = data->reads;
	float threshold = data->threshold;

	sem_t *fragments_sem = sem_open("fragments_sem", O_CREAT, S_IRWXU, 1);
	sem_t *overlaps_sem = sem_open("overlaps_sem", O_CREAT, S_IRWXU, 1);

	while (true) {
		sem_wait(fragments_sem);

			if (fragments->empty()) {
				sem_post(fragments_sem);

				sem_close(fragments_sem);
				sem_close(overlaps_sem);

				return NULL;
			}

			auto it = fragments->begin();

			int_pair indices = it->first;

			int ind1 = indices.first;
			int ind2 = indices.second;

			fragment f = it->second;

			int s1 = f.s1;
			int s2 = f.s2;
			int len = f.len;

			fragments->erase(it);

		sem_post(fragments_sem);

		int e1 = s1 + len;
		int e2 = s2 + len;

		overlap o;

		o.len1 = o.len2 = len;
		o.edist = 0;

		const unsigned char *seq1, *seq2;
		int len1, len2;
		int edist, end;

		len1 = reads[ind1].length - e1;
		len2 = reads[ind2].length - e2;

		if (len1 != 0 && len2 != 0) {
			seq1 = reads[ind1].sequence + e1;
			seq2 = reads[ind2].sequence + e2;

			if (len1 < len2) {
				myersCalcEditDistance(seq1, len1, seq2, len2, 4, edist_k(len1), MYERS_MODE_SHW, &edist, &end);
				o.len1 += len1;
				o.len2 += end + 1;
			} else {
				myersCalcEditDistance(seq2, len2, seq1, len1, 4, edist_k(len2), MYERS_MODE_SHW, &edist, &end);
				o.len1 += end + 1;
				o.len2 += len2;
			}

			if (edist == -1)
				continue;

			o.edist += edist;
		}

		len1 = s1;
		len2 = s2;

		if (len1 != 0 && len2 != 0) {
			seq1 = reads[ind1].reversed + reads[ind1].length - s1;
			seq2 = reads[ind2].reversed + reads[ind2].length - s2;

			if (len1 < len2) {
				myersCalcEditDistance(seq1, len1, seq2, len2, 4, edist_k(len1), MYERS_MODE_SHW, &edist, &end);
				o.len1 += len1;
				o.len2 += end + 1;

				if (end + 1 <= len2) {
					o.type = "EB";
				} else {
					o.type = "BE";
				}
			} else {
				myersCalcEditDistance(seq2, len2, seq1, len1, 4, edist_k(len2), MYERS_MODE_SHW, &edist, &end);
				o.len1 += end + 1;
				o.len2 += len2;

				if (end + 1 <= len1) {
					o.type = "EB";
				} else {
					o.type = "BE";
				}
			}

			if (edist == -1)
				continue;

			o.edist += edist;
		}

		if (o.edist / ((o.len1 + o.len2) / 2.0) < threshold) {
			int_pair ids = std::make_pair(reads[ind1].id, reads[ind2].id);

			sem_wait(overlaps_sem);

			overlaps->insert(std::make_pair(ids, o));

			sem_post(overlaps_sem);
		}
	}

	return NULL;
}
