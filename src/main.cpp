#include <cstdio>
#include <cstdlib>
#include <ctime>

#include "io.hpp"
#include "sais.hpp"
#include "lcp.hpp"
#include "overlapper.hpp"

#define FRAGMENT_THR 10
#define EDIST_THR 0.35

#define time_diff(s,f) (((f) - (s)) / (double) CLOCKS_PER_SEC)

int main(int argc, char **argv) {
	clock_t start, finish;

	read_data *reads;

	start = clock();
	int reads_cnt = parse_reads_readsim(&reads, argv[1]);
	finish = clock(); 
	printf("Parsing reads: %.4f sec\n", time_diff(start, finish));

	reads_sequence reads_seq;

	start = clock();
	int reads_seq_len = generate_reads_sequence(&reads_seq, reads, reads_cnt);
	finish = clock();
	printf("Generating reads sequence: %.4f sec\n", time_diff(start, finish));

	int *SA = (int*) malloc(reads_seq_len * sizeof(int));

	start = clock();
	sa_is(SA, reads_seq.sequence, reads_seq_len, reads_cnt + 4);
	finish = clock(); 
	printf("Generating SA: %.4f sec\n", time_diff(start, finish));

	int *LCP = (int*) malloc(reads_seq_len * sizeof(int));

	start = clock();
	lcp(LCP, SA, reads_seq.sequence, reads_seq_len);
	finish = clock();
	printf("Generating LCP: %.4f sec\n", time_diff(start, finish));

	fragments_map fragments;

	start = clock();
	preprocess_overlaps(&fragments, reads_seq, reads_seq_len, SA, LCP, FRAGMENT_THR);
	finish = clock();
	printf("Preprocessing overlaps: %.4f sec\n", time_diff(start, finish));

	free_reads_sequence(&reads_seq);
	free(SA);
	free(LCP);

	overlaps_map final_overlaps;

	start = clock();
	get_final_overlaps(&final_overlaps, fragments, reads, EDIST_THR);
	finish = clock();
	printf("Getting final overlaps: %.4f sec\n", time_diff(start, finish));

	printf("Saving final overlaps...\n");
	save_final_overlaps(final_overlaps, argv[2]);

	overlaps_map true_overlaps;

	start = clock();
	get_true_overlaps(&true_overlaps, reads, reads_cnt);
	finish = clock();
	printf("Getting true overlaps: %.4f sec\n", time_diff(start, finish));

	printf("\nEvaluating...\n");
	evaluate(final_overlaps, true_overlaps, reads_cnt);

	free_read_data(&reads, reads_cnt);

	return 0;
}
