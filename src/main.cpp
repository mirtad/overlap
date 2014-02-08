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
	time_t start, finish;

	read_data *reads;

	printf("Parsing reads...\n");
	int reads_cnt = parse_reads_readsim(&reads, argv[1]);

	reads_sequence reads_seq;

	printf("Generating reads sequence...\n");
	int reads_seq_len = generate_reads_sequence(&reads_seq, reads, reads_cnt);

	int *SA = (int*) malloc(reads_seq_len * sizeof(int));

	time(&start);
	sa_is(SA, reads_seq.sequence, reads_seq_len, reads_cnt + 4);
	time(&finish); 
	printf("Generated SA: %.4f sec\n", difftime(finish, start));

	int *LCP = (int*) malloc(reads_seq_len * sizeof(int));

	time(&start);
	lcp(LCP, SA, reads_seq.sequence, reads_seq_len);
	time(&finish); 
	printf("Generated LCP: %.4f sec\n", difftime(finish, start));

	fragments_map fragments;

	time(&start);
	preprocess_overlaps(&fragments, reads_seq, reads_seq_len, SA, LCP, FRAGMENT_THR);
	time(&finish); 
	printf("Preprocessed overlaps: %.4f sec\n", difftime(finish, start));

	free_reads_sequence(&reads_seq);
	free(SA);
	free(LCP);

	overlaps_map final_overlaps;

	time(&start);
	get_final_overlaps(&final_overlaps, &fragments, reads, EDIST_THR, atoi(argv[3]));
	time(&finish); 
	printf("Got final overlaps: %.4f sec\n", difftime(finish, start));

	printf("Saving final overlaps...\n");
	save_final_overlaps(final_overlaps, argv[2]);

	overlaps_map true_overlaps;

	printf("Getting true overlaps...\n");
	get_true_overlaps(&true_overlaps, reads, reads_cnt);

	printf("\nEvaluating...\n");
	evaluate(final_overlaps, true_overlaps, reads_cnt);

	free_read_data(&reads, reads_cnt);

	return 0;
}
