#include <cstdio>
#include <cstdlib>
#include <ctime>

#include "io.hpp"
#include "sais.hpp"
#include "lcp.hpp"

int main(int argc, char **argv) {
	clock_t start, finish;

	read_data *reads;

	start = clock();
	int reads_cnt = parse_reads_readsim(&reads, argv[1]);
	finish = clock(); 
	printf("Parsing reads: %.4f sec\n", (double) (finish - start) / (double) CLOCKS_PER_SEC);

	reads_sequence reads_seq;

	start = clock();
	int reads_seq_len = generate_reads_sequence(&reads_seq, reads, reads_cnt);
	finish = clock();
	printf("Generating reads sequence: %.4f sec\n", (double) (finish - start) / (double) CLOCKS_PER_SEC);

	int *SA = (int*) malloc(reads_seq_len * sizeof(int));

	start = clock();
	sa_is(SA, reads_seq.sequence, reads_seq_len, reads_cnt + 4);
	finish = clock(); 
	printf("Generating SA: %.4f sec\n", (double) (finish - start) / (double) CLOCKS_PER_SEC);

	int *LCP = (int*) malloc(reads_seq_len * sizeof(int));

	start = clock();
	lcp(LCP, SA, reads_seq.sequence, reads_seq_len);
	finish = clock();
	printf("Generating LCP: %.4f sec\n", (double) (finish - start) / (double) CLOCKS_PER_SEC);

	printf("Checking SA and LCP...\n");

	for (int i = 1; i < reads_seq_len; i++) {
		for (int j = 0; j < LCP[i]; j++) {
			if (reads_seq.sequence[SA[i-1]+j] != reads_seq.sequence[SA[i]+j]) {
				printf("LCP error at position %d.\n", i);
			}
		}

		if (SA[i-1] + LCP[i] < reads_seq_len && SA[i] + LCP[i] < reads_seq_len) {
			if (reads_seq.sequence[SA[i-1]+LCP[i]] >= reads_seq.sequence[SA[i]+LCP[i]]) {
				printf("SA error at positions %d-%d.\n", i - 1, i);
			}
		}
	}

	printf("Done!\n");

	free_read_data(&reads, reads_cnt);

	free_reads_sequence(&reads_seq);

	free(SA);

	free(LCP);

	return 0;
}
