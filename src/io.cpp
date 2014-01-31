#define _IO_CPP_

#include <cstdio>
#include <cstdlib>
#include <cstring>

#include "io.hpp"

int parse_reads_readsim(read_data **reads, const char *filename) {
	FILE *fp = fopen(filename, "r");

	if (fp == NULL) {
		perror("fopen");
		exit(-1);
	}

	int reads_cnt = 0;

	while (fscanf(fp, "%*s\n%*s\n") != EOF) reads_cnt++;

	rewind(fp);

	*reads = (read_data*) malloc(reads_cnt * sizeof(read_data));

	int too_short_cnt = 0;

	for (int i = 0; i < reads_cnt; i++) {
		if (fgetc(fp) != '>') {
			fprintf(stderr, "Invalid readsim file: %s\n", filename);
			exit(-1);
		}

		char *name;
		read_line(&name, &fp, 0);

		int start, end;
		if (sscanf(name, "%*[^:]:%d-%d:F\n", &start, &end) != 2) {
			fprintf(stderr, "Invalid readsim file: %s\n", filename);
			exit(-1);
		}

		char *sequence;
		int length = read_line(&sequence, &fp, 1);

		char *reversed = (char*) malloc(length * sizeof(char));
		for (int j = length - 1, k = 0; j >= 0; j--) {
			reversed[k++] = sequence[j];
		}

		if (length < READ_LEN_THR) {
			free(name);
			free(sequence);
			too_short_cnt++;

			continue;
		}

		(*reads)[i-too_short_cnt].id = i;
		(*reads)[i-too_short_cnt].name = name;
		(*reads)[i-too_short_cnt].sequence = (unsigned char*) sequence;
		(*reads)[i-too_short_cnt].reversed = (unsigned char*) reversed;
		(*reads)[i-too_short_cnt].length = length;
		(*reads)[i-too_short_cnt].start = start;
		(*reads)[i-too_short_cnt].end = end;
	}

	fclose(fp);

	reads_cnt -= too_short_cnt;

	*reads = (read_data*) realloc(*reads, reads_cnt * sizeof(read_data));

	return reads_cnt;
}

void free_read_data(read_data **reads, int reads_cnt) {
	for (int i = 0; i < reads_cnt; i++) {
		free((*reads)[i].name);
		free((*reads)[i].sequence);
		free((*reads)[i].reversed);
	}

	free(*reads);
}

int generate_reads_sequence(reads_sequence *reads_seq, const read_data *reads, int reads_cnt) {
	int reads_len = 0;

	for (int i = 0; i < reads_cnt; i++) {
		reads_len += reads[i].length;
	}

	int seq_len = reads_len + reads_cnt;

	reads_seq->sequence = (int*) malloc(seq_len * sizeof(int));
	reads_seq->read_index = (int*) malloc(seq_len * sizeof(int));
	reads_seq->read_pos = (int*) malloc(seq_len * sizeof(int));

	int pos = 0;

	for (int i = 0; i < reads_cnt; i++) {
		for (int j = 0; j < reads[i].length; j++) {
			(reads_seq->sequence)[pos] = reads[i].sequence[j] + reads_cnt;
			(reads_seq->read_index)[pos] = i;
			(reads_seq->read_pos)[pos] = j;

			pos++;
		}

		(reads_seq->sequence)[pos] = reads_cnt - i - 1;
		(reads_seq->read_index)[pos] = i;
		(reads_seq->read_pos)[pos] = reads[i].length;

		pos++;
	}

	return seq_len;
}

void free_reads_sequence(reads_sequence *reads_seq) {
	free(reads_seq->sequence);
	free(reads_seq->read_index);
	free(reads_seq->read_pos);
}

static int read_line(char **line, FILE **fp, int flag) {
	*line = NULL;

	int length = 0;

	do {
		*line = (char*) realloc(*line, (length + BUF_SIZE) * sizeof(char));
		fgets(*line + length, BUF_SIZE, *fp);

		length = strlen(*line);
	} while ((*line)[length-1] != '\n' && !feof(*fp));

	length -= 1;

	if (!flag) {
		*line = (char*) realloc(*line, (length + 1) * sizeof(char));
		(*line)[length] = '\0';
	} else {
		*line = (char*) realloc(*line, length * sizeof(char));

		for (int i = 0; i < length; i++) {
			switch ((*line)[i]) {
				case 'A': (*line)[i] = 0; break;
				case 'C': (*line)[i] = 1; break;
				case 'G': (*line)[i] = 2; break;
				case 'T': (*line)[i] = 3; break;
				default: fprintf(stderr, "Invalid sequence character %c", (*line)[i]); exit(-1);
			}
		}
	}

	return length;
}
