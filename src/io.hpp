#ifndef _IO_HPP_
#define _IO_HPP_

typedef struct {
	int id;
	char *name;
	unsigned char *sequence;
	unsigned char *reversed;
	int length;

	int start;
	int end;
} read_data;

typedef struct {
	int *sequence;
	int *read_index;
	int *read_pos;
} reads_sequence;

int parse_reads_readsim(read_data **reads, const char *filename);
void free_read_data(read_data **reads, int reads_cnt);

int generate_reads_sequence(reads_sequence *reads_seq, const read_data *reads, int reads_cnt);
void free_reads_sequence(reads_sequence *reads_seq);

#ifdef _IO_CPP_

#define READ_LEN_THR 15

#define BUF_SIZE 4096

static int read_line(char **line, FILE **fp, int flag);

#endif

#endif
