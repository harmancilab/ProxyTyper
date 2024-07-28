#ifndef _UTILS_
#define _UTILS_

#include <stdio.h>
#include <vector>

#define CR (0x0D)
#define LF (0x0A)

using namespace std;

#define __MAX_PATH (10000)

#define DATAPATH_ENV_VAR "DATAPATH"
#define LOCAL_DATA_PATH "data"

struct t_file_buffer
{
	char* file_buffer;
	long int f_ptr;
	long int l_file;
};

class t_string;

int compressFile(const char* inFile, const char * const outFileName);
int compressSaveBuffer(const unsigned char* outBuffer, const size_t n_bytes_2_write, const char* outFileName);
int decompressGzipFile(const char* source, unsigned char* outBuffer, size_t& n_bytes_read);
void concatenateGzipFiles(const char* outputFile, vector<char*>* gzip_files);

char* load_header(char* fp);

void delete_file(const char* fp);

bool does_file_exist(char* path);
void get_current_working_directory(char* buffer);
void create_directory(char* directory_path);
void remove_directory(char* directory_path);

void get_absolute_fp(char* fp, char* abs_fp);

void copy_file(char* src_fp, char* dest_fp);
void read_bin_double_array(FILE* f_bin, double* buff, int len);
char read_bin_char(FILE* f_bin);
int read_bin_int(FILE* f_bin);
unsigned char read_bin_uchar(FILE* f_bin);
short read_bin_short(FILE* f_bin);
unsigned short read_bin_ushort(FILE* f_bin);

double read_bin_double(FILE* f_bin);
void read_bin_int_array(FILE* f_bin, int* buff, int len);
void read_bin_char_array(FILE* f_bin, char* buff, int len);
void read_bin_uchar_array(FILE* f_bin, unsigned char* buff, int len);
void read_bin_uint_array(FILE* f_bin, unsigned int* buff, int len);

char* resolve_data_dir();
bool check_file(const char* fp);
void validate_file(char* fp);
char* x_fgets(char* buff, int size, FILE* file);
FILE* open_f(const char* fp, const char* mode);
void close_f(FILE* f, const char* fp);
vector<char*>* load_directory_files(char* root_dir, char* extension);
char* get_file_name(char* fp);
char* get_file_extension(char* fp);

long int get_file_size(char* fp);
void save_lines(vector<char*>* lines, const char* fp);
vector<char*>* buffer_file(const char* fp);
t_file_buffer* load_file(char* fp);
void unload_file(t_file_buffer* file_buff);
int get_next_string(FILE* file, char* buff, int buff_size);

// Loads a line from a file, without the buffer size. In principle, reads line of any length and buffers it.
char* getline(FILE* file);
bool get_next_token_per_file_till_newline(FILE* file, char* buffer, int l_buffer, const char* delim_char, bool& new_line_check, bool& eof_check);
char* getline_per_file_buffer(t_file_buffer* file_buffer);

bool get_next_char_per_file_buffer(t_file_buffer* file_buffer, char& char_val);

int get_n_lines(FILE* file);
bool line_empty(char* line);

int get_n_non_empty_lines(char* fp);
void subsample_file_lines_no_buffer(char* fp, int n_lines_to_subsample);

#endif // _UTILS_

