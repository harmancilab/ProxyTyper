#include "prxytypr_utils.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <vector>
#include "prxytypr_ansi_string.h"
#include "prxytypr_rng.h"
#include "prxytypr_seed_manager.h"
#include "prxytypr_vector_macros.h"

#if defined(_WIN32) || defined(WIN32) || defined(__WIN32__) 
	#include <Windows.h>
#endif

#if defined(unix) || defined(__unix__)  || defined(linux) || defined(__linux__) 
	#include <sys/types.h>
	#include <sys/stat.h>
	#include <unistd.h>

	#include <zlib.h>
	#include <limits.h> /* for PATH_MAX */
#endif

using namespace std;

bool __DUMP_UTIL_MESSAGES__ = false;

short read_bin_short(FILE* f_bin)
{
	short val;
	if (fread(&(val), sizeof(short), 1, f_bin) != 1)
	{
		fprintf(stderr, "Could not read integer @ %s(%d)\n", __FILE__, __LINE__);
		exit(1);
	}

	return(val);
}

unsigned short read_bin_ushort(FILE* f_bin)
{
	unsigned short val;
	if (fread(&(val), sizeof(unsigned short), 1, f_bin) != 1)
	{
		fprintf(stderr, "Could not read integer @ %s(%d)\n", __FILE__, __LINE__);
		exit(1);
	}

	return(val);
}

unsigned char read_bin_uchar(FILE* f_bin)
{
	unsigned char val;
	if (fread(&(val), sizeof(unsigned char), 1, f_bin) != 1)
	{
		fprintf(stderr, "Could not read integer @ %s(%d)\n", __FILE__, __LINE__);
		exit(1);
	}

	return(val);
}

char read_bin_char(FILE* f_bin)
{
	char val;
	if (fread(&(val), sizeof(char), 1, f_bin) != 1)
	{
		fprintf(stderr, "Could not read integer @ %s(%d)\n", __FILE__, __LINE__);
		exit(1);
	}

	return(val);
}

int read_bin_int(FILE* f_bin)
{
	int val;
	if (fread(&(val), sizeof(int), 1, f_bin) != 1)
	{
		fprintf(stderr, "Could not read integer @ %s(%d)\n", __FILE__, __LINE__);
		exit(1);
	}

	return(val);
}

double read_bin_double(FILE* f_bin)
{
	double val = -1;
	if (fread(&(val), sizeof(double), 1, f_bin) != 1)
	{
		fprintf(stderr, "Could not read integer @ %s(%d)\n", __FILE__, __LINE__);
		exit(1);
	}

	return(val);
}

void read_bin_int_array(FILE* f_bin, int* buff, int len)
{
	if (fread(buff, sizeof(int), len, f_bin) != (size_t)len)
	{
		fprintf(stderr, "Could not read int array of length %d @ %s(%d)\n", len, __FILE__, __LINE__);
		exit(1);
	}
}

void read_bin_char_array(FILE* f_bin, char* buff, int len)
{
	if (fread(buff, sizeof(char), len, f_bin) != (size_t)len)
	{
		fprintf(stderr, "Could not read char string of length %d @ %s(%d)\n", len, __FILE__, __LINE__);
		exit(1);
	}
}

void read_bin_uchar_array(FILE* f_bin, unsigned char* buff, int len)
{
	if (fread(buff, sizeof(unsigned char), len, f_bin) != (size_t)len)
	{
		fprintf(stderr, "Could not read char string of length %d @ %s(%d)\n", len, __FILE__, __LINE__);
		exit(1);
	}
}

void read_bin_uint_array(FILE* f_bin, unsigned int* buff, int len)
{
	if (fread(buff, sizeof(unsigned int), len, f_bin) != (size_t)len)
	{
		fprintf(stderr, "Could not read integer array of length %d @ %s(%d)\n", len, __FILE__, __LINE__);
		exit(1);
	}
}

void read_bin_double_array(FILE* f_bin, double* buff, int len)
{
	if (fread(buff, sizeof(double), len, f_bin) != (size_t)len)
	{
		fprintf(stderr, "Could not read double array of length %d @ %s(%d)\n", len, __FILE__, __LINE__);
		exit(1);
	}
}

void delete_file(const char* fp)
{
#ifdef __unix__
	unlink(fp);
#endif

#ifdef _WIN32
	remove(fp);
#endif
}

int compressFile(const char* inFile, const char * const outFileName)
{
#ifdef __unix__
	FILE *in = open_f(inFile, "rb");

	char buf[BUFSIZ] = { 0 };
	size_t bytes_read = 0;
	gzFile out = gzopen(outFileName, "wb");
	if (!out)
	{
		/* Handle error */
		fprintf(stderr, "Unable to open %s for writing\n", outFileName);
		return -1;
	}
	bytes_read = fread(buf, 1, BUFSIZ, in);
	while (bytes_read > 0)
	{
		int bytes_written = gzwrite(out, buf, bytes_read);
		if (bytes_written == 0)
		{
			int err_no = 0;
			fprintf(stderr, "Error during compression: %s", gzerror(out, &err_no));
			gzclose(out);
			return -1;
		}
		bytes_read = fread(buf, 1, BUFSIZ, in);
	}
	gzclose(out);

	fclose(in);
#endif 
	return 0;
}

void concatenateGzipFiles(const char* outputFile, vector<char*>* gzip_files)
{
	FILE* outFile = fopen(outputFile, "wb");
	if (outFile == NULL) {
		perror("Failed to open output file");
		exit(EXIT_FAILURE);
	}

	for (int i = 0; i < vecsize(gzip_files); ++i)
	{
		FILE* inFile = fopen(gzip_files->at(i), "rb");
		if (inFile == NULL) {
			perror("Failed to open input file");
			fclose(outFile);
			exit(EXIT_FAILURE);
		}

		char buffer[BUFSIZ];
		size_t bytesRead;
		while ((bytesRead = fread(buffer, 1, sizeof(buffer), inFile)) > 0) {
			fwrite(buffer, 1, bytesRead, outFile);
		}

		fclose(inFile);
	}

	fclose(outFile);
}

int compressSaveBuffer(const unsigned char* outBuffer, const size_t n_bytes_2_write, const char* outFileName)
{
	//size_t bytes_read = 0;
	gzFile out = gzopen(outFileName, "wb");
	if (!out)
	{
		/* Handle error */
		fprintf(stderr, "Unable to open %s for writing\n", outFileName);
		return 1;
	}

	int bytes_written = gzwrite(out, outBuffer, n_bytes_2_write);
	if (bytes_written == 0)
	{
		int err_no = 0;
		fprintf(stderr, "Error during compression: %s", gzerror(out, &err_no));
		gzclose(out);
		return 1;
	}
	
	gzclose(out);

	return(0);
}

#define CHUNK (16*1024*1024)

// Function to decompress and read a gzip file
//int decompressGzipFile(const char* source, unsigned char* outBuffer, size_t outSize) {
int decompressGzipFile(const char* source, unsigned char* outBuffer, size_t& n_bytes_read) 
{
	gzFile file = gzopen(source, "rb");
	if (!file) 
	{
		perror("gzopen");
		return -1;
	}

	unsigned char* buffer = new unsigned char[CHUNK];
	size_t totalSize = 0;
	//size_t bufferSize = CHUNK;
	//unsigned char* data = (unsigned char*)malloc(bufferSize);
	//unsigned char* data = new unsigned chat[bufferSize];
	unsigned char* data = outBuffer;

	if (data == NULL) 
	{
		perror("malloc");
		gzclose(file);
		return -1;
	}

	int bytesRead;
	while ((bytesRead = gzread(file, buffer, CHUNK)) > 0) 
	{
		memcpy(data + totalSize, buffer, bytesRead);
		totalSize += bytesRead;
	}

	if (bytesRead < 0) 
	{
		int err;
		const char* error_string = gzerror(file, &err);
		if (err) 
		{
			fprintf(stderr, "Error reading from file: %s\n", error_string);
			free(data);
			gzclose(file);
			return -1;
		}
	}

	gzclose(file);

	n_bytes_read = totalSize;

	delete[] buffer;

	//*outBuffer = data;
	//*outSize = totalSize;

	return 0;
}

char* get_file_extension(char* fp)
{
	t_string_tokens* tokens = t_string::tokenize_by_chars(fp, "/\\.");
	if(tokens->size() > 0)
	{
		char* fn = t_string::copy_me_str(tokens->back()->str());
		t_string::clean_tokens(tokens);
		return(fn);
	}
	else
	{
		return(NULL);
	}
}

bool does_file_exist(char* path)
{
#ifdef unix
	struct stat stat_var;
	int ret = stat(path, &stat_var);

	if(ret == 0)
	{
		return(true);
	}
	else
	{
		return(false);
	}
#endif 

#ifdef _WIN32
	fprintf(stderr, "no stat() implemented in Windows.\n");
	return(false);
#endif

	return(false);
}

void create_directory(char* directory_path)
{
#ifdef _WIN32
	CreateDirectory((LPCWSTR)directory_path, 0);
#endif

#ifdef unix
	mkdir(directory_path, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
#endif
}

void copy_file(char* src_fp, char* dest_fp)
{
	FILE* f_src = fopen(src_fp, "rb");
	FILE* f_dest = fopen(dest_fp, "wb");
	const size_t l_buff = (size_t)1024 * (size_t)1024;
	unsigned char buff[l_buff];
	while (1)
	{
		size_t n_read_bytes = fread(buff, 1, l_buff, f_src);
		fwrite(buff, 1, n_read_bytes, f_dest);
		if (n_read_bytes < l_buff)
		{
			break;
		}
	}
	fclose(f_dest);
	fclose(f_src);
}


void get_absolute_fp(char* fp, char* abs_fp)
{
	char cwd[1000];
	get_current_working_directory(cwd);

	if(fp[0] == '/')
	{
		strcpy(abs_fp, fp);
	}
	else
	{
		sprintf(abs_fp, "%s/%s", cwd, fp);
	}
}

void get_current_working_directory(char* buffer)
{
#ifdef _WIN32
	#error "Cannot compile getcwd.";
#endif

#ifdef unix
	char wd_buff[1000];
	char* ret_buff = getcwd(wd_buff, 1000);
	strcpy(buffer, ret_buff);
#endif
}

void remove_directory(char* directory_path)
{
}

char* get_file_name(char* fp)
{
	t_string_tokens* tokens = t_string::tokenize_by_chars(fp, "/\\");
	if(tokens->size() > 0)
	{
		char* fn = t_string::copy_me_str(tokens->back()->str());
		t_string::clean_tokens(tokens);
		return(fn);
	}
	else
	{
		return(NULL);
	}
}

char* get_directory_per_fp(char* fp)
{
	t_string_tokens* fp_tokens = t_string::tokenize_by_chars(fp, "\\/");

	if(fp_tokens->size() == 1)
	{
		t_string::clean_tokens(fp_tokens);
		return(NULL);
	}
	else
	{
		char* directory = new char[strlen(fp) + 2];
		directory[0] = 0;

		for(int i = 0; i < (int)fp_tokens->size() - 1; i++)
		{
			strcat(directory, fp_tokens->at(i)->str());
		} // i loop.

		t_string::clean_tokens(fp_tokens);

		return(directory);
	}
}

vector<char*>* load_directory_files(char* root_dir, char* extension)
{
        char ls_cmd[1000];

		if(extension != NULL)
		{
			sprintf(ls_cmd, "ls -d %s/*.%s | xargs -Ifiles basename files > cmd_op.txt", root_dir, extension);
			//int ret_val = system(ls_cmd);
		}
		else
		{
			sprintf(ls_cmd, "ls -d %s/* | xargs -Ifiles basename files > cmd_op.txt", root_dir);
			//int ret_val = system(ls_cmd);
		}

        vector<char*>* dir_files = new vector<char*>();

        FILE* f_cmd_op = fopen("cmd_op.txt", "r");
        char current_fn[1000];
        while(fscanf(f_cmd_op, "%s", current_fn) == 1)
        {
                char* new_fn = new char[strlen(current_fn) + 3];
                strcpy(new_fn, current_fn);
                dir_files->push_back(new_fn);
        }
        fclose(f_cmd_op);

		// Erase the file.
        sprintf(ls_cmd, "rm -f cmd_op.txt");
        //int ret_val = system(ls_cmd);

        return(dir_files);
}

char* load_header(char* fp)
{
	FILE* f = open_f(fp, "r");
	char* header_str = getline(f);
	fclose(f);

	return(header_str);
}

long int get_file_size(char* fp)
{
	FILE* f = NULL;

	if(strcmp(fp, "stdin") == 0 ||
		t_string::ends_with(fp, "gz"))
	{
		fprintf(stderr, "Cannot buffer stdin.\n");
		return(-1);
	}
	else
	{
		// Open the file if it is not the stdin.
		f = fopen(fp, "rb");
		if (f == NULL) 
		{ 
			fprintf(stderr, "Could not open %s\n", fp);
			return(-1);
		} 
	}

	fseek(f, 0, SEEK_END);
	long int file_size = ftell(f);
	fprintf(stderr, "Loading %s of %ld bytes.\n", fp, file_size);

	return(file_size);
}

// Loads a file into memory.
t_file_buffer* load_file(char* fp)
{
	FILE* f = NULL;

	// If the file is "stdin", load the standard input.
	if (strcmp(fp, "stdin") == 0 ||
		t_string::ends_with(fp, "gz"))
	{
		fprintf(stderr, "Cannot buffer stdin.\n");
		return(NULL);
	}
	else
	{
		// Open the file if it is not the stdin.
		f = fopen(fp, "rb");
		if (f == NULL) 
		{ 
			fprintf(stderr, "Could not open %s\n", fp);
			return(NULL);
		} 
	}

	fseek(f, 0, SEEK_END);
	long int file_size = ftell(f);
if(__DUMP_UTIL_MESSAGES__)
{
	fprintf(stderr, "Loading %s of %ld bytes.\n", fp, file_size);
}

	fseek(f, 0, SEEK_SET);
	char* file_buff = new char[file_size + 1];

if(__DUMP_UTIL_MESSAGES__)
{
	fprintf(stderr, "Allocated buffer memory.\n");
}

	if (file_size != (long int)fread(file_buff, sizeof(char), file_size, f)) 
	{ 
		delete [] file_buff;
		return NULL;
	}

if(__DUMP_UTIL_MESSAGES__)
{
	fprintf(stderr, "Read successfully.\n");
}

	fclose(f);
	file_buff[file_size] = 0;

	t_file_buffer* file_buffer = new t_file_buffer();
	file_buffer->file_buffer = file_buff;
	file_buffer->f_ptr = 0;
	file_buffer->l_file = file_size;

	return file_buffer;
}

void unload_file(t_file_buffer* file_buff)
{
	delete [] file_buff->file_buffer;
	delete file_buff;
}

int get_next_string(FILE* file, char* buff, int buff_size)
{
	memset(buff, 0, buff_size);

	int i_s = 0;
	int ret = 0;
	ret = getc(file);
	while(ret != ' ' && ret != EOF && ret != '\n' && ret != '\t')
	{		
		buff[i_s] = ret;
		i_s++;
		ret = getc(file);
	} // file reading loop.

	return(ret);
}

void save_lines(vector<char*>* lines, const char* fp)
{
	FILE* f_lines = open_f(fp, "w");
	for (int i_l = 0; i_l < vecsize(lines); i_l++)
	{
		fprintf(f_lines, "%s\n", lines->at(i_l));
	}
	close_f(f_lines, fp);
}

vector<char*>* buffer_file(const char* fp)
{
	vector<char*>* file_lines = new vector<char*>();
	//printf("Buffering %s.\n", fp);

	int n_lines = 0;
	FILE* f  = open_f(fp, "r");

	if(f == NULL)
	{
		//fprintf(stderr, "Could not open %s\n", fp);
		return(NULL);
	}

	while(1)
	{
		char* new_line = getline(f);
		if(new_line == NULL)
		{
			break;
		}
		file_lines->push_back(new_line);

		n_lines++;

if(__DUMP_UTIL_MESSAGES__)
{
		if((n_lines % 10000) == 0)
		{
			printf("Loaded %d lines.               \r", n_lines);
		}
}
	} // file buffering loop.

if(__DUMP_UTIL_MESSAGES__)
{
	fprintf(stderr, "\n");
}

	close_f(f, fp);

	return(file_lines);
}

bool check_file(const char* fp)
{
	FILE* f_temp = fopen(fp, "r");
	if(f_temp == NULL)
	{		
		return(false);
	}

	fclose(f_temp);
	return(true);
}

// Check for valid CR LF's depending on OS,
// Should be run for all ASCII input files.
void validate_file(char* fp)
{
#ifdef _WIN32
	char cur_char;
	// Open file in binary.
	FILE* f_ip_bin = open_f(fp, "rb");
	while(fread(&cur_char, 1, 1, f_ip_bin) == 1)
	{
		if(cur_char == CR)
		{
			if(fread(&cur_char, 1, 1, f_ip_bin) == 1)
			{
				if(cur_char != LF)
				{
					// Just a warning here.
					printf("%s is not compatible with dos ascii files. CR+LF problem at %s(%d).\n", fp, __FILE__, __LINE__);
					//exit(0);
				}
			}
			else
			{
				// Just a warning here.
				printf("%s is not compatible with dos ascii files. CR+LF problem at %s(%d).\n", fp, __FILE__, __LINE__);
				//exit(0);
			}
		}
		else if(cur_char == LF) // If there is an immediate LF before seeing a CR, this is a linux file.
		{
			// Just a warning here.
			printf("%s is not compatible with dos ascii files. CR+LF problem at %s(%d).\n", fp, __FILE__, __LINE__);
			//exit(0);
		}

	}
	fclose(f_ip_bin);
#endif

#ifdef __unix__
	char cur_char;
	// Open file in binary.
	FILE* f_ip_bin = open_f(fp, "rb");
	while(fread(&cur_char, 1, 1, f_ip_bin) == 1)
	{
		// Linux files do not contain CR's.
		// They only contain LF's.
		if(cur_char == CR)
		{
			// Just a warning here.
			printf("%s is not compatible with Linux ascii files. CR+LF problem at %s(%d).\n", fp, __FILE__, __LINE__);
			//exit(0);
		}
	}
	fclose(f_ip_bin);
#endif

#ifdef __APPLE__
	char cur_char;

	// Open file in binary.
	FILE* f_ip_bin = open_f(fp, "rb");
	while(fread(&cur_char, 1, 1, f_ip_bin) == 1)
	{
		// MAC files do not contain LF's.
		if(cur_char == LF)
		{
			// Just a warning here.
			printf("%s is not compatible with Linux ascii files. CR+LF problem at %s(%d).\n", fp, __FILE__, __LINE__);
			//exit(0);
		}
	}
	fclose(f_ip_bin);
#endif

}

FILE* open_f(const char* fp, const char* mode)
{
	if (fp == NULL || mode == NULL)
	{
		printf("Invalid arguments to open_f: NULL\n");
		exit(1);
	}

	FILE* f = NULL;
	if (strcmp(fp, "stdin") == 0)
	{
		f = stdin;
	}
	else if (strcmp(fp, "stdout") == 0)
	{
		f = stdout;
	}
	else if (strcmp(fp, "stderr") == 0)
	{
		f = stderr;
	}
	else if (t_string::ends_with(fp, ".bsam") ||
		t_string::ends_with(fp, ".bsam.gz"))
	{
		if (mode[0] == 'r')
		{
			char sbam_view_cmd[1000];
			sprintf(sbam_view_cmd, "mapped_read_tools -convert_bSAM_2_SAM_tagged \"%s\" stdout", fp);
#ifdef _WIN32
			f = _popen(sbam_view_cmd, "r");
#else 
			f = popen(sbam_view_cmd, "r");
#endif
		}
	}
	else if (t_string::ends_with(fp, ".bam"))
	{
		if (mode[0] == 'r')
		{
			char samtools_view_cmd[1000];
			sprintf(samtools_view_cmd, "samtools view -h \"%s\"", fp);
#ifdef _WIN32
			f = _popen(samtools_view_cmd, "r");
#else 
			f = popen(samtools_view_cmd, "r");
#endif
		}
	}
	else if (t_string::ends_with(fp, ".gz"))
	{
		if (mode[0] == 'r')
		{
			char ungzip_cmd[1000];
			sprintf(ungzip_cmd, "gzip -cd \"%s\"", fp);
#ifdef _WIN32
			f = _popen(ungzip_cmd, "r");
#else 
			f = popen(ungzip_cmd, "r");
#endif
		}
		else if (mode[0] == 'w')
		{
			char gzip_cmd[1000];
			sprintf(gzip_cmd, "gzip - -c -f | tee \"%s\" > /dev/null", fp);
#ifdef _WIN32
			f = _popen(gzip_cmd, "w");
#else 
			f = popen(gzip_cmd, "w");
#endif
		}
		else if (mode[0] == 'a')
		{
			char gzip_cmd[1000];
			sprintf(gzip_cmd, "gzip - -c -f | tee -a \"%s\" > /dev/null", fp);
#ifdef _WIN32
			f = _popen(gzip_cmd, "w");
#else 
			f = popen(gzip_cmd, "w");
#endif
		}
	}
	else
	{
		f = fopen(fp, mode);
	}

	if (f == NULL)
	{
		if (mode[0] == 'r')
		{
			printf("Could not open %s for reading.\n", fp);
			exit(1);
		}
		else if (mode[0] == 'w')
		{
			printf("Could not open %s for writing.\n", fp);
			exit(1);
		}
		else
		{
			printf("Could not open %s for requested operation.\n", fp);
			exit(1);
		}
	}

	return(f);
}

void close_f(FILE* f, const char* fp)
{
	if (fp == NULL)
	{
		fclose(f);
		return;
	}

	if (strcmp(fp, "stdin") == 0)
	{
	}
	else if (strcmp(fp, "stdout") == 0)
	{
	}
	else if (strcmp(fp, "stderr") == 0)
	{
	}
	else if (t_string::ends_with(fp, "gz"))
	{
#ifdef _WIN32
		_pclose(f);
#else 
		pclose(f);
#endif
	}
	else
	{
		fclose(f);
	}
}

int get_n_non_empty_lines(char* fp)
{
	int n_lines = 0;
	FILE* f = open_f(fp, "r");
	while(1)
	{
		char* cur_line = getline(f);
		if(cur_line == NULL)
		{
			break;
		}

		if(!line_empty(cur_line))
		{
			n_lines++;
		}

		delete [] cur_line;
	}
	fclose(f);

	return(n_lines);
}

void subsample_file_lines_no_buffer(char* fp, int n_lines_to_subsample)
{
	fprintf(stderr, "Counting the number of lines.\n");
	int n_lines = get_n_non_empty_lines(fp);
	fprintf(stderr, "%d lines in the file, subsampling %d lines.\n", n_lines, n_lines_to_subsample);

	//double p_select = (double)n_lines_to_subsample/n_lines;

	int n_sampled_lines = 0;
	//t_rng* rng = new t_rng(t_seed_manager::seed_me());

	// Open the file.
	FILE* f_sub = open_f("subsampled.txt", "w");
	FILE* f = open_f(fp, "r");

	// Loop to sample required number of lines.
	while(n_sampled_lines != n_lines_to_subsample)
	{
		fseek(f, 0, SEEK_SET);

		while(1)
		{
			char* cur_line = getline(f);
			if(cur_line == NULL)
			{
				break;
			}
			
			// Make sure that the line is not empty.
			if(!line_empty(cur_line))
			{
				fprintf(f_sub, "%s\n", cur_line);

				n_sampled_lines++;

				if(n_sampled_lines % 10000 == 0)
				{
					fprintf(stderr, "Sampled %d. line.         \r", n_sampled_lines);
				}
			}

			delete [] cur_line;
		} // file reading loop.
	} // sampling loop.

	fclose(f);
	fclose(f_sub);
}

bool line_empty(char* line)
{
	int l = strlen(line);
	for(int i = 0; i < l; i++)
	{
		if(line[i] != ' ' &&
			line[i] != '\t' &&
			line[i] != '\n')
		{
			return(false);
		}
	}

	return(true);
}

char* getline_per_file_buffer(t_file_buffer* file_buffer)
{
	//t_string* line_str = new t_string();	
	int i_buff = 0;
	int l_buff = 500;
	char* cur_buff = new char[l_buff];
	memset(cur_buff, 0, l_buff);

	char ret = 0;
	while(1)
	{
		//ret = getc(file);
		char cur_char;
		ret = get_next_char_per_file_buffer(file_buffer, cur_char);

		if(ret == false)
		{
			break;
		}
		
#ifdef _WIN32
		if(cur_char == CR)
		{
			// Make sure that the CR is followed by the correct character on the current system.
			ret = get_next_char_per_file_buffer(file_buffer, cur_char);

			if(ret == false)
			{
				fprintf(stderr, "Could not read a character after CR.\n");
				exit(1);
			}
			else if(cur_char != LF)
			{
				fprintf(stderr, "CR is not followed by LF.\n");
				exit(1);
			}

			break;
		}
#elif defined(__unix__)
		// There is nothing to check in UNIX since a carriage return is a new line.
		if(cur_char == CR)
		{
			fprintf(stderr, "Encountered CR character in __unix__.\n");
			exit(1);
		}
		else if(cur_char == LF)
		{
			// Unix has LF as new line. 
			break;
		}
#elif defined(__APPLE__)
		if(cur_char == LF)
		{
			fprintf(stderr, "Encountered CR character in __unix__.\n");
			exit(1);
		}
		else if(cur_char == CR)
		{
			// Unix has LF as new line. 
			break;
		}		
#else
#error "Neither _WIN32 nor __unix__ nor __APPLE__ is not defined.\n";
#endif

		// Update the buffer if it became too long.
		if(i_buff > (l_buff - 5))
		{
			l_buff *= 2;
			char* new_buff = new char[l_buff];
			memset(new_buff, 0, l_buff);
			strcpy(new_buff, cur_buff);
			delete [] cur_buff;
			cur_buff = new_buff;
		}

		cur_buff[i_buff] = cur_char;
		i_buff++;
	} // file reading loop.

	// If the end-of-file is encountered and there was nothing read, return NULL to indicate the there is nothing left to read and EOF is reached.
	if(ret == false && 
		//line_str->length() == 0)
		i_buff == 0)
	{
		//delete(line_str);
		delete [] cur_buff;
		return(NULL);
	}

	return(cur_buff);	
}

bool get_next_char_per_file_buffer(t_file_buffer* file_buffer, char& char_val)
{
	if(file_buffer->f_ptr < file_buffer->l_file)
	{
		char_val = file_buffer->file_buffer[file_buffer->f_ptr];
		file_buffer->f_ptr++;
		return(true);
	}
	else
	{
		char_val = 0;
		return(false);
	}
}

int get_n_lines(FILE* file)
{
	int n_lines = 0;

	char ret = 0;
	int n_chars_in_cur_line = 0;
	while(1)
	{
		ret = getc(file);

		if(ret == EOF)
		{
			break;
		}
		else if(ret == '\n')
		{
			if(n_chars_in_cur_line)
				n_lines++;

			n_chars_in_cur_line = 0;
		}
		else
		{
			n_chars_in_cur_line = 1;
		}
	} // file reading loop.

	return(n_lines);
}

bool get_next_token_per_file_till_newline(FILE* file, char* buffer, int l_buffer, const char* delim_char, bool& new_line_check, bool& eof_check)
{
	int i_buff = 0;
	char ret = 0;
	new_line_check = false;
	eof_check = false;
	while(1)
	{
		ret = getc(file);

		// If we found the end-of-file, break, then check whether there is anything in the buffer.
		if(ret == EOF)
		{
			eof_check = true;
			break;
		}

		// If we found a new line, break. The file pointer points to the new line's beginning.
		if(ret == '\n')
		{
			new_line_check = true;
			break;
		}

		// The current character is not a new line, is it the delimiter?
		if(ret == delim_char[0])
		{
			break;
		}
		else if(i_buff < l_buffer)
		{
			// Copy the read character to the buffer, if the buffer did not overflow, otherwise does not copy but keeps reading the file till the next token is found.
			buffer[i_buff] = ret;
			i_buff++;
		}
	} // file reading loop.

	if(i_buff < l_buffer)
	{
		buffer[i_buff] = 0;
	}
	else
	{
		buffer[l_buffer-1] = 0;
	}
	
	if(i_buff == 0)
	{
		return(false);
	}
	else
	{
		return(true);
	}
}

char* getline(FILE* file)
{
	vector<char>* line_vec = new vector<char>();

	char ret = 0;
	while(1)
	{
		ret = getc(file);

		if(ret == EOF)
		{
			break;
		}
		else if(ret == '\n')
		{
			break;
		}

		// Add this character to the vector of characters.
		line_vec->push_back(ret);
	} // file reading loop.

	// If the end-of-file is encountered and there was nothing read, return NULL to indicate the there is nothing left to read and EOF is reached.
	if(ret == EOF && 
		line_vec->size() == 0)
	{
		delete(line_vec);
		return(NULL);
	}

	char* buff = new char[line_vec->size() + 2];
	memset(buff, 0, sizeof(char) * (line_vec->size() + 2));

	// Copy the vector to buffer.
	for(int i = 0; i < (int)line_vec->size(); i++)
	{
		buff[i] = line_vec->at(i);
	} // i loop.

	delete(line_vec);

	return(buff);
}

char* resolve_data_dir()
{
	// try to resolve the DATAPATH_ENV_VAR.
	char* data_dir_from_env = getenv(DATAPATH_ENV_VAR);

	if(data_dir_from_env != NULL)
	{
		char* data_dir = t_string::copy_me_str(data_dir_from_env);
		return(data_dir);
	}
	else
	{
		char* data_dir = t_string::copy_me_str(LOCAL_DATA_PATH);
		return(data_dir);
	}

	printf("Could not resolve thermodynamics data directory.\n");
	exit(1);
}

char* x_fgets(char* buff, int size, FILE* file)
{
	if(fgets(buff, size, file) == NULL)
	{
		return(NULL);
	}

	if(buff[strlen(buff) - 1] == '\n')
	{
		int i_new_line_char = (int)strlen(buff) - 1;
		buff[i_new_line_char] = 0;
	}

	return(buff);
}

