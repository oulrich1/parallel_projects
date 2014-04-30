#include "util.h"



FILE*
file_create(char* filename, char* mode);


FILE*
file_write(FILE* ifp, char* buffer);


FILE*
file_read(FILE* ifp, char* buffer);


void
file_kill(FILE* fp);


int
file_check(FILE* fp);


int
file_exist(char* filename);

