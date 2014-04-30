
#include "file_io.h"

FILE*
file_create(char* filename, char* mode) {
    if(!filename)
        return NULL;
    FILE *ofp;

    char _mode[2];
    if(!mode) {
        memcpy(_mode, "a+", 2);
        mode = _mode;
    }
    ofp = fopen(filename, mode);

    if (!ofp) {
        fprintf(stdout, "Cannot open the output file..\n");
        exit(1);
    }
    return ofp;
}


FILE*
file_write(FILE* ofp, char* buffer) {
    if( !ofp ) 
        return NULL;
    if (buffer) {
        fprintf(ofp, "%s\n", buffer);
    }
    return ofp;
}


FILE*
file_read(FILE* ifp, char* buffer) {
    if( !ifp )
        return NULL;
    if (buffer) {
        fscanf(ifp, "%s\n", buffer);
    }
    return ifp;
}

void
file_kill(FILE* fp) {
    if(fp){
        fclose(fp);
    }
}

int
file_check(FILE* fp) {
    if (fp) {
        return 1;
    }
    return 0;
}

int
file_exist(char* filename) {
    if (filename) {
        return 1;
    }
    return 0;
}