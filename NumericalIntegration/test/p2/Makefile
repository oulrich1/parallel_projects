CC  = mpicc
CXX = mpic++

ROOT_DIR = ../../

TARGET1	 = find_min_trapezoids
TARGET2  = parallel_integration

BIN_DIR  = $(ROOT_DIR)bin/

INC_DIR  = $(ROOT_DIR)include/

LIB_DIR  = $(ROOT_DIR)lib/
LIB_SRCS = $(wildcard $(LIB_DIR)/*.c)
LIB_OBJS = $(LIB_SRCS:.c=.o)

LIBS  = -lpthread -lm

CFLAGS  = -I$(INC_DIR) -g \
        -Wall -Wextra -Werror \
        -Wno-unused-result -Wno-missing-field-initializers -Wno-unused-parameter \
        -Wno-unused-but-set-parameter -Wno-unused-function \
        -std=c99 $(LIBS)

LDFLAGS = -L$(LIB_DIR)

CS_DIR  = .
CS_SRCS = $(wildcard $(CS_DIR)/*.c)
CS_OBJS = $(CS_SRCS:.c=.o)




DEPS_DIR = $(INC_DIR)
# _DEPS = approx_methods.h serial_integrations.h util.h
DEPS = $(wildcard $(DEPS_DIR)/*.h)

OBJ_DIR = obj/
_OBJ    = approx_methods.o functional_expression.o integral_operations.o util.o file_io.o
OBJ     = $(patsubst %,$(OBJ_DIR)%,$(_OBJ))  


# hopefully this makes two executables..
all: $(TARGET1) $(TARGET2)

$(TARGET1):	$(OBJ) $(OBJ_DIR)$(TARGET1).o
	$(CC) -o $@ $^ $(CFLAGS)

$(TARGET2): $(OBJ) $(OBJ_DIR)$(TARGET2).o
	$(CC) -o $@ $^ $(CFLAGS)
# $(OBJ_DIR)%.o: %.c 
# 	$(CC) -c -o $@ $< $(CFLAGS)

$(OBJ_DIR)$(TARGET1).o: $(TARGET1).c
	$(CC) -c -o $(OBJ_DIR)$(TARGET1).o $< $(CFLAGS)

$(OBJ_DIR)$(TARGET2).o: $(TARGET2).c
	$(CC) -c -o $(OBJ_DIR)$(TARGET2).o $< $(CFLAGS)

$(OBJ_DIR)functional_expression.o: $(LIB_DIR)functional_expression.c
	$(CC) -c -o $(OBJ_DIR)functional_expression.o $< $(CFLAGS)

$(OBJ_DIR)integral_operations.o: $(LIB_DIR)integral_operations.c
	$(CC) -c -o $(OBJ_DIR)integral_operations.o $< $(CFLAGS)

$(OBJ_DIR)approx_methods.o: $(LIB_DIR)approx_methods.c
	$(CC) -c -o $(OBJ_DIR)approx_methods.o $< $(CFLAGS)

$(OBJ_DIR)file_io.o: $(LIB_DIR)file_io.c
	$(CC) -c -o $(OBJ_DIR)file_io.o $< $(CFLAGS)

$(OBJ_DIR)util.o: $(LIB_DIR)util.c
	$(CC) -c -o $(OBJ_DIR)util.o $< $(CFLAGS)


.PHONY: clean

clean:
	rm -f $(OBJ_DIR)*.o 
	rm -f *~ $(INC_DIR)*~ 
	rm -f $(TARGET1)
	rm -f $(TARGET2)

# References 
# https://www.gnu.org/software/make/manual/make.html

# http://stackoverflow.com/questions/7004702/
# how-can-i-create-a-makefile-for-c-projects-
# with-src-obj-and-bin-subdirectories



# $(TARGET1): $(CS_OBJS) $(LIB_OBJS)
#     $(CC) $(LDFLAGS) -o $(BIN_DIR)/$(TARGET1) $(CS_OBJS) $(LIBS)

# %.o: %.c
#     $(CC) $(CFLAGS) -c $< -o $@
