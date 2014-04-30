# Searches tree on a network of SGI's running the mpich implementation
# of MPI

F77       =  f77
# Additional FFLAGS: -DSTATS, -DDEBUG
FFLAGS   =  -g -fullwarn -DSTATS
LDFLAGS  =
INCLUDE  =  -I/usr/local/mpich/include -I.
LIB      =  -L/usr/local/mpich/lib/IRIX/ch_p4 -lmpi

SOURCE = \
        main.f \
	par_tree_search.f \
	par_dfs.f \
	service_requests.f \
	work_remains.f \
	node_stack.f \
	solution.f \
	queue.f \
	terminate.f \
	stats.f \
	cio.f

OBJS = \
        main.o \
	par_tree_search.o \
	par_dfs.o \
	service_requests.o \
	work_remains.o \
	node_stack.o \
	solution.o \
	queue.o \
	terminate.o \
        stats.o \
	cio.o


tree: $(OBJS)
	$(F77) -o tree $(OBJS) $(INCLUDE) $(LIB)

clean:
	rm -f tree *.o core

main.o: ciof.h mainf.h node_stackf.h par_tree_searchf.h solutionf.h \
        terminatef.h statsf.h queuef.h

par_tree_search.o: ciof.h par_tree_searchf.h node_stackf.h par_dfsf.h \
        service_requestsf.h work_remainsf.h solutionf.h statsf.h

par_dfs.o: par_dfsf.h node_stackf.h queuef.h solutionf.h statsf.h

service_requests.o: service_requestsf.h node_stackf.h queuef.h terminatef.h \
	statsf.h

work_remains.o: work_remainsf.h node_stackf.h terminatef.h queuef.h \
	service_requestsf.h statsf.h

terminate.o: terminatef.h

node_stack.o: node_stackf.h

solution.o: solutionf.h ciof.h

queue.o: queuef.h node_stackf.h

stats.o: statsf.h

cio.o: ciof.h 

.f.o:
	$(F77) -c $(CFLAGS) $*.c $(INCLUDE)
