C 
C  service_requestsf.h
C 
C  Definitions and declarations used by the service_requests functions 
C
         integer send_stack_mpi_t
        common /SndStkMpi/ send_stack_mpi_t
C
         integer  allocated_size
         integer  block_lengths 
         integer  displacements 
         integer  node_list, node_list_size
         parameter (node_list_size = 100)
         integer  node_count
         common /Svc/ allocated_size, block_lengths(0:200), 
     +          displacements(0:200), node_list(0:200),
     +          node_count
C
C  functions called in Service_requests   
       integer   Allocate_type_arrays 
       integer   Nodes_available 
C
