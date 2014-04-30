C  cyclic_iof.h -- header file for cyclic array I/O functions
C 
C  See Chap 8, pp. 158 & ff in PPMPI
C  
      integer MAX_ENTRIES
      integer LOCAL_ENTRIES
      integer PROC_MAX
      parameter (MAX_ENTRIES = 1024, LOCAL_ENTRIES = 1024,
     +          PROC_MAX = 64)
C
C  Create a pseudo-structure for cylic array
C     integer  CYCLIC_ARRAY_T(8 + 1024 + 1024)   
      integer  arr_comm, arr_p, arr_myrank, global_order, padded_size,
     +         local_size, arr_stride, type, entry, local_entry
      parameter (arr_comm = 0,         arr_p = 1,       
     +           arr_myrank = 2,       global_order = 3,
     +           padded_size = 4,      local_size = 5,
     +           arr_stride = 6,       type = 7,
     +           entry = 8,
     +           local_entry=1031)
C   
 
