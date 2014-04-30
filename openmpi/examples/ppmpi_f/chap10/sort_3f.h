C  sort_3.h -- header file for sort_3.f 
C  3. Add prototype for Key_compare
C  3. Add macro for LIST_BUF_SIZE and MAX_KEY_STRING
C
      integer KEY_MIN
      integer KEY_MAX
      integer KEY_MOD
      integer LIST_BUF_SIZE
      integer MAX_KEY_STRING 
      parameter (KEY_MIN = 0, KEY_MAX = 32767,
     +           KEY_MOD = 32769 ,
     +           LIST_BUF_SIZE = 128, MAX_KEY_STRING = 10)
C
      integer list_size, allocated_size, keys
      parameter (list_size = 1, allocated_size = 0,
     +           keys = 2)
C
C functions
       integer Get_list_size 
       integer Allocate_list 
       integer Key_compare
