C  statsf.h -- definitions and declarations for stats.f
C
C
      integer MAX_TESTS
      integer MEMBERS
      integer INT_MEMBERS
      integer DOUBLE_MEMBERS 
      parameter (MAX_TESTS = 100, MEMBERS = 9, INT_MEMBERS = 6,
     +           DOUBLE_MEMBERS = 3)
C
C typedef struct  for statsi
       integer       nodes_expanded 
       integer       requests_sent
       integer       rejects_sent 
       integer       work_sent 
       integer       rejects_recd 
       integer       work_recd
C typedef struct for statsd 
       integer       par_dfs_time 
       integer       svc_req_time 
       integer       work_rem_time 
       parameter (nodes_expanded = 0, requests_sent = 1,
     +    rejects_sent = 2, work_sent = 3, rejects_recd = 4,
     +    work_recd = 5,  par_dfs_time = 0, svc_req_time = 1,
     +    work_rem_time = 2)
C
       integer          statsi
       double precision overhead_time
       double precision statsd
       double precision start_time
       common /StatInts/  statsi(0:5)
       common /StatReals/  overhead_time, statsd(0:2),
     +          start_time
C

