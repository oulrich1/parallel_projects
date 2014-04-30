C  stats.f -- Functions for use in gathering performance information on
C      parallel tree search.  
C 
C  Note:  All times are local!
C 
C  See Chap 14, pp. 328 & ff, in PPMPI for a discussion of parallel tree
C      search.
C
C ******************************************************************  
      subroutine Get_overhead() 
      integer i
      include 'statsf.h'
C
      statsd(par_dfs_time) = 0.0D
      do 100 i = 0 , MAX_TESTS - 1
          call Start_time_rtn( )
          call Finish_time(par_dfs_time)
 100  continue
      overhead_time = statsd(par_dfs_time)/MAX_TESTS
      statsd(par_dfs_time) = 0.0D
      return
      end
C
C ******************************************************************  
      subroutine Set_stats_to_zero(statini,statind )
      integer statini(0:5)
      double precision    statind(0:2)
      include 'statsf.h'   
C
      statini(nodes_expanded) = 0
      statini(requests_sent) = 0
      statini(rejects_sent) = 0
      statini(work_sent) = 0
      statini(rejects_recd) = 0
      statini(work_recd) = 0
      statind(par_dfs_time) = 0.0D
      statind(svc_req_time) = 0.0D
      statind(work_rem_time) = 0.0D
      return
      end
C
C ******************************************************************  
      subroutine Update_totals(totalsi, totalsd, newi, newd)
      integer  totalsi(0:5), newi(0:5)   
      double precision totalsd(0:2) ,newd(0:2)
C
      include 'statsf.h'
C      
      totalsi(nodes_expanded) = totalsi(nodes_expanded) +
     +                          newi(nodes_expanded)
      totalsi(requests_sent) =  totalsi(requests_sent) +
     +                          newi(requests_sent)
      totalsi(rejects_sent) =   totalsi(rejects_sent) +
     +                          newi(rejects_sent)
      totalsi(work_sent) =      totalsi(work_sent) +
     +                          newi(work_sent)
      totalsi(rejects_recd) =   totalsi(rejects_recd) +
     +                          newi(rejects_recd)
      totalsi(work_recd) =      totalsi(work_recd) +
     +                           newi(work_recd)
      if (totalsd(par_dfs_time) .LT. newd(par_dfs_time)) then
          totalsd(par_dfs_time) = newd(par_dfs_time)
      endif
      if (totalsd(svc_req_time) .LT. newd(svc_req_time)) then
          totalsd(svc_req_time) = newd(svc_req_time)
      endif
      if (totalsd(work_rem_time) .LT. newd(work_rem_time)) then
          totalsd(work_rem_time) = newd(work_rem_time)
      endif
      return
      end 
C
C
C ******************************************************************  
      subroutine  Print_title() 
C
      print *,' '
      print *,'                Performance Statistics'
      print *,'     (Totals are sums for counts and maxima
     + for times)'
      print *,'      Nodes  Reqs  Rejs  Work  Rejs  Work
     +    DFS     Svc       Wk_rm'
      print *,'Proc   Exp   Sent  Sent  Sent  Recd  Recd
     +    Time    Time       Time'
      print *,'----  -----  ----  ----  ----  ----  ----
     +    ----    ----      -----'
      return
      end
C
C
C ******************************************************************  
      subroutine Print_ind_stats(rank, statini, statind)
      integer       rank
      integer       statini(0:5)
      double precision statind(0:2)
      include       'statsf.h' 
C  
      if (rank  .LT.  0) then
            print 10, 'TOT',
     +          statini(nodes_expanded),
     +          statini(requests_sent),
     +          statini(rejects_sent),
     +          statini(work_sent),
     +          statini(rejects_recd),
     +          statini(work_recd),
     +          statind(par_dfs_time),
     +          statind(svc_req_time),
     +          statind(work_rem_time)
10          format(A,4x,I5,1X,I5,1x,I5,1x,I5,1x,
     +       I5,1x,I5,6x,
     +       F7.4,3x,F7.4,3x,F7.4)
      else
           print 20, rank,
     +          statini(nodes_expanded),
     +          statini(requests_sent),
     +          statini(rejects_sent),
     +          statini(work_sent),
     +          statini(rejects_recd),
     +          statini(work_recd),
     +          statind(par_dfs_time),
     +          statind(svc_req_time),
     +          statind(work_rem_time)
20        format(I3,4x,I5,1X,I5,1x,I5,1x,I5,1x,
     +        I5,1x,I5,6x,
     +        F7.4,3x,F7.4,3x,F7.4)
      endif
      return
      end
C
C
C ******************************************************************  
      subroutine Print_stats(io_comm)
      integer io_comm
      include 'statsf.h'
      include 'mpif.h'
      include 'ciof.h'
      integer recd_statsi(0:5), totalsi(0:5)
      double precision recd_statsd(0:2),
     +                 totalsd(0:2)
      integer         p
      integer         io_rank
      integer         my_rank
      integer         q
      integer         status(MPI_STATUS_SIZE)
      integer         ierr
      integer         retval
      integer         DTAG, ITAG
      parameter       (DTAG=500,ITAG=600)
C
       call MPI_COMM_SIZE(io_comm,  p, ierr)
       call MPI_COMM_RANK(io_comm,  my_rank, ierr)
       retval = Get_io_rank(io_comm,  io_rank)
C
       call Set_stats_to_zero( totalsi, totalsd)
       call Set_stats_to_zero( recd_statsi, recd_statsd)
C
      if (my_rank .EQ. io_rank)  then
          call Print_title()
          do 100 q = 0  ,io_rank -1
              call MPI_RECV( recd_statsi, INT_MEMBERS,
     +            MPI_INTEGER, q, ITAG, io_comm,
     +            status, ierr)
              call MPI_RECV( recd_statsd, DOUBLE_MEMBERS,
     +            MPI_DOUBLE_PRECISION, q, DTAG, io_comm,
     +            status, ierr)
              call Print_ind_stats(q,  recd_statsi,
     +                                recd_statsd)
              call Update_totals( totalsi, totalsd,
     +                     ,  recd_statsi, recd_statsd)
 100      continue
          call Print_ind_stats(io_rank,  statsi, statsd)
          call Update_totals( totalsi,totalsd,
     +                         statsi, statsd)
          do 200 q = io_rank+1 , p-1   
              call MPI_RECV( recd_statsi, INT_MEMBERS, 
     +            MPI_INTEGER, q, ITAG, io_comm,
     +            status, ierr)
              call MPI_RECV( recd_statsd, DOUBLE_MEMBERS,
     +            MPI_DOUBLE_PRECISION, q, DTAG, io_comm,
     +            status, ierr)
              call Print_ind_stats(q,
     +                       recd_statsi,recd_statsd)
              call Update_totals( totalsi, totalsd,
     +                    recd_statsi, recd_statsd)
 200      continue
          call Print_ind_stats(-1,  totalsi, totalsd)
      else  
           call MPI_SEND( statsi, INT_MEMBERS, MPI_INTEGER, 
     +             io_rank, ITAG, io_comm, ierr)
           call MPI_SEND(statsd, DOUBLE_MEMBERS,
     +            MPI_DOUBLE_PRECISION, io_rank,
     +            DTAG,io_comm,ierr)
      endif
      return
      end
C****************************************************************
      subroutine Start_time_rtn() 
C
      include 'statsf.h'
      include 'mpif.h'
      integer ierr
      Start_time = MPI_Wtime(ierr)
      return
      end
C
C****************************************************************
      subroutine Finish_time(time) 
C
      integer time
      include 'statsf.h'
      include 'mpif.h'
      double precision finish, diff
      integer ierr
      finish = MPI_Wtime(ierr) 
      diff = finish - start_time - overhead_time
      statsd(time) = statsd(time) + diff
      return                  
      end
C
C****************************************************************
      subroutine  Incr_stat(member)
      integer member 
      include 'statsf.h'
C
      if (member .LE. 5) then
         statsi(member) = statsi(member) + 1
      endif
      return
      end
C
