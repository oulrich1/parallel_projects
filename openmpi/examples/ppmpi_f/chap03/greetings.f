c  greetings.f -- greetings program
c 
c  Send a message from all processes with rank != 0 to process 0.
c     Process 0 prints the messages received.
c 
c  Input: none.
c  Output: contents of messages received by process 0.
c
c  Note:  Due to the differences in character data in Fortran and char
c      in C, their may be problems in MPI_Send/MPI_Recv
c 
c  See Chapter 3, pp. 41 & ff in PPMPI.
c 
	program greetings
c
	include 'mpif.h'
c
	integer my_rank
	integer p
	integer source
	integer dest
	integer tag
	character*100 message
	character*10 digit_string
	integer size
	integer status(MPI_STATUS_SIZE)
	integer ierr
c
c   function
	integer string_len
c
	call MPI_Init(ierr)
c
	call MPI_Comm_rank(MPI_COMM_WORLD, my_rank, ierr)
	call MPI_Comm_size(MPI_COMM_WORLD, p, ierr)
c
	if (my_rank.ne.0) then
	    call to_string(my_rank, digit_string, size)
	    message = 'Greetings from process ' // digit_string(1:size)
     +                 // '!'
	    dest = 0
	    tag = 0
	    call MPI_Send(message, string_len(message), MPI_CHARACTER,
     +                    dest, tag, MPI_COMM_WORLD, ierr)
	else
	    do 200 source = 1, p-1
                tag = 0
	        call MPI_Recv(message, 100, MPI_CHARACTER, source,
     +                        tag, MPI_COMM_WORLD, status, ierr)
                call MPI_Get_count(status, MPI_CHARACTER, size, ierr)
                write(6,100) message(1:size)
  100           format(' ',a)
  200	    continue
	endif
c
	call MPI_Finalize(ierr)
	end
c
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c   Converts the integer stored in number into an ascii 
c   string.  The string is returned in string.  The number of 
c   digits is returned in size.
	subroutine to_string(number, string, size)
	integer number
	character *(*) string
	integer size

	character*100 temp
	integer local
	integer last_digit
	integer i

	local = number
        i = 0

c   strip digits off starting with least significant
c   do-while loop
 100	    last_digit = mod(local,10)
	    local = local/10
	    i = i + 1
	    temp(i:i) = char(last_digit + ichar('0'))
	if (local.ne.0) go to 100
	
	size = i

c   reverse digits
	do 200 i = 1, size
	    string(size-i+1:size-i+1) = temp(i:i)
 200	continue
c
	return
	end
c   to_string
c
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c   Finds the number of characters stored in a string
c
	integer function string_len(string)
	character*(*) string
c
	character*1 space
	parameter (space = ' ')
	integer i
c
	i = len(string)

c   while loop
  100   if ((string(i:i).eq.space).and.(i.gt.1)) then
	    i = i - 1
            go to 100
	endif
c
	if ((i.eq.1).and.(string(i:i).eq.space)) then
	    string_len = 0
	else
	    string_len = i
	endif
c
	return
	end
c   end of string_len
