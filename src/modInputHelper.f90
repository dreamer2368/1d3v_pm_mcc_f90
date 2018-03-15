module modInputHelper

  use, intrinsic :: iso_fortran_env, only : output_unit
  use constants

  implicit none

  interface getOption
     procedure getOptionInteger_
     procedure getOptionLogical_
     procedure getOptionScalar_
     procedure getOptionString_
  end interface getOption

  interface getRequiredOption
     procedure getRequiredOptionInteger_
     procedure getRequiredOptionLogical_
     procedure getRequiredOptionScalar_
     procedure getRequiredOptionString_
  end interface getRequiredOption

!  private :: getOptionInteger_, getOptionLogical_, getOptionScalar_, getOptionString_,       &
!       getRequiredOptionInteger_, getRequiredOptionLogical_, getRequiredOptionScalar_,       &
!       getRequiredOptionString_

  integer, parameter :: STRING_LENGTH = 256

  type t_DictElement
     character(len = STRING_LENGTH) :: key, val
  end type t_DictElement

  type(t_DictElement), allocatable :: dict(:)

contains

  function split(str, leftSubstring, rightSubstring, separator) result(extraSeparatorCount)

    ! <<< Arguments >>>
    character(len = *), intent(in) :: str
    character(len = *), intent(out) :: leftSubstring, rightSubstring
    character(len = 1), intent(in) :: separator

    ! <<< Result >>>
    integer :: extraSeparatorCount

    ! <<< Local variables >>>
    integer :: i, j

    i = 1

    do j = 1, len_trim(str)
       if (str(j:j) == separator) exit
       leftSubstring(i:i) = str(j:j)
       i = i + 1
    end do

    leftSubstring(i:) = " "
    rightSubstring = str(j+1:len_trim(str))
    call trimAll(leftSubstring)
    call trimAll(rightSubstring)

    extraSeparatorCount = scan(rightSubstring, separator)

  end function split

  subroutine trimAll(str)

    ! <<< Arguments >>>
    character(len = *), intent(inout) :: str

    ! <<< Local variables >>>
    integer :: i

    do i = 1, len(str)
       if (ichar(str(i:i)) > 32) exit
       if (ichar(str(i:i)) <= 32) str(i:i) = " "
    end do

    do i = len(str), 1, -1
       if (ichar(str(i:i)) > 32) exit
       if (ichar(str(i:i)) <= 32) str(i:i) = " "
    end do

    str = trim(adjustl(str))

  end subroutine trimAll

  subroutine sort()

    ! <<< Local variables >>>
    integer :: i, j
    type(t_DictElement) :: temp

    do i = 2, size(dict)
       j = i - 1
       temp = dict(i)
       do while (dict(j)%key > temp%key)
          dict(j+1) = dict(j)
          j = j - 1
          if (j < 1) exit
       end do
       dict(j+1) = temp
    end do

  end subroutine sort

  subroutine find(key, index)

    ! <<< Arguments >>>
    character(len = *), intent(in) :: key
    integer, intent(out) :: index

    ! <<< Local variables >>>
    integer :: iStart, iEnd, iMid

    index = -1

    if (allocated(dict)) then

       iStart = 1
       iEnd = size(dict)

       do while (iEnd >= iStart) !... `dict` is sorted; use binary search.
          iMid = (iStart + iEnd) / 2
          if (dict(iMid)%key == trim(key)) then
             index = iMid
             return
          end if
          if (dict(iMid)%key > trim(key)) then
             iEnd = iMid - 1
          else
             iStart = iMid + 1
          end if
       end do

    end if

  end subroutine find

function getFreeUnit(fileUnit) result(freeUnit)

  ! <<< Arguments >>>
  integer, intent(out), optional :: fileUnit

  ! <<< Result >>>
  integer :: freeUnit

  ! <<< Local variables >>>
  logical :: isOpened

  do freeUnit = 10, 10000
     inquire(unit = freeUnit, opened = isOpened)
     if (.not. isOpened) exit
  end do

  if (present(fileUnit)) fileUnit = freeUnit

end function getFreeUnit

subroutine parseInputFile(filename, commentMarker, separator)

  implicit none

  ! <<< Arguments >>>
  character(len = *), intent(in) :: filename
  character(len = 1), intent(in), optional :: commentMarker, separator

  ! <<< Local variables >>>
  character(len = 1) :: commentMarker_, separator_
  integer :: i, fileUnit, procRank, dictSize, lineNo, istat, ierror
  character(len = STRING_LENGTH) :: line, message

  if( .not. (len_trim(filename) > 0) ) then
    write(message,"(A)") "no filename??"
    write(output_unit,"(A)") message
    flush(output_unit)
    stop
  end if

  ! Assign default values to optional arguments.
  commentMarker_ = '#'
  if (present(commentMarker)) commentMarker_ = commentMarker
  separator_ = '='
  if (present(separator)) separator_ = separator

  ! Check if file exists.
  open(unit = getFreeUnit(fileUnit), file = trim(filename), action = 'read',              &
       status = 'old', iostat = istat)
  if (istat /= 0) then
     write(message, "(2A)") trim(filename), ": File not found or permission denied!"
    write(output_unit,"(A)") message
    flush(output_unit)
     stop
  end if

  ! Only the root process reads the input file.
  dictSize = 0
  do !... read once to find input dictionary size.
     read(fileUnit, '(A)', iostat = istat) line
     if (istat < 0) exit
     call stripComments(line, commentMarker_) !... skip comments.
     if (len_trim(line) == 0) cycle !... skip empty lines.
     dictSize = dictSize + 1
  end do
  close(fileUnit)

  ! Broadcast the input dictionary size to all processes.
  if (dictSize == 0) then
     write(message, "(2A)") trim(filename), ": File is empty or does not contain any input!"
    write(output_unit,"(A)") message
    flush(output_unit)
     stop
  end if

  ! Allocate memory to hold the dictionary.
  if ( allocated(dict) ) deallocate(dict)
  allocate(dict(dictSize), stat = istat)
  if (istat /= 0) then
     write(message, "(A)") "Insufficient memory: Could not allocate storage for input!"
    write(output_unit,"(A)") message
    flush(output_unit)
     stop
  end if

  ! Again, only the root process reads the input file.
  i = 0 ; lineNo = 0 ; istat = 0
  open(unit = getFreeUnit(fileUnit), file = trim(filename),                               &
       action = 'read', status = 'old')
  do !... read again to fill input dictionary.

     read(fileUnit, '(A)', iostat = istat) line
     if (istat < 0) then
        istat = 0
        exit
     end if
     lineNo = lineNo + 1 !... need line number for reporting errors.
    call stripComments(line, commentMarker_)
    if (len_trim(line) == 0) cycle
    i = i + 1

    ! Parse in 'key <separator> value' format.
    istat = split(line, dict(i)%key, dict(i)%val, separator_)

    if (istat /= 0) then
       write(message, "(2A,I0.0,A)") trim(filename), ": Failed to parse input on line ", &
            lineNo, "!"
       exit
    end if

    if (len_trim(dict(i)%key) == 0) then
       istat = 1
       write(message, "(2A,I0.0,A)") trim(filename), ": Empty parameter key on line ",   &
            lineNo, "!"
       exit
    end if

  end do
  close(fileUnit)

  ! Broadcast istat to collectively return on error.

  ! Check if an error occurred.
  if (istat /= 0) then
    print *, 'error'
    stop
  end if

  ! Sort the dictionary and broadcast it to all processes.
  call sort()

end subroutine parseInputFile

subroutine stripComments(str, commentMarker)

  ! <<< Arguments >>>
  character(len = *), intent(inout) :: str
  character(len = 1), intent(in) :: commentMarker

  ! <<< Local variables >>>
  integer :: i

  call trimAll(str)

  do i = 1, len(str)
     if (str(i:i) == commentMarker) then
        str(i:) = " "
        exit
     end if
  end do

end subroutine stripComments

function getOptionInteger_(key, defaultValue) result(val)

  ! <<< Arguments >>>
  character(len = *), intent(in) :: key
  integer, intent(in) :: defaultValue

  ! <<< Result >>>
  integer :: val

  ! <<< Local variables >>>
  integer :: index, stat

  call find(key, index)

  if (index == -1) then
     val = defaultValue
     return
  end if

  read(dict(index)%val, *, iostat = stat) val
  if (stat /= 0) val = defaultValue

end function getOptionInteger_

function getOptionLogical_(key, defaultValue) result(val)

  ! <<< Arguments >>>
  character(len = *), intent(in) :: key
  logical, intent(in) :: defaultValue

  ! <<< Result >>>
  logical :: val

  ! <<< Local variables >>>
  integer :: index

  call find(key, index)

  if (index == -1) then
     val = defaultValue
     return
  end if

  if (trim(dict(index)%val) == "true" .or. trim(dict(index)%val) == "TRUE" .or.              &
       trim(dict(index)%val) == "True") then
     val = .true.
  else if (trim(dict(index)%val) == "false" .or. trim(dict(index)%val) == "FALSE" .or.       &
       trim(dict(index)%val) == "False") then
     val = .false.
  else
     val = defaultValue
  end if

end function getOptionLogical_

function getOptionScalar_(key, defaultValue) result(val)

  ! <<< External modules >>>

  ! <<< Private members >>>

  ! <<< Arguments >>>
  character(len = *), intent(in) :: key
  real(mp), intent(in) :: defaultValue

  ! <<< Result >>>
  real(mp) :: val

  ! <<< Local variables >>>
  integer :: index, stat

  call find(key, index)

  if (index == -1) then
     val = defaultValue
     return
  end if

  read(dict(index)%val, *, iostat = stat) val
  if (stat /= 0) then
     val = defaultValue
  end if

end function getOptionScalar_

function getOptionString_(key, defaultValue) result(val)

  ! <<< Private members >>>

  ! <<< Arguments >>>
  character(len = *), intent(in) :: key
  character(len = *), intent(in) :: defaultValue

  ! <<< Result >>>
  character(len = STRING_LENGTH) :: val

  ! <<< Local variables >>>
  integer :: index

  call find(key, index)

  if (index == -1) then
     if (len_trim(defaultValue) == 0) then
        val = ""
     else
        read(defaultValue, '(A)') val
     end if
     return
  end if

  read(dict(index)%val, *) val

end function getOptionString_

subroutine getRequiredOptionInteger_(key, val, comm)

  ! <<< Arguments >>>
  character(len = *), intent(in) :: key
  integer, intent(out) :: val
  integer, intent(in), optional :: comm

  ! <<< Local variables >>>
  integer :: comm_, index, stat
  character(len = STRING_LENGTH) :: message

  if (present(comm)) comm_ = comm

  call find(key, index)
  if (index /= -1) read(dict(index)%val, *, iostat = stat) val

  if (index == -1 .or. stat /= 0) then
     write(message, "(3A)") "Required parameter '", trim(key), "' not found or invalid!"
     stop
  end if

end subroutine getRequiredOptionInteger_

subroutine getRequiredOptionLogical_(key, val, comm)

  ! <<< Arguments >>>
  character(len = *), intent(in) :: key
  logical, intent(out) :: val
  integer, intent(in), optional :: comm

  ! <<< Local variables >>>
  integer :: comm_, index, stat
  character(len = STRING_LENGTH) :: message

  if (present(comm)) comm_ = comm

  call find(key, index)
  if (index /= -1) then
     if (trim(dict(index)%val) == "true" .or. trim(dict(index)%val) == "TRUE" .or.           &
          trim(dict(index)%val) == "True") then
        val = .true.
     else if (trim(dict(index)%val) == "false" .or. trim(dict(index)%val) == "FALSE" .or.    &
          trim(dict(index)%val) == "False") then
        val = .false.
     else
        stat = -1
     end if
  end if

  if (index == -1 .or. stat /= 0) then
     write(message, "(3A)") "Required parameter '", trim(key), "' not found or invalid!"
     stop
  end if

end subroutine getRequiredOptionLogical_

subroutine getRequiredOptionScalar_(key, val, comm)

  ! <<< Arguments >>>
  character(len = *), intent(in) :: key
  real(mp), intent(out) :: val
  integer, intent(in), optional :: comm

  ! <<< Local variables >>>
  integer :: comm_, index, stat
  character(len = STRING_LENGTH) :: message

!  comm_ = MPI_COMM_WORLD
  if (present(comm)) comm_ = comm

  call find(key, index)
  if (index /= -1) then

     read(dict(index)%val, *, iostat = stat) val

  end if

  if (index == -1 .or. stat /= 0) then
     write(message, "(3A)") "Required parameter '", trim(key), "' not found or invalid!"
     stop
  end if

end subroutine getRequiredOptionScalar_

subroutine getRequiredOptionString_(key, val, comm)

  ! <<< Arguments >>>
  character(len = *), intent(in) :: key
  character(len = STRING_LENGTH), intent(out) :: val
  integer, intent(in), optional :: comm

  ! <<< Local variables >>>
  integer :: comm_, index
  character(len = STRING_LENGTH) :: message

  if (present(comm)) comm_ = comm

  call find(key, index)
  if (index /= -1) then
     read(dict(index)%val, *) val
  else
     write(message, "(3A)") "Required parameter '", trim(key), "' not found or invalid!"
     stop
  end if

end subroutine getRequiredOptionString_

    function determineMachinePrecision(filename) result(mp)

        !<< external modules >>
        use, intrinsic :: iso_fortran_env, only : output_unit

        !<< input variables >>
        character(len=STRING_LENGTH), intent(in) :: filename

        !<< output variables >>
        integer :: mp

        !<< internal variables >>
        character(len=STRING_LENGTH) :: MPOption

        call parseInputFile(filename)
        MPOption = GetOption('precision','64bit')
!        write(output_unit,'(A)') 'The machine precision is '//MPOption
!        flush(output_unit)
        
        select case (MPOption)
            case('32bit')
                mp = SELECTED_REAL_KIND(6,37)
            case('64bit')
                mp = SELECTED_REAL_KIND(15,307)
            case('128bit')
                mp = SELECTED_REAL_KIND(33,4931)
        end select
    end function

end module modInputHelper
