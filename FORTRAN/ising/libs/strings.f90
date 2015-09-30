module strings

  use globales,     only: dp
  
  implicit none

contains

!===============================================================================
! STR_TO_REAL converts an arbitrary string to a real(8)
!===============================================================================

  function str_to_real(string) result(num)

    character(*), intent(in) :: string
    real(8)                  :: num

    integer :: ioError

    ! Read string
    read(UNIT=string, FMT=*, IOSTAT=ioError) num
    if (ioError > 0) num = huge(0) 

  end function str_to_real

!===============================================================================
! STARTS_WITH determines whether a string starts with a certain
! sequence of characters
!===============================================================================

  logical function starts_with(str, seq)

    character(*) :: str ! string to check
    character(*) :: seq ! sequence of characters

    integer :: i
    integer :: i_start
    integer :: str_len
    integer :: seq_len

    str_len = len_trim(str)
    seq_len = len_trim(seq)

    ! determine how many spaces are at beginning of string
    i_start = 0
    do i = 1, str_len
      if (str(i:i) == ' ' .or. str(i:i) == achar(9)) cycle
      i_start = i
      exit
    end do

    ! Check if string starts with sequence using INDEX intrinsic
    if (index(str(1:str_len), seq(1:seq_len)) == i_start) then
      starts_with = .true.
    else
      starts_with = .false.
    end if

  end function starts_with

!===============================================================================
! STR_TO_INT converts a string to an integer.
!===============================================================================

  function str_to_int(str) result(num)

    character(*), intent(in) :: str
    integer(8) :: num

    character(5) :: fmt
    integer      :: w
    integer      :: ioError

    ! Determine width of string
    w = len_trim(str)

    ! Create format specifier for reading string
    write(UNIT=fmt, FMT='("(I",I2,")")') w

    ! read string into integer
    read(UNIT=str, FMT=fmt, IOSTAT=ioError) num
    if (ioError > 0) num = huge(0) 

  end function str_to_int


end module strings
