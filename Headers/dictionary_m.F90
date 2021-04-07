!> \file dictionary_m.f90
!! \brief Module file for dictionary_t

!> Dictionary type that uses strings for the keys and values
!!
!! Design:
!!  - djb2 hash function (D. J. Bernstein, see http://www.cse.yorku.ca/~oz/hash.html)
!!  - The strings are all "character(len=:), allocatable" variables
!!  - There is no linked list nor pointers, only allocatable arrays for the dynamic data structure
!!  - set rewrites existing entries without complaining
!!
!!%%%BMY
!!%%%BMY  Modified to return integer values (see !!%%%BMY comments for
!!%%%BMY  places where the code was changed) -- Bob Yantosca (06 Dec 2019)
!!%%%BMY

module dictionary_m
  implicit none

  private

  public :: dictionary_t

  !> Single entry in the dictionary
  type entry_t
     character(len=:), allocatable :: key
!!%%%BMY
!!%%%BMY Change value type from string to integer
!!%%%BMY     character(len=:), allocatable :: value
!!%%%BMY
     INTEGER :: value
  end type entry_t

  !> A bucket contains several entries
  type bucket_t
     type(entry_t), allocatable :: entries(:)
     integer :: current_size = 0
     integer :: current_idx = 0
   contains
     procedure :: find
  end type bucket_t

  !> The dictionary contains dict_size buckets (defined at run time)
  type dictionary_t
     type(bucket_t), allocatable :: buckets(:)
     integer :: dict_size = 0
   contains
     procedure :: djb2
     procedure :: set
     procedure :: get
     procedure :: init
     procedure :: show
     procedure :: destroy
  end type dictionary_t

  integer, parameter :: BUCKET_EMPTY = -2
  integer, parameter :: BUCKET_ENTRY_NOT_FOUND = -4

contains

  !> djb2 hash function
  !!
  !! \param this the dictionary_t object
  !! \param s a string
  !!
  !! \return the hash value between 0 and dict_size-1
  function djb2(this, s) result(r)
    class(dictionary_t), intent(in) :: this
    character(len=*), intent(in) :: s
    integer :: r

    integer :: i, l

    l = len(s)

    r = 5381

    do i = 1, l
       r = r*33 + ichar(s(i:i))
    end do

    r = modulo(r, this%dict_size)

  end function djb2

  !> Add or replace an entry in the dictionary
  !!
  !! \param this the dictionary_t object
  !! \param k the key
  !! \param v the value
  subroutine set(this, k, v)
    class(dictionary_t), intent(inout) :: this
    character(len=*), intent(in) :: k
!!%%%BMY
!!%%%BMY Change value type from string to integer
!!%%%BMY    character(len=*), intent(in) :: v
!!%%%BMY
    INTEGER, INTENT(IN) :: v

    type(bucket_t) :: tmp_bucket

    integer :: h, i, b_idx

    h = this%djb2(k) + 1

    b_idx = this%buckets(h)%find(k)

    if (b_idx == BUCKET_EMPTY) then
       ! allocate bucket for 1 entry
       ! also, means we can take the first entry
       allocate(this%buckets(h)%entries(1))
       this%buckets(h)%current_size = 1
       this%buckets(h)%current_idx = 1
       b_idx = 1
       this%buckets(h)%entries(1)%key = trim(k)
!!%%%BMY
!!%%%BMY Change value type from string to integer
!!%%%BMY       this%buckets(h)%entries(1)%value = trim(v)
!!%%%BMY
       this%buckets(h)%entries(1)%value = v
       ! the values are registered, exit
       return
    end if

    if (b_idx == BUCKET_ENTRY_NOT_FOUND) then
       ! copy and grow bucket entries

       allocate(tmp_bucket%entries(this%buckets(h)%current_size + 1))
       tmp_bucket%current_size = this%buckets(h)%current_size + 1
       tmp_bucket%current_idx = this%buckets(h)%current_idx + 1

       do i = 1, this%buckets(h)%current_size
          tmp_bucket%entries(i)%key = this%buckets(h)%entries(i)%key
          tmp_bucket%entries(i)%value = this%buckets(h)%entries(i)%value
       end do

       deallocate(this%buckets(h)%entries)
       allocate(this%buckets(h)%entries, source=tmp_bucket%entries)
       deallocate(tmp_bucket%entries)

       this%buckets(h)%current_size = tmp_bucket%current_size
       this%buckets(h)%current_idx = tmp_bucket%current_idx
       b_idx = this%buckets(h)%current_idx
    end if

    if (b_idx > 0) then
       this%buckets(h)%entries(b_idx)%key = k
       this%buckets(h)%entries(b_idx)%value = v
    end if

  end subroutine set

  !> Initialize a dictionary object
  !!
  !! \param this the dictionary_t object
  !! \param dict_size the size of the hash table
  subroutine init(this, dict_size)
    class(dictionary_t), intent(out) :: this
    integer, intent(in) :: dict_size

    allocate(this%buckets(dict_size))
    this%dict_size = dict_size

  end subroutine init

  !> Display the content of a dictionary
  !!
  !! \param this the dictionary_t object
  subroutine show(this)
    class(dictionary_t), intent(in) :: this

    integer :: i, j, s
    integer :: n

    n = 0
    do i = 1, this%dict_size
       s = this%buckets(i)%current_idx
       if (s > 0) then
             write(*,*) 'bucket   : ', i, ' size ', s
          do j = 1, s
             write(*,*) 'key      : ', this%buckets(i)%entries(j)%key
             write(*,*) 'value    : ', this%buckets(i)%entries(j)%value
          end do
       end if
    end do

  end subroutine show

  !> Find the "in-bucket" index for a given key
  !!
  !! Negative return values correspond to module-defined return codes.
  !!
  !! \param this the bucket_t object
  !! \param k the key
  !!
  !! \return the index (1-based) of the key in the bucket or a return code
  function find(this, k) result(r)
    class(bucket_t), intent(in) :: this
    character(len=*), intent(in) :: k
    integer :: r

    integer :: i

    if (this%current_size == 0) then
       r = BUCKET_EMPTY
       return
    end if

    r = BUCKET_ENTRY_NOT_FOUND
    do i = 1, this%current_size
       if (this%entries(i)%key == trim(k)) then
          r = i
          exit
       end if
    end do

  end function find

  !> Fetch an entry in the dictionary.
  !!
  !! \param this the dictionary_t object
  !! \param k the key
  !!
  !! \return the value if found, an empty string else
  function get(this, k) result(r)
    class(dictionary_t), intent(in) :: this
    character(len=*), intent(in) :: k
!!%%%BMY
!!%%%BMY Change result type from string to integer
!!%%%BMY    character(len=:), allocatable :: r
!!%%%BMY
    INTEGER :: r

    integer :: h, b_idx

    h = this%djb2(k) + 1

    b_idx = this%buckets(h)%find(k)

    if ( (b_idx == BUCKET_EMPTY) .or. &
         (b_idx == BUCKET_ENTRY_NOT_FOUND) ) then
!!%%%BMY
!!%%%BMY Set default return to -1, to match Ind_ function
!!%%%BMY       r = ''
!!%%%BMY
       r = -1
       return
    end if

    if (b_idx>0) then
       r = this%buckets(h)%entries(b_idx)%value
    end if

  end function get

!%%%BMY !> Free a dictionary object
!%%%BMY !!
!%%%BMY !! \param this the dictionary_t object
!%%%BMY !! \param dict_size the size of the hash table
  subroutine destroy(this)
    class(dictionary_t), intent(inout) :: this

    if ( allocated( this%buckets ) ) then
       deallocate( this%buckets )
    endif

  end subroutine destroy

end module dictionary_m
