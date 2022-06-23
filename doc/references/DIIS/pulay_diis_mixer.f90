!****h* FHI-aims/pulay_diis_mixer
!  NAME
!    pulay_diis_mixer - provides the following mixing routines:
!                 
!                 o Pulay/DIIS mixer
!
!  SYNOPSIS

module pulay_diis_mixer

!  PURPOSE
!    Provides mixing routines based on the DIIS mixing scheme:
!      P. Pulay, Convergence acceleration of iterative sequences. the case of
!      scf iteration,
!      Chemical Physics Letters, Volume 73, Issue 2, 1980, Pages 393-398,
!      ISSN 0009-2614, http://dx.doi.org/10.1016/0009-2614(80)80396-4.
!
!    The following routines are intended to be universal. They are designed
!    to have no nasty side effects as the only dependencies are
!    LAPACK/BLAS and everything is stored in an instance of the corresponding 
!    mixer type owned by the caller.
!
!  USAGE
!    (illustrated at the example of the double precision vector mixer)
!    A. First, an instance of a mixer type has to be created:
!       >>> type(DIISdpVectorMixer) :: mixer
!       >>> mixer = diis_vector_mixer_init(...)
!
!    B1. Then, in every (e.g. SCF) step, a new vector can be added to the mixer
!       >>> diis_vector_mixer_add(...)
!
!    B2. In order to calculated the mixed vector, simply call
!       >>> diis_vector_mixer_solve(...)
!       Attention: Input vector is overwritten!
!
!    Done. :)
!
!    B1 and B2 can also be done together with a single call to
!       >>> diis_vector_mixer_step(...)
!
!    In order to use the mixer in case of distributed vectors, one has to
!    supply a proper function "vdot_func" to the "diis_*_add" or "diis_*_step"
!    routine which calculates the dot-product of two such distributed vectors
!    (which involves some sort of communication). As an example, see the
!    definition and usage of "dpvector_dotfunc_local" here.
!
!    A note on the employed Fortran standard:
!      I would have liked to store the vector dot-product function in the 
!      derived type of the mixer as a type-bound procedure 
!      once during initialization. 
!      However, this is not possible within the F95 standard and the usage of 
!      Fortran2003 standard is discouraged in the project at the time of 
!      writing this piece of code.
!
!  USES
   use types, only: dp
   implicit none

!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  SEE ALSO
!    Volker Blum, Ralf Gehrke, Felix Hanke, Paula Havu, Ville Havu,
!    Xinguo Ren, Karsten Reuter, and Matthias Scheffler,
!    "Ab initio simulations with Numeric Atom-Centered Orbitals: FHI-aims",
!    Computer Physics Communications 180 (2009), 2175-2196.
!  COPYRIGHT
!    Max-Planck-Gesellschaft zur Foerderung der Wissenschaften
!    e.V. Please note that any use of the "FHI-aims-Software" is subject to
!    the terms and conditions of the respective license agreement."
!  HISTORY
!    Development version, FHI-aims (2016).
!  SOURCE

   private  ! make everything private by default

   ! PUBLIC ROUTINES

   public :: diis_vector_mixer_init
   public :: diis_vector_mixer_step
   public :: diis_vector_mixer_add
   public :: diis_vector_mixer_solve

   ! Generic interfaces

   interface diis_vector_mixer_init
      module procedure diis_dpvector_mixer_init
   end interface

   interface diis_vector_mixer_step
      module procedure diis_dpvector_mixer_step
      module procedure diis_dpvector_mixer_step_callback
   end interface

   interface diis_vector_mixer_add
      module procedure diis_dpvector_mixer_add
      module procedure diis_dpvector_mixer_add_callback
   end interface

   interface diis_vector_mixer_solve
      module procedure diis_dpvector_mixer_solve
   end interface

   ! Return codes

   type :: DIISReturnCodes
      integer :: SUCCESS=0
      integer :: UNINIT=-100
      integer :: LENGTHMISMATCH=-102
   end type
   type(DIISReturnCodes), parameter, public :: &
      diis_return_codes=DIISReturnCodes()

   ! Types

   type, public :: DIISdpVectorMixer
      real(dp) :: damping_factor=0
      ! private
      integer, private :: n_store_max=-1
      integer, private :: i_cont=-1
      integer, private :: vec_length=-1
      real(dp), allocatable, private :: storage_vec(:,:)
      real(dp), allocatable, private :: error_matrix(:)
      real(dp), allocatable, private :: weights(:)
   end type


contains


!-------------------------------------------------------------------------------
!****s* pulay_diis_mixer/dpvector_dotfunc_local
!  NAME
!    dpvector_dotfunc_local
!  SYNOPSIS

function dpvector_dotfunc_local(v1, v2) result(d)

!  PURPOSE
!    This is an example function to calculate the dot-product of two vectors
!    locally, i.e. without communication (e.g. to other MPI threads)
!
!  USES
   implicit none
!  ARGUMENTS
   real(dp), intent(in) :: v1(:), v2(:)
   real(dp) :: d
!  INPUTS
!   o v1, v2 -- vectors of same length N
!  RETURNS
!   o d -- local dot product: \sum_{i=1}^{N} ( v1(i)*v2(i) )
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Development version, FHI-aims (2016).
!  SOURCE

   real(dp), external :: DDOT

   integer :: N

   ! missing assertion: size(v1) == size(v2)
   N = size(v1)
   d = DDOT(N, v1, 1, v2, 1)
end function dpvector_dotfunc_local


!-------------------------------------------------------------------------------
!****s* pulay_diis_mixer/diis_dpvector_mixer_init
!  NAME
!    diis_dpvector_mixer_init
!  SYNOPSIS

subroutine diis_dpvector_mixer_init(mixer, n_max_mixing, length, damping_factor)

!  PURPOSE
!    This subroutine initializes an instance of type DIISdpVectorMixer
!
!  USES
   implicit none
!  ARGUMENTS
   type(DIISdpVectorMixer), intent(out) :: mixer
   integer, intent(in) :: n_max_mixing, length
   real(dp), intent(in), optional :: damping_factor
!  INPUTS
!   o n_max_mixing -- max number of vectors to be mixed
!   o length -- length of vectors for storage (i.e. the part of interest)
!   o damping_factor -- Valid range is 0..1.
!     Be df the damping factor, then the return vector v_r
!     of the mixer when adding a new vector v_new will be
!        v_r = df * v_new + (1-df) * v_mix
!     where v_mix is the mixed vector.
!     A value of 1 means no mixing but only new vector, whereas a value 
!     of 0 means full mixing.
!  OUTPUTS
!   o mixer -- initialized instance of the corresponding mixer type
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Development version, FHI-aims (2016).
!  SOURCE

   mixer % n_store_max = n_max_mixing
   mixer % vec_length = length
   mixer % i_cont = -1

   allocate(mixer % storage_vec(mixer % vec_length, mixer % n_store_max))
   allocate(mixer % error_matrix( &
      (mixer % n_store_max*(mixer % n_store_max+1))/2) &
   )
   allocate(mixer % weights(mixer % n_store_max))

   if (present(damping_factor)) &
      mixer % damping_factor = min(1.0_dp,max(0.0_dp,damping_factor))

end subroutine diis_dpvector_mixer_init



!-------------------------------------------------------------------------------
!****s* pulay_diis_mixer/diis_dpvector_mixer_step_callback
!  NAME
!    diis_dpvector_mixer_step_callback
!  SYNOPSIS

subroutine diis_dpvector_mixer_step_callback(mixer, length, new_vector, vdot_func, &
   returncode)

!  PURPOSE
!    This procedure accepts a new vector, calculates the proper mixing
!    coefficients and overwrites the new vector with the mixed one.
!
!  USES
   implicit none
!  ARGUMENTS
   type(DIISdpVectorMixer), intent(inout) :: mixer
   integer, intent(in) :: length
   real(dp), intent(inout) :: new_vector(length)
   interface
      real(dp) function vdot_func(v1, v2)
         use types, only: dp
         real(dp), intent(in) :: v1(:), v2(:)
      end function
   end interface
   integer, intent(out) :: returncode
!  INPUTS
!   o mixer -- instance of the corresponding mixer type
!   o length -- length of the vector (i.e. the part of interest)
!   o new_vector -- on entry, vector to be added to the mixer
!  OUTPUT
!   o mixer -- updated instance of the corresponding mixer type
!   o new_vector -- on exit, overwritten by mixed vector
!   o returncode -- see DIISReturnCodes
!           SUCCESS: everything fine
!           UNINIT: mixer not initialized
!           LENGTHMISMATCH: length of vector does not match
!           else: error code returned by solver routine
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Development version, FHI-aims (2016).
!  SOURCE

   returncode = diis_return_codes % SUCCESS
   if (mixer % n_store_max .le. 0) then
      returncode = diis_return_codes % UNINIT
   elseif (mixer % vec_length .ne. length) then
      returncode = diis_return_codes % LENGTHMISMATCH
   endif

   if (returncode .eq. diis_return_codes % SUCCESS) then
      if (mixer % i_cont .lt. 0) then
         ! this is the initial vector. no differences yet
         mixer % storage_vec(:,mixer % n_store_max) = new_vector
         mixer % i_cont = 0
      else
         call diis_dpvector_mixer_add_callback( &
            mixer, length, new_vector, vdot_func)
         call diis_dpvector_mixer_solve(mixer, length, new_vector, returncode)
      endif
   endif

end subroutine diis_dpvector_mixer_step_callback


!-------------------------------------------------------------------------------
!****s* pulay_diis_mixer/diis_dpvector_mixer_step
!  NAME
!    diis_dpvector_mixer_step
!  SYNOPSIS

subroutine diis_dpvector_mixer_step(mixer, length, new_vector, returncode)

!  PURPOSE
!    This procedure accepts a new vector, calculates the proper mixing
!    coefficients and overwrites the new vector with the mixed one.
!    This function assumes that the input vector is not distributed.
!    Thus, the local vector dot-product "dpvector_dotfunc_local" is used and 
!    no communication is performed.
!    For the meaning of the arguments, see "diis_dpvector_mixer_step_callback"

   implicit none
   type(DIISdpVectorMixer), intent(inout) :: mixer
   integer, intent(in) :: length
   real(dp), intent(inout) :: new_vector(length)
   integer, intent(out) :: returncode

   call diis_dpvector_mixer_step_callback( &
      mixer=mixer, &
      length=length, &
      new_vector=new_vector, &
      vdot_func=dpvector_dotfunc_local, &
      returncode=returncode &
   )
end subroutine diis_dpvector_mixer_step



!-------------------------------------------------------------------------------
!****s* pulay_diis_mixer/diis_dpvector_mixer_add_callback
!  NAME
!    diis_dpvector_mixer_add_callback
!  SYNOPSIS

subroutine diis_dpvector_mixer_add_callback(mixer, length, new_vector, vdot_func)

!  PURPOSE
!    This procedure updates the mixer instance with the given vector.
!
!  USES
   implicit none
!  ARGUMENTS
   type(DIISdpVectorMixer), intent(inout) :: mixer
   integer, intent(in) :: length
   real(dp), intent(in) :: new_vector(length)
   interface
      real(dp) function vdot_func(v1, v2)
         use types, only: dp
         real(dp), intent(in) :: v1(:), v2(:)
      end function
   end interface
!  INPUTS
!   o mixer -- instance of the corresponding mixer type
!   o length -- length of the vector (i.e. the part of interest)
!   o new_vector -- vector to be added to the mixer
!   o vdot_func -- function that performs dot product operation on
!                  possibly distributed data structure of new_vector
!  OUTPUT
!   o mixer -- updated instance of the corresponding mixer type
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Development version, FHI-aims (2016).
!  SOURCE

   integer :: i_new, i_old, i_sym, i_last, max_store
   real(dp), allocatable :: diff_new(:), diff_var(:)
   real(dp) :: weight_new

   allocate(diff_new(length))
   allocate(diff_var(length))

   mixer % i_cont = mixer % i_cont + 1
   max_store = min(mixer % i_cont, mixer % n_store_max)

   ! ATTENTION:
   ! overwrite old vector at the end! (needed for otf differences)

   ! this is the index where the new entry goes in the storage
   i_new = modulo( mixer % i_cont-1, mixer % n_store_max ) + 1
   ! this is the last vector added
   i_last = modulo( mixer % i_cont-2, mixer % n_store_max ) + 1
   diff_new = new_vector - mixer % storage_vec(:,i_last)
   weight_new = get_weight(diff_new)
   diff_new = diff_new * weight_new
   ! update Pulay matrix
   !  (indices here are formulated in an upper symmetric manner)
   ! get index of (1,i_new) in Pulay matrix
   i_sym = (i_new*(i_new-1))/2 + 1
   ! walk down to (i_new,i_new) (column)
   i_last = mixer % n_store_max
   do i_old = 1, i_new-1, 1
      diff_var = mixer % weights(i_old) * &
            (mixer % storage_vec(:,i_old) - mixer % storage_vec(:,i_last))
      mixer % error_matrix(i_sym) = vdot_func(diff_new, diff_var)
      i_last = i_old
      i_sym = i_sym + 1
   enddo ! walk column
   ! Place diagonal element at (i_new,i_new)
   mixer % error_matrix(i_sym) = vdot_func(diff_new, diff_new)
   ! walk from i_new to the right (row)
   do i_old = i_new + 1, max_store ! limit to occupied size
      i_sym = i_sym + i_old - 1
      diff_var = mixer % weights(i_old) * &
            (mixer % storage_vec(:,i_old) - mixer % storage_vec(:,i_old-1))
      mixer % error_matrix(i_sym) = vdot_func(diff_new, diff_var)
   enddo ! walk row
   ! now put vector in storage
   mixer % storage_vec(:,i_new) = new_vector
   mixer % weights(i_new) = weight_new

   deallocate(diff_var)
   deallocate(diff_new)

   contains

   function vector_norm(v) result(norm)
      real(dp), intent(in) :: v(:)
      real(dp) :: norm

      norm = sqrt(vdot_func(v,v))
   end function vector_norm

   function get_weight(v) result(w)
      real(dp), intent(in) :: v(:)
      real(dp) :: w

      real(dp), parameter :: weight_threshold=1.e-18_dp

      w = vector_norm(v)
      if (w.lt.weight_threshold) then
         w = 1
      else
         w = 1/w
      endif
   end function get_weight

end subroutine diis_dpvector_mixer_add_callback


!-------------------------------------------------------------------------------
!****s* pulay_diis_mixer/diis_dpvector_mixer_add
!  NAME
!    diis_dpvector_mixer_add
!  SYNOPSIS

subroutine diis_dpvector_mixer_add(mixer, length, new_vector)

!  PURPOSE
!    This procedure updates the mixer instance with the given vector.
!    This function assumes that the input vector is not distributed.
!    Thus, the local vector dot-product "dpvector_dotfunc_local" is used and 
!    no communication is performed.
!    For the meaning of the arguments, see "diis_dpvector_mixer_add_callback"

   implicit none
   type(DIISdpVectorMixer), intent(inout) :: mixer
   integer, intent(in) :: length
   real(dp), intent(in) :: new_vector(length)

   call diis_dpvector_mixer_add_callback( &
      mixer=mixer, &
      length=length, &
      new_vector=new_vector, &
      vdot_func=dpvector_dotfunc_local &
   )
end subroutine diis_dpvector_mixer_add



!-------------------------------------------------------------------------------
!****s* pulay_diis_mixer/diis_dpvector_mixer_solve
!  NAME
!    diis_dpvector_mixer_solve
!  SYNOPSIS

subroutine diis_dpvector_mixer_solve(mixer, length, &
               mixed_vector, returncode)

!  PURPOSE
!    This procedure returns the mixed value calculated from the mixer's storage.
!
!  USES
   implicit none
!  ARGUMENTS
   type(DIISdpVectorMixer), intent(in) :: mixer
   integer, intent(in) :: length
   real(dp), intent(out) :: mixed_vector(length)
   integer, intent(out) :: returncode
!  INPUTS
!   o mixer -- instance of the corresponding mixer type
!   o length -- length of the vector (i.e. the part of interest)
!  OUTPUT
!   o mixed_vector -- self-explanatory
!   o returncode -- LAPACK return value
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Development version, FHI-aims (2016).
!  SOURCE

   logical, parameter :: complete_orthogonal_factorization=.true.
   real(dp), parameter :: sum_of_coefficients_threshold=1.e-16_dp
   real(dp), parameter :: conditioning_threshold=1.e-16_dp

   character(64) :: fmtstr

   integer :: max_store, N, i_col, i_row, i_sym, lwork, rank
   real(dp) :: work_dummy(1)
   integer, allocatable :: ipiv(:)
   real(dp), allocatable :: A(:,:), b(:), work(:)
   real(dp) :: sum_of_coefficients, mixing_factor

   max_store = min(mixer%i_cont, mixer%n_store_max)
   N = max_store + 1 ! Laplacian

   allocate(A(N,N))
   allocate(ipiv(N))
   allocate(b(N))

   ipiv = 0

   b(1:max_store) = 0
   b(N) = -1

   ! copy upper triangular part from error matrix
   i_sym = 0
   do i_col = 1, max_store
      do i_row = 1, i_col
         i_sym = i_sym + 1
         A(i_row,i_col) = mixer % error_matrix(i_sym)
      enddo ! i_row
   enddo ! i_col
   ! fill last column
   do i_row = 1, max_store
      A(i_row,N) = -mixer % weights(i_row)
   enddo
   A(N,N) = 0

   ! mirror matrix along diagonal if necessary
   if (complete_orthogonal_factorization) then
      do i_col = 1, max_store
         do i_row = i_col+1, N
            A(i_row,i_col) = A(i_col,i_row)
         enddo ! i_row
      enddo ! i_col
   endif

   lwork = -1
   if (complete_orthogonal_factorization) then
      call DGELSY(N, N, 1, A, N, b, N, ipiv, conditioning_threshold, rank, &
         work_dummy, lwork, returncode)
   else
      call DSYSV('Upper', N, 1, A, N, ipiv, b, N, work_dummy, lwork, returncode)
   endif

   lwork = int(work_dummy(1))
   allocate(work(lwork))
   if (complete_orthogonal_factorization) then
      call DGELSY(N, N, 1, A, N, b, N, ipiv, conditioning_threshold, rank, &
         work, lwork, returncode)
   else
      call DSYSV('Upper', N, 1, A, N, ipiv, b, N, work, lwork, returncode)
      rank = N
   endif
   deallocate(work)
   ! now, coefficients := b(1:max_store)

   sum_of_coefficients = 0
   do i_row = 1, max_store
      ! re-weight coefficients
      b(i_row) = b(i_row) * mixer % weights(i_row)
      ! sum up coeffs
      sum_of_coefficients = sum_of_coefficients + b(i_row)
   enddo
   ! check sum of coefficients, set mixing factor
   if (abs(1-sum_of_coefficients) .gt. sum_of_coefficients_threshold) then
      ! return only current vector
      mixing_factor = 0
   else
      mixing_factor = 1 - mixer % damping_factor
   endif

   ! apply mixing factor
   do i_row = 1, max_store
      b(i_row) = b(i_row) * mixing_factor
   enddo
   ! include damping factor in coefficient
   b(mixer % i_cont) = b(mixer % i_cont) + 1 - mixing_factor

   call DGEMV('Normal', mixer % vec_length, max_store, &
               1, mixer % storage_vec, mixer % vec_length, &
               b, 1, &
               0, mixed_vector, 1)

   deallocate(b)
   deallocate(ipiv)
   deallocate(A)

end subroutine diis_dpvector_mixer_solve



end module pulay_diis_mixer

