!************************************************************
!
!  This example shows how to recursively traverse a file
!  using H5Ovisit.  The program prints all of
!  the objects in the file specified in FILE.  The default
!  file used by this example implements the structure described
!  in the User's Guide, chapter 4, figure 26.
!
!  This file is intended for use with HDF5 Library version 1.8
!  with --enable-fortran2003
!
!************************************************************
Module HDF5_LOFAR_Read
  USE HDF5
  USE ISO_C_BINDING
  use FitParams, only :  ImagingRun
  IMPLICIT NONE
  Integer, parameter :: Group_max=10
  Integer           :: Group_nr
  Character(len=50) ::  Group_names(Group_max)
  CHARACTER(LEN=140) :: filename
  Integer, parameter :: DSet_max=24
  Integer           :: DSet_nr
  Character(len=50) ::  DSet_names(DSet_max)
  INTEGER(HID_T), save :: file_ID ! File Handle
  INTEGER(HID_T), save :: group_id, dset_id ! handles
  INTEGER :: hdferr
  !
  Real*8 :: ANTENNA_POSITION(3)
  Integer :: DATA_LENGTH, SAMPLE_NUMBER_first, Absolute_TIME, Ant_ID, STATION_ID
  Real*8 :: DIPOLE_CALIBRATION_DELAY
  Integer, parameter :: List_dim=1500
  INTEGER(HID_T),save :: file_id_list(List_dim)       ! File identifier
  INTEGER(HID_T),save :: group_id_list(List_dim), dset_id_list(List_dim) ! handles
  Character(len=50),save ::  Group_name_list(0:List_dim)
  CHARACTER(LEN=140),save :: filename_list(0:List_dim)
  Character(len=50),save ::  DSet_list(0:List_dim)
  Integer,save           :: list_max=1
  !
CONTAINS
!------------------------
Subroutine GetFileName(filenr, nxx)
   use DataConstants, only : DataFolder
  IMPLICIT NONE
  Integer, intent(in) :: filenr
  Integer, intent(out) :: nxx
  !  Character(len=25), intent(in) :: DataFolder
  integer :: i
  !Character(len=100) :: Fdir
  !Character(len=100) :: Fname
  !  Fdir='/exp_app2/appexp1/public/raw_data/2018/'
  !  Fname='D20180813T153001.413Z/L664246_D20180813T153001.413Z_CS002_R000_tbb.h5'
  !  filename=trim(Fdir)//trim(Fname)
  Open(Unit=11,STATUS='old',ACTION='READ',FILE='directory.out',FORM ="formatted")
  Do i=1,filenr
    read(11,"(A)",IOSTAT=nxx) filename
    if(nxx.ne.0) return  ! stop('EOF reached in directory')
    !Write(2,*) 'file=',filename
  Enddo
  close(Unit=11)
  WRITE(2,'(A,i2,A,A)') 'Filenr ',filenr,': ',trim(filename)
End Subroutine GetFileName
!=========================================
Subroutine ListGroups
! based on MODULE g_visit from example "h5ex_g_visit_F03.f90" supplied by the HDF-group
! Open the file for future reading &
! obtain a list of Group_Names on this file
!
  IMPLICIT NONE
  !
  !INTEGER(HID_T) :: file_ID ! File Handle
  INTEGER :: status
  TYPE(C_FUNPTR) :: funptr
  TYPE(C_PTR) :: ptr
  INTEGER :: ret_value
  !CHARACTER(LEN=140) :: filename
  !
  write(2,*) '===================='
  !
  !
   CALL h5open_f(status)
   CALL H5Fopen_f(filename, H5F_ACC_RDONLY_F, file_ID, status)
  ! Begin iteration using H5Ovisit
  !
  WRITE(2,'(A,a)') "Objects infile:", trim(filename)
  !
  funptr = C_FUNLOC(op_func)
  CALL H5Ovisit_f(file_ID, H5_INDEX_NAME_F, H5_ITER_NATIVE_F, funptr, C_NULL_PTR, ret_value, status)
  !
  Return
End Subroutine ListGroups
!================================
Subroutine CloseFile
! Close the file
  IMPLICIT NONE
  INTEGER :: status
  CALL H5Fclose_f(file_ID, status)
End Subroutine CloseFile
!=====================================================
  INTEGER FUNCTION op_func(loc_id, name, info, cptr) bind(C)
!************************************************************
!
!  Operator function for H5Ovisit.  This function prints the
!  name and type of the object passed to it.
!
!************************************************************
    IMPLICIT NONE
    !
    INTEGER(HID_T), VALUE :: loc_id
    CHARACTER(LEN=1), DIMENSION(1:50) :: name ! We must have LEN=1 for bind(C) strings
                                              ! in order to be standard compliant
    TYPE(H5O_info_t) :: info
    CHARACTER(LEN=50) :: name_string = ' '
    TYPE(C_PTR) :: cptr
    INTEGER   :: i
    !
    name_string = ' '
    DO i = 1, 50
       IF(name(i)(1:1).EQ.C_NULL_CHAR) EXIT ! Read up to the C NULL termination
       name_string(i:i) = name(i)(1:1)
    EndDO
    !
    WRITE(2,"('/')",ADVANCE="NO")  !  Print root group in object path
    !
    ! Check if the current object is the root group, and if not print
    ! the full path name and type.
    !
    IF(name(1)(1:1) .EQ. '.')THEN         ! Root group, do not print '.'
        WRITE(2,"('  (Group)')")
        Group_nr=0
    ELSE
       IF(info%type.EQ.H5O_TYPE_GROUP_F)THEN
          WRITE(2,'(A,"  (Group)")') TRIM(name_string)
          Group_nr=Group_nr+1
          If(Group_nr.le.Group_max) then
            Group_names(Group_nr)=name_string
          Else
            write(2,"(A,i0)") 'More groups than expected, Group_nr=',Group_nr
          Endif
          !
       ELSE IF(info%type.EQ.H5O_TYPE_DATASET_F)THEN
          WRITE(2,'(A,"  (Dataset)")') TRIM(name_string)
       ELSE IF(info%type.EQ.H5O_TYPE_NAMED_DATATYPE_F)THEN
          WRITE(2,'(A,"  (Datatype)")') TRIM(name_string)
       ELSE
          WRITE(2,'(A,"  (Unknown)")') TRIM(name_string)
       EndIF
    EndIF
    !
    op_func = 0 ! return successful
  End FUNCTION op_func
!=======================================================
Subroutine ListGroupStructure(GroupName)
! This produces a list of DSet_names & opens the group handle
  IMPLICIT NONE
  CHARACTER(LEN=50) :: dsetname, groupname
  !INTEGER(HID_T) :: file_id       ! File identifier
  !INTEGER(HID_T) :: group_id, dset_id ! handles
  !INTEGER :: hdferr
  INTEGER :: storage_type ! Type of storage for links in group:
                          !   H5G_STORAGE_TYPE_COMPACT: Compact storage
                          !   H5G_STORAGE_TYPE_DENSE: Indexed storage
                          !   H5G_STORAGE_TYPE_SYMBOL_TABLE: Symbol tables
  INTEGER :: nlinks       ! Number of links in group
  INTEGER :: max_corder   ! Current maximum creation order value for group
  INTEGER(SIZE_T) :: size  ! Size of name
  INTEGER(HSIZE_T) :: i     ! Index
  CHARACTER(LEN=80) :: name ! Output buffer
  logical :: prn
!-----------
  INTEGER     ::   error ! Error flag
  !INTEGER     ::   j,k,m
  !
  !CALL h5open_f(error)
  !CALL h5fopen_f (trim(filename), H5F_ACC_RDONLY_F, file_id, error)
  CALL h5gopen_f(file_id, GroupName,group_id,error)
  !
  ! Get group info.
  CALL H5Gget_info_f(group_id, storage_type, nlinks, max_corder, hdferr)
  !
  ! Traverse links in the primary group using alphabetical indices
  ! (H5_INDEX_NAME).
  Write(2,'("Traversing group=",A)') trim(GroupName)
  DSet_nr=0
  DO i = 0, nlinks-1
    !
    ! Get name and size of name
    CALL H5Lget_name_by_idx_f(file_id, GroupName, H5_INDEX_NAME_F, H5_ITER_INC_F, i, name, hdferr, size)
    Write(2,'("Index ",i2,": ",A)') INT(i), TRIM(name)
    DSet_nr=DSet_nr+1
    If(DSet_nr.le.DSet_max) then
        DSet_names(DSet_nr)=name
    Else
        write(2,"(A,i0)") 'More Data_sets than expected, DSet_nr=',DSet_nr
    Endif
    !
  EndDO
  Call AttrRead(file_id, group_id,GroupName,prn)
  !
  !CALL h5gclose_f(group_id, error)
  !CALL h5fclose_f(file_id, error)
  !CALL h5close_f(error)
  return
End Subroutine ListGroupStructure
! ========================================
Subroutine CloseGroup
  IMPLICIT NONE
  INTEGER     ::   error ! Error flag
  CALL h5gclose_f(group_id, error)
  !CALL h5fclose_f(file_id, error)
  !CALL h5close_f(error)
End Subroutine CloseGroup
! ===============================================
Subroutine ListDataAtt(GroupName,DSetName, prnt)
  IMPLICIT NONE
  CHARACTER(LEN=50), intent(in) :: DSetName, groupname     ! Dataset name
  logical, intent(in), optional :: prnt
  !
  !INTEGER(HID_T) :: file_id       ! File identifier
  !INTEGER(HID_T) :: group_id, dset_id ! handles
  !INTEGER :: hdferr
  logical :: prn
!-----------
  !
  prn=.true.
  If(present(prnt)) then
    prn=prnt
  Endif
  !CALL h5open_f(error)
  !CALL h5fopen_f (trim(filename), H5F_ACC_RDONLY_F, file_id, error)
  !CALL h5gopen_f(file_id, GroupName,group_id,error)
  ! <<<<<<<<<<<<<<<<<<<<<<<<<<<
  CALL h5dopen_f (group_id, DSetName, dset_id, hdferr)
  Call AttrRead(group_id, dset_id,DSetName,prn)
    !Write(*,*) 'data_len=',data_len
  !CALL h5dclose_f(dset_id, error)
  !CALL h5gclose_f(group_id, error)
  !CALL h5fclose_f(file_id, error)
  !CALL h5close_f(error)
  !
  Return
End Subroutine ListDataAtt
! ========================================
Subroutine CloseDSet
  IMPLICIT NONE
  INTEGER     ::   error ! Error flag
  CALL h5dclose_f(dset_id, error)
  !CALL h5gclose_f(group_id, error)
  !CALL h5fclose_f(file_id, error)
  !CALL h5close_f(error)
End Subroutine CloseDSet
! ===============================================
! ===============================================
Subroutine GetData(Chunk, DSet_offset, DSet_dim)
  IMPLICIT NONE
  Integer, intent(in) :: DSet_offset, DSet_dim
  Integer*2, intent(out) :: Chunk(:)
  !INTEGER(HID_T) :: file_id       ! File identifier
  !INTEGER(HID_T) :: group_id, dset_id ! handles
  !
    call DataRead(dset_id, Chunk, DSet_offset, DSet_dim)
  !
  Return
End Subroutine GetData
! ========================================
Subroutine AttrRead(file_id,Obj_id,Obj_name,prn)
!
  IMPLICIT NONE
  ! This should map to REAL*8 on most modern processors
  INTEGER, PARAMETER :: real_kind_15 = SELECTED_REAL_KIND(15,307)
  INTEGER(HID_T), intent(in) :: Obj_id       ! File identifier
  CHARACTER(LEN=50), intent(in) :: Obj_name   ! Group or Dataset name
  Logical, intent(in) :: prn
  !Integer, intent(out), optional :: data_len
  Character(len=50) :: AttName
  INTEGER(HID_T) :: filetype_id, memtype_id    ! ????? identifier
  INTEGER(HID_T) :: attr_id,space_id ! handles
  Integer(HID_T) :: n, file_id
  INTEGER(HSIZE_T)  :: npoints
  INTEGER(HSIZE_T), DIMENSION(1:1)   :: dims
  INTEGER(HSIZE_T), DIMENSION(1:1)   :: maxdims
  !
  INTEGER, DIMENSION(:), ALLOCATABLE, TARGET :: rdata    ! Read buffer
  Real(KIND=real_kind_15), DIMENSION(:), ALLOCATABLE, TARGET :: xdata    ! Read buffer
  INTEGER(SIZE_T)  , PARAMETER :: sdim      = 9
  CHARACTER(LEN=sdim), DIMENSION(:), ALLOCATABLE, TARGET :: chdata
  TYPE(C_PTR), DIMENSION(:), ALLOCATABLE, TARGET :: vlCHdata ! Read buffer
  integer, parameter :: vclen=10
  CHARACTER(len = vclen, kind=c_char),  POINTER :: vcdata ! A pointer to a Fortran string
  TYPE(C_PTR) :: f_ptr
  !----------------
  !INTEGER :: hdferr
  INTEGER(SIZE_T) :: size  ! Size of name
  INTEGER(HSIZE_T) :: i     ! Index
  !-----------
  INTEGER :: F_class
  INTEGER     ::   error ! Error flag
  INTEGER     ::   j,k,m
  integer   :: attr_num
!
  write(2,*) ' =================================, get attributes for object=',trim(Obj_name)
  Call h5aget_num_attrs_f(Obj_id, attr_num, error)
  !write(2,'(3A,I4,i3)') 'Number of attributes for (object=',trim(Obj_name),')=',attr_num, error
  !
  !Initialize all to zero
  ANTENNA_POSITION=0.
  DATA_LENGTH=0; SAMPLE_NUMBER_first=0; Absolute_TIME=0; Ant_ID=0; STATION_ID=0
  DIPOLE_CALIBRATION_DELAY=0.
  !
  !obj_name=trim(SourceFile)
  DO i = 0, attr_num-1
    !Write(2,*) '==================, attribute #=',i
    Call h5aget_name_by_idx_f(file_id, Obj_name, H5_INDEX_NAME_F, H5_ITER_INC_F, i, AttName, error)
    !Write(2,'("h5aget_name_by_idx,file:",i2)') error
    !Write(2,'("==================, Attribute # & Name:",i2,": ",A)') i, trim(AttName)
    CALL h5aopen_f(Obj_id, AttName, attr_id, hdferr)
    !Write(2,'("h5aopen:",i2)') hdferr
    !
    ! Get the datatype of attributes and its dimensions.
    !
    CALL H5Aget_type_f(attr_id, filetype_id, hdferr)
    CALL h5tget_class_f(filetype_id, F_class, error)
    !
    If(F_class .eq. H5T_STRING_F) goto 9
    ! Get dataspace and allocate memory for read buffer.  This is a
    ! three dimensional attribute when the array datatype is included.
    !
    CALL H5Aget_space_f(attr_id, space_id, hdferr)
    !  Write(2,'("h5aget_space:",i2)') hdferr
    CALL H5Sget_simple_extent_ndims_f(space_id, k, hdferr)
    if(k.ne.1) then
        Write(2,*) 'H5Sget_simple_extent_ndims returns rank of the array=',k,', should be unity!!!'
        stop
    Endif
    !CALL H5Sget_simple_extent_npoints_f(space_id, npoints, hdferr)
    !  Write(2,*) 'H5Sget_simple_extent_npoints returns npoints=',npoints, 'total number of elements'
    CALL H5Sget_simple_extent_dims_f(space_id, dims, maxdims, hdferr)
    !Write(2,'("h5sget_simple_extent_dims:",i2)') hdferr
    !Write(2,*) 'idims:',dims
    !Write(2,*) 'ndims:',Maxdims
    !
    If(prn) Write(2,'(2a,i0,a)', ADVANCE='NO') trim(AttName),'(',dims(k),'): ['
    if(F_class .eq. H5T_FLOAT_F) then
      ALLOCATE(xdata(1:dims(k)) )
      CALL H5Aread_f(attr_id, filetype_id, xdata, dims, hdferr)
      DO j=1, dims(k)
        If(prn) Write(2,'(g14.5)', ADVANCE='NO') xdata(j)
      EndDO
      if(trim(AttName).eq."ANTENNA_POSITION_VALUE") ANTENNA_POSITION(1:3)=xdata(1:3)
      if(trim(AttName).eq."DIPOLE_CALIBRATION_DELAY_VALUE") DIPOLE_CALIBRATION_DELAY=xdata(1)
      DEALLOCATE(xdata)
    !
!  Real*8 :: ANTENNA_POSITION(3)
!  Integer :: DATA_LENGTH, SAMPLE_NUMBER, Absolute_TIME
!  Real*8 :: DIPOLE_CALIBRATION_DELAY_VALUE
    else if(F_class .eq. H5T_INTEGER_F) then
      ALLOCATE(rdata(1:dims(k)) )
      f_ptr = C_LOC(rdata)
      CALL H5Aread_f(attr_id, filetype_id, f_ptr, hdferr)
      DO j=1, dims(k)
        !Write(2,'(2a,i0,a,i13)') trim(AttName),'(',j,'):',rdata(j)
        If(prn) Write(2,'(i13)', ADVANCE='NO') rdata(j)
      EndDO
      !if(trim(AttName).eq."DATA_LENGTH" .and. present(data_len) ) data_len=rdata(1)
      if(trim(AttName).eq."DATA_LENGTH") DATA_LENGTH=rdata(1)
      if(trim(AttName).eq."SAMPLE_NUMBER") SAMPLE_NUMBER_first=rdata(1)
      if(trim(AttName).eq."TIME") Absolute_TIME=rdata(1)
      if(trim(AttName).eq."RCU_ID") Ant_ID=rdata(1)
      if(trim(AttName).eq."STATION_ID") STATION_ID=rdata(1)
      DEALLOCATE(rdata)
      !
    else if(F_class .eq. H5T_STRING_F) then
      ! same for variable length
      ALLOCATE(vlCHdata(1:dims(1)))
      f_ptr = C_LOC(vlCHdata(1))
      CALL h5Aread_f(attr_id, H5T_STRING, f_ptr, hdferr)
      DO m = 1, dims(1)
        CALL C_F_POINTER(vlCHdata(m), vcdata)
        j = 0
        DO
            IF(vcDATA(j+1:j+1).EQ.C_NULL_CHAR .OR. j.GE.vclen) EXIT
            j = j + 1
        EndDO
        !Write(2,'(A,i2,2A,"(",I0,"): ",A)') 'length=',j,', ',trim(AttName), m, vcdata(1:j)
        If(prn) Write(2,"(A,', ')", ADVANCE='NO') vcdata(1:j)
      End DO
      DEALLOCATE(vlCHdata)
      !
    else
        Write(2,*) "Stored datatype is of a different class"
    End if
    If(prn) Write(2,"(A)") ']'
    CALL H5Sclose_f(space_id, hdferr)
    !Write(*,'("h5sclose:",i2)') hdferr
    !
9   continue    !
    CALL H5Aclose_f(attr_id, hdferr)
    !Write(*,'("h5aclose:",i2)') hdferr
    CALL H5Tclose_f(filetype_id, hdferr)
    !Write(*,'("h5tclose:",i2)') hdferr
  Enddo  ! attr_num
End Subroutine AttrRead
! ===============================================
Subroutine DataRead(dset_id, Chunk, DSet_offset, DSet_dim)
!
  IMPLICIT NONE
  INTEGER(HID_T), intent(in) :: dset_id       !  identifier
  Integer, intent(in) :: DSet_offset, DSet_dim
  Integer*2, intent(out) :: Chunk(DSet_dim)
  !
  !INTEGER, PARAMETER :: dim     = 2048
  !INTEGER*2, DIMENSION(1:dim) :: rHyData    ! Read buffer, Hypercubed
  INTEGER(HSIZE_T), DIMENSION(1:1)   :: dims
  INTEGER(HSIZE_T), DIMENSION(1:1)   :: start, count, offset
  INTEGER :: closedferr
  !
  INTEGER(HID_T)  :: space_id, dcpl_id, LocalMem_id ! Handles
  !INTEGER*2, DIMENSION(1:S) :: rdata    ! Read buffer
  INTEGER :: i,j,k
  !
  ! set similar blocksize in local array:
  dims(1)=DSet_dim
  !write(*,*) 'Call h5screate_simple_f'
  Call h5screate_simple_f(1, dims, LocalMem_id, hdferr)
  If(hdferr.ne.0) write(2,*) 'error in call to h5screate_simple_f'
  start = 0  ! off-set (0 is minimum)
  count = DSet_dim     ! number of blocks
  !write(*,*) 'Call h5sselect_hyperslab_f'
  CALL h5sselect_hyperslab_f (LocalMem_id, H5S_SELECT_SET_F, start, count, hdferr)
  If(hdferr.ne.0) write(2,*) 'error in call to h5sselect_hyperslab_f'
  !
  ! Define and select the hyperslab to use for reading part of the data.
  !
  CALL h5dget_space_f(dset_id, space_id, hdferr)
  If(hdferr.ne.0) write(2,*) 'error in call to h5dget_space_f'
  offset = DSet_offset  ! off-set (0 is minimum)
  count = DSet_dim     ! number of blocks
  CALL h5sselect_hyperslab_f (space_id, H5S_SELECT_SET_F, offset, count, hdferr)
  If(hdferr.ne.0) write(2,*) 'error in call to h5sselect_hyperslab_f',DSet_offset, DSet_dim
  !
  ! Read the data using the previously defined hyperslabs.
  !
  CALL h5dread_f(dset_id, H5T_STD_I16LE, Chunk, dims, hdferr, file_space_id=space_id, mem_space_id=LocalMem_id)
  !
  !Write(2, "('Data as read from disk with offset=',i0)") offset
  CALL h5sclose_f(space_id, closedferr)
  !Write(2,*) 'close space_id error=',hdferr
  CALL h5sclose_f(LocalMem_id, closedferr)
  return
End Subroutine DataRead
! ===============================================
Subroutine GetDataChunk(GroupName,DSetName, Chunk, DSet_offset, DSet_dim, prnt, DataReadErr)
  IMPLICIT NONE
  CHARACTER(LEN=50), intent(in) :: DSetName, groupname     ! Dataset name
  Integer, intent(in) :: DSet_offset, DSet_dim
  Integer*2, intent(out) :: Chunk(:)
  logical, intent(in), optional :: prnt
  integer, intent(out), optional :: DataReadErr
  !
  !INTEGER :: hdferr
  logical :: Dset_get, prn, CloseH5=.false.
!-----------
  INTEGER     ::   error=0 ! Error flag
  INTEGER     ::   List_nr, j,k,m
  !
  prn=.true.
  If(present(prnt)) then
    prn=prnt
  Endif
   !
   If(list_max.eq.1) filename_list='.'
   !write(2,*) '-----',list_max,trim(filename),' ; ',trim(GroupName),' ; ',trim(DSetName)
   error=0
   Do List_nr=1,list_max
      If((filename_list(List_nr).eq.trim(filename)) .and. (Group_name_list(List_nr).eq.trim(GroupName)) .and. &
            (DSet_list(List_nr).eq.trim(DSetName))) then
         exit
      Elseif(filename_list(List_nr).eq.'.') then  ! add a new entry to the list
         filename_list(List_nr) = trim(filename)
         Group_name_list(List_nr) = trim(GroupName)
         DSet_list(List_nr) = trim(DSetName)
         !write(2,*) 'GetDataChunk:List_nr=',List_nr,trim(filename),' ; ',trim(GroupName),' ; ',trim(DSetName)
         If(list_max.eq.1) CALL h5open_f(error)
         If(.not.(filename_list(List_nr-1).eq.trim(filename))) then ! a new file is opened
            !write(2,*) 'file'
            CALL h5fopen_f (trim(filename), H5F_ACC_RDONLY_F, file_id_list(List_nr), error)
            CALL h5gopen_f(file_id_list(List_nr), GroupName,group_id_list(List_nr),error)
         ElseIf(.not.(Group_name_list(List_nr-1).eq.trim(GroupName))) then ! a new group on an existing file
            !write(2,*) 'group'
            file_id_list(List_nr)=file_id_list(List_nr-1)
            CALL h5gopen_f(file_id_list(List_nr), GroupName,group_id_list(List_nr),error)
         Else ! only a new dataset is opened
            !write(2,*) 'dset'
            file_id_list(List_nr)=file_id_list(List_nr-1)
            group_id_list(List_nr)=group_id_list(List_nr-1)
         Endif
         If(error.ne.0) then
            write(2,*) '****** Error in GetDataChunk: ',trim(filename),' ; ',trim(GroupName), ' ; ',trim(DSetName)
            write(2,*) 'Non-existing file?'
            write(*,*) 'Non-existing file????'
            Stop 'GetDataChunk: file-open problem'
         endif
         CALL h5dopen_f (group_id_list(List_nr), DSetName, dset_id_list(List_nr), hdferr)
         list_max=list_max+1
         If(list_max.gt.List_dim) then
            Write(2,*) 'ERROR: list_max .gt. List_dim:',list_max,' .gt.',List_dim
            Stop 'GetDataChunk: data-list exceeded'
         Endif
         Exit
      Endif
   Enddo
   If(hdferr.ne.0) then
      write(2,*) '****** Error in GetDataChunk, h5dopen_f: ', trim(DSetName)
      write(2,*) 'Probably due to a broken SSHFS connection to the data repository; retry or restore'
      write(*,*) 'Lost SSHFS connection to data repository???'
      Stop 'GetDataChunk: data-open problem'
   Endif
   !write(2,*) 'Call DataRead'
   !write(*,*) 'Call DataRead'
   DataReadErr=0
    Call DataRead(dset_id_list(List_nr), Chunk, DSet_offset, DSet_dim)
    If(hdferr.ne.0) then ! set in call to 'h5dread_f' in 'GetDataChunk'
        If(.not. ImagingRun) write(*,*) '!!HDF5 error captured for antenna ',trim(DSetName),', no problem!!!'
        If(.not. ImagingRun) write(2,*) 'error in call to h5dread_f, data for antenna ',trim(DSetName),' are zeroed'
        Chunk=0
        DataReadErr=-1
    endif

   If(CloseH5) then
      Do List_nr=1,list_max-1
         CALL h5dclose_f(dset_id_list(List_nr), error)
         CALL h5gclose_f(group_id_list(List_nr), error)
         CALL h5fclose_f(file_id_list(List_nr), error)
      Enddo
      CALL h5close_f(error)
  Endif
  !
  Return
End Subroutine GetDataChunk
! ===============================================
Subroutine CloseDataFiles()
   IMPLICIT NONE
   INTEGER     ::   error=0 ! Error flag
   INTEGER     ::   List_nr
   !
   List_nr=list_max-1
   Do List_nr=list_max-1,2,-1
      CALL h5dclose_f(dset_id_list(List_nr), error)
      !write(2,*) 'closing file handle#',List_nr,trim(filename_list(List_nr)), &
      !   '; ',trim(Group_name_list(List_nr)),'; ',trim(DSet_list(List_nr))
      If(filename_list(List_nr-1).ne.filename_list(List_nr)) then ! close group and file
         CALL h5gclose_f(group_id_list(List_nr), error)
         CALL h5fclose_f(file_id_list(List_nr), error)
         !write(2,*) 'closing file and group',error
         !write(2,*) List_nr,trim(filename_list(List_nr))
         !write(2,*) List_nr-1,trim(filename_list(List_nr-1))
      Else If(Group_name_list(List_nr-1).ne.Group_name_list(List_nr)) then ! close group
         CALL h5gclose_f(group_id_list(List_nr), error)
         !write(2,*) 'closing group',error,trim(Group_name_list(List_nr))
      EndIF
      !CALL h5close_f(error)
      If(error.ne.0) then
         write(2,*) 'H5 closing error code=',error
         exit
      EndIf
   Enddo
   list_max=1
   List_nr=list_max
   CALL h5dclose_f(dset_id_list(List_nr), error)
   CALL h5gclose_f(group_id_list(List_nr), error)
   CALL h5fclose_f(file_id_list(List_nr), error)
   CALL h5close_f(error)
   flush(unit=2)
   Return
End Subroutine CloseDataFiles
!==========================================
End Module HDF5_LOFAR_Read
