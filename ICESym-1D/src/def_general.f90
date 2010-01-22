module def_general

  use, intrinsic :: ISO_C_BINDING

contains

  subroutine open_unit(nu,namefile,nfsize) BIND(C)
    use, intrinsic :: ISO_C_BINDING
    integer(C_INT) :: nu,nfsize
    character(C_CHAR),dimension(nfsize) :: namefile
    integer :: io
    !character,dimension(:),pointer :: nfAux
    character(LEN=nfsize) :: nf
    do i=1,(nfsize+1)
       nf(i:i) = namefile(i)
    enddo
    write(*,*) namefile
    write(*,*) nf
    open(UNIT=nu,FILE=nf,ACTION='WRITE',STATUS='REPLACE',IOSTAT=io)
    if(io /= 0) then
      write(*,*) "El archivo "//nf//" no pudo ser abierto"
    endif
    return
  end subroutine open_unit

  subroutine close_unit(nu) BIND(C)
    use, intrinsic :: ISO_C_BINDING
    integer(C_INT) :: nu
    integer:: io
    close(nu)
    write(*,*) "El archivo se cerró"
    return
  end subroutine close_unit

end module def_general
