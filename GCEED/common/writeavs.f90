!
!  Copyright 2017 SALMON developers
!
!  Licensed under the Apache License, Version 2.0 (the "License");
!  you may not use this file except in compliance with the License.
!  You may obtain a copy of the License at
!
!      http://www.apache.org/licenses/LICENSE-2.0
!
!  Unless required by applicable law or agreed to in writing, software
!  distributed under the License is distributed on an "AS IS" BASIS,
!  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
!  See the License for the specific language governing permissions and
!  limitations under the License.
!
subroutine writeavs(fp,suffix,tmatbox_l)
  use scf_data
  implicit none
  integer, intent(IN) :: fp
  real(8), intent(IN) :: tmatbox_l(lg_sta(1):lg_end(1),lg_sta(2):lg_end(2),  &
                                   lg_sta(3):lg_end(3))
  character(30),intent(in):: suffix
  character(30):: filename
  integer::ix,iy,iz
  integer::jj
  integer::jsta,jend
  character(8)  :: filenumber_data
  
  if(numfile_movie>=2)then
    if(myrank<numfile_movie)then
      write(filenumber_data, '(i8)') myrank
      filename = trim(suffix)//"."//adjustl(filenumber_data)
      open(fp,file=filename)
      jsta=myrank*(lg_num(1)*lg_num(2)*lg_num(3))/numfile_movie+1
      jend=(myrank+1)*(lg_num(1)*lg_num(2)*lg_num(3))/numfile_movie
      do jj=jsta,jend
        if(abs(tmatbox_l(icoo1d(1,jj),icoo1d(2,jj),icoo1d(3,jj)))>=1.0d-10)then
          write(fp,'(e20.8)') tmatbox_l(icoo1d(1,jj),icoo1d(2,jj),icoo1d(3,jj))
        else
          write(fp,'(a1)') "0"
        end if
      enddo
      close(1)
    end if
  else
    if(myrank==0)then
      filename=trim(suffix)
      open(fp,file=filename)
      do iz=lg_sta(3),lg_end(3)
      do iy=lg_sta(2),lg_end(2)
      do ix=lg_sta(1),lg_end(1)
        if(abs(tmatbox_l(ix,iy,iz))>=1.0d-10)then
          write(fp,'(e20.8)') tmatbox_l(ix,iy,iz)
        else
          write(fp,'(a1)') "0"
        end if
      enddo
      enddo
      enddo
      close(fp)
    endif
  endif
  
end subroutine writeavs

