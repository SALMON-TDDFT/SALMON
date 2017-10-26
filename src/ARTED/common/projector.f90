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
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120-------130
module projector
  implicit none

contains
  subroutine update_projector(kac_in)
    use Global_Variables
    real(8),intent(in) :: kac_in(NK,3)
    integer :: ik, ia, j, i, ix, iy, iz
    real(8) :: kr

    do ik=NK_s,NK_e
      do ia=1,NI
        do j=1,Mps(ia)
          i=Jxyz(j,ia); ix=Jxx(j,ia); iy=Jyy(j,ia); iz=Jzz(j,ia)
          kr=kac_in(ik,1)*(Lx(i)*Hx-ix*aLx)+kac_in(ik,2)*(Ly(i)*Hy-iy*aLy)+kac_in(ik,3)*(Lz(i)*Hz-iz*aLz)
          ekr_omp(j,ia,ik)=exp(zI*kr)
         end do
       end do
    end do

    do ilma=1,Nlma
      ia=a_tbl(ilma)
      do j=1,Mps(ia)
        zproj(j,ilma,ik) = conjg(ekr_omp(j,ia,ik))*uv(j,ilma)
      end do
    end do


  end subroutine update_projector
end module projector
