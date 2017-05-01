! Copyright 2017 Katsuyuki Nobusada, Masashi Noda, Kazuya Ishimura, Kenji Iida, Maiku Yamaguchi
!
! Licensed under the Apache License, Version 2.0 (the "License");
! you may not use this file except in compliance with the License.
! You may obtain a copy of the License at
!
!     http://www.apache.org/licenses/LICENSE-2.0
!
! Unless required by applicable law or agreed to in writing, software
! distributed under the License is distributed on an "AS IS" BASIS,
! WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
! See the License for the specific language governing permissions and
! limitations under the License.

subroutine calc_Mps3rd
use scf_data
implicit none
integer :: iatom,jj,ix,iy,iz
integer :: icount

numatom_ps=0
do iatom=1,MI
  do jj=1,Mps(iatom)
    if(Jxyz(1,jj,iatom)>=mg_sta(1).and.Jxyz(1,jj,iatom)<=mg_end(1).and.  &
       Jxyz(2,jj,iatom)>=mg_sta(2).and.Jxyz(2,jj,iatom)<=mg_end(2).and.  &
       Jxyz(3,jj,iatom)>=mg_sta(3).and.Jxyz(3,jj,iatom)<=mg_end(3))then
      numatom_ps(Jxyz(1,jj,iatom),Jxyz(2,jj,iatom),Jxyz(3,jj,iatom))=  &
        numatom_ps(Jxyz(1,jj,iatom),Jxyz(2,jj,iatom),Jxyz(3,jj,iatom))+1
    end if
  end do
end do

icount=0
MImax=0
do iz=mg_sta(3),mg_end(3)
do iy=mg_sta(2),mg_end(2)
do ix=mg_sta(1),mg_end(1)
  if(numatom_ps(ix,iy,iz)>0) then
    icount=icount+1
    if(numatom_ps(ix,iy,iz)>MImax)then
      MImax=numatom_ps(ix,iy,iz)
    end if
  end if
end do
end do
end do

allocate(iatomnum_ps(MImax,mg_sta(1):mg_end(1),mg_sta(2):mg_end(2),mg_sta(3):mg_end(3)))
allocate(Jxyz_all(3,icount))

numatom_ps=0
do iatom=1,MI
  do jj=1,Mps(iatom)
    if(Jxyz(1,jj,iatom)>=mg_sta(1).and.Jxyz(1,jj,iatom)<=mg_end(1).and.  &
       Jxyz(2,jj,iatom)>=mg_sta(2).and.Jxyz(2,jj,iatom)<=mg_end(2).and.  &
       Jxyz(3,jj,iatom)>=mg_sta(3).and.Jxyz(3,jj,iatom)<=mg_end(3))then
      numatom_ps(Jxyz(1,jj,iatom),Jxyz(2,jj,iatom),Jxyz(3,jj,iatom))=  &
        numatom_ps(Jxyz(1,jj,iatom),Jxyz(2,jj,iatom),Jxyz(3,jj,iatom))+1
      iatomnum_ps(numatom_ps(Jxyz(1,jj,iatom),Jxyz(2,jj,iatom),Jxyz(3,jj,iatom)),  &
                             Jxyz(1,jj,iatom),Jxyz(2,jj,iatom),Jxyz(3,jj,iatom))=iatom
    end if
  end do
end do

maxMps_all=0
MImax=0
do iz=mg_sta(3),mg_end(3)
do iy=mg_sta(2),mg_end(2)
do ix=mg_sta(1),mg_end(1)
  if(numatom_ps(ix,iy,iz)>0)then
    maxMps_all=maxMps_all+1
    Jxyz_all(1,maxMps_all)=ix
    Jxyz_all(2,maxMps_all)=iy
    Jxyz_all(3,maxMps_all)=iz
   if(numatom_ps(ix,iy,iz)>MImax) MImax=numatom_ps(ix,iy,iz)
  end if
end do
end do
end do

allocate(Mps3rd(MImax,mg_sta(1):mg_end(1),mg_sta(2):mg_end(2),mg_sta(3):mg_end(3)))
numatom_ps=0
do iatom=1,MI
  do jj=1,Mps(iatom)
    if(Jxyz(1,jj,iatom)>=mg_sta(1).and.Jxyz(1,jj,iatom)<=mg_end(1).and.  &
       Jxyz(2,jj,iatom)>=mg_sta(2).and.Jxyz(2,jj,iatom)<=mg_end(2).and.  &
       Jxyz(3,jj,iatom)>=mg_sta(3).and.Jxyz(3,jj,iatom)<=mg_end(3))then
      numatom_ps(Jxyz(1,jj,iatom),Jxyz(2,jj,iatom),Jxyz(3,jj,iatom))=  &
        numatom_ps(Jxyz(1,jj,iatom),Jxyz(2,jj,iatom),Jxyz(3,jj,iatom))+1
      Mps3rd(numatom_ps(Jxyz(1,jj,iatom),Jxyz(2,jj,iatom),Jxyz(3,jj,iatom)),  &
                        Jxyz(1,jj,iatom),Jxyz(2,jj,iatom),Jxyz(3,jj,iatom))=jj
    end if
  end do
end do


allocate(numatom_ps_2nd(maxMps_all))
do jj=1,maxMps_all
  numatom_ps_2nd(jj)=numatom_ps(Jxyz_all(1,jj),Jxyz_all(2,jj),Jxyz_all(3,jj))
end do

jamax=0
do jj=1,maxMps_all
  do iatom=1,numatom_ps(Jxyz_all(1,jj),Jxyz_all(2,jj),Jxyz_all(3,jj))
    jamax=jamax+1
  end do
end do

allocate(Jxyz_all_2nd(3,jamax))
allocate(iatomnum_ps_2nd(jamax))
allocate(Mps3rd_2nd(jamax))
allocate(jja(MImax,maxMps_all))

jamax=0
do jj=1,maxMps_all
  do iatom=1,numatom_ps(Jxyz_all(1,jj),Jxyz_all(2,jj),Jxyz_all(3,jj))
    jamax=jamax+1
    Jxyz_all_2nd(1,jamax)=Jxyz_all(1,jj)
    Jxyz_all_2nd(2,jamax)=Jxyz_all(2,jj)
    Jxyz_all_2nd(3,jamax)=Jxyz_all(3,jj)
    iatomnum_ps_2nd(jamax)=iatomnum_ps(iatom,Jxyz_all(1,jj),Jxyz_all(2,jj),Jxyz_all(3,jj))
    Mps3rd_2nd(jamax)=Mps3rd(iatom,Jxyz_all(1,jj),Jxyz_all(2,jj),Jxyz_all(3,jj))
    jja(iatom,jj)=jamax
  end do
end do
deallocate(Mps3rd)

end subroutine calc_Mps3rd
