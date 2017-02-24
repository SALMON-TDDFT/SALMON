!
!  Copyright 2016 ARTED developers
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
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120--------130
Subroutine write_GS_data
  use Global_Variables
  use communication, only: comm_is_root
  implicit none
  integer ik,ib,ia,iter,j

  if (comm_is_root()) then
    open(403,file=file_GS)
    write(403,*) '#---------------------------------------------------------'
    write(403,*) '#grid information-----------------------------------------'
    write(403,*) '#aL =',aL
    write(403,*) '#ax,ay,az =',ax,ay,az
    write(403,*) '#aLx,aLy,aLz =',aLx,aLy,aLz
    write(403,*) '#bLx,bLy,bLz =',bLx,bLy,bLz
    write(403,*) '#Nd =',Nd
    write(403,*) '#NLx,NLy,NLz=',NLx,NLy,NLz
    write(403,*) '#NL =',NL
    write(403,*) '#Hx,Hy,Hz =',Hx,Hy,Hz
    write(403,*) '#(pi/max(Hx,Hy,Hz))**2 =',(pi/max(Hx,Hy,Hz))**2
    write(403,*) '#(pi/Hx)**2+(pi/Hy)**2+(pi/Hz)**2 =',(pi/Hx)**2+(pi/Hy)**2+(pi/Hz)**2
    write(403,*) '#Hxyz =',Hxyz
    write(403,*) '#NKx,NKy,NKz=',NKx,NKy,NKz
    write(403,*) '#NKxyz =',NKxyz
    write(403,*) '#Sym=',Sym
    write(403,*) '#NK =',NK
    write(403,*) '#NEwald, aEwald =',NEwald, aEwald 
    write(403,*) '#---------------------------------------------------------'
    write(403,*) '#GS calc. option------------------------------------------'
    write(403,*) '#FSset_option =',FSset_option
    write(403,*) '#Ncg=',Ncg
    write(403,*) '#Nmemory_MB,alpha_MB =',Nmemory_MB,alpha_MB
    write(403,*) '#NFSset_start,NFSset_every =',NFSset_start,NFSset_every
    write(403,*) '#Nscf=',Nscf
    write(403,*) '#NI,NE=',NI,NE
    write(403,*) '#Zatom=',(Zatom(j),j=1,NE)
    write(403,*) '#Lref=',(Lref(j),j=1,NE)
    write(403,*) '#i,Kion(ia)','(Rion(j,a),j=1,3)'
    do ia=1,NI
      write(403,*) '#',ia,Kion(ia)
      write(403,*) '#',(Rion(j,ia),j=1,3)
    end do
    write(403,*) '#---------------------------------------------------------'
    write(403,*) '#GS information-------------------------------------------'
    write(403,*) '#NB,Nelec=',NB,Nelec
    write(403,*) '#Eall =',Eall
    write(403,*) '#dns_diff(iter = Nscf)',dns_diff(Nscf)
    write(403,*) '#esp_var_ave(iter = Nscf)',esp_var_ave(Nscf)
    write(403,*) '#esp_var_max(iter = Nscf)',esp_var_max(Nscf)
    write(403,*) '#NBoccmax is ',NBoccmax
    write(403,*) '#---------------------------------------------------------'
    write(403,*) '#band information-----------------------------------------'
    write(403,*) '#Bottom of VB',minval(esp_vb_min(:))
    write(403,*) '#Top of VB',maxval(esp_vb_max(:))
    write(403,*) '#Bottom of CB',minval(esp_cb_min(:))
    write(403,*) '#Top of CB',maxval(esp_cb_max(:))
    write(403,*) '#Fundamental gap',minval(esp_cb_min(:))-maxval(esp_vb_max(:))
    write(403,*) '#Fundamental gap[eV]',(minval(esp_cb_min(:))-maxval(esp_vb_max(:)))*2.0*Ry
    write(403,*) '#BG between same k-point',minval(esp_cb_min(:)-esp_vb_max(:))
    write(403,*) '#BG between same k-point[eV]',(minval(esp_cb_min(:)-esp_vb_max(:)))*2.0*Ry
    write(403,*) '#Physicaly upper bound of CB for DOS',minval(esp_cb_max(:))
    write(403,*) '#Physicaly upper bound of CB for eps(omega)',minval(esp_cb_max(:)-esp_vb_min(:))
    write(403,*) '#---------------------------------------------------------'
    do iter=1,Nscf
      write(403,'(1x,i5,4e20.10)') iter, Eall_GS(iter),dns_diff(iter),esp_var_ave(iter),esp_var_max(iter)
    end do
    close(403)

    open(404,file=file_DoS)
    open(405,file=file_band)
    write(404,*) '#Occupation distribution at Ground State'
    write(404,*) '#(NK,NB)=','(',NK,NB,')'
    write(405,*) '#Bandmap at Ground State'
    write(405,*) '#(NK,NB)=','(',NK,NB,')'
    do ik=1,NK
      do ib=1,NB
        write(404,'(1x,2i5,2f15.10)') ik,ib,esp(ib,ik),occ(ib,ik)*NKxyz
        write(405,'(1x,2i5,3f15.8,3x,f15.10)') ik,ib,kAc(ik,1),kAc(ik,2),kAc(ik,3),esp(ib,ik)
      enddo
    enddo
    close(404)
    close(405)

  end if

  return
End Subroutine write_GS_data
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120--------130
