C      Included with permission from W. Bachman.
      subroutine wrtab
C      subroutine wrtab5 5/10/01 wjb
      save
C      include '/usr/nmv/nm/NSIZES'
      include 'c:\nmv\nm\NSIZES'
      common /cm1/ ntheta,nth
      common /cm2/ neta,neps,nre2
      common /cm12/ cov(lpar3),invcov(lpar3),sthta(lth),
     1              sen(lvr,lvr),cor(lpar3),nvls,vls(lpar)
      common /cm18/ spec
      common /cm19/ avail(11)
      common /cm23/ theta(lth),varn(lvr,lvr)
      common /cm43/ omnm(lvr,lvr),sgnm(lvr,lvr)
C     N99 is number of estimated parameters
      COMMON /CM10/ S(LPAR),N99,IL,IU
      
      double precision S
      real invcov
      real*8 theta,varn
      character*6 S6(100)
      character*12 dots
      dots='............'
      unitp=42
      unitc=43
      unitv=44
      if(avail(2).eq.0)then
        do 25 i=1,ntheta
          sthta(i)=spec
   25   continue
        do 35 i=1,neta+neps
        do 35 j=1,i
          sen(i,j)=spec
   35   continue
C      write(unitp,*)'Covariance step not implemented or unsuccessful'
C      write(unitc,*)'Covariance step not implemented or unsuccessful'
C      write(unitv,*)'Covariance step not implemented or unsuccessful'
      endif 
      write(unitp,*)
     1 ' Parameter    Estimate       StandardError '
      do 100 i=1,ntheta
         if(i.lt.10)then
           if(sthta(i).ne.SPEC)then
             write(unitp,110)i,theta(i),sthta(i)
           else
             write(unitp,111)i,theta(i),dots
           endif
         else
           if(sthta(i).ne.SPEC)then
             write(unitp,112)i,theta(i),sthta(i)
           else
             write(unitp,113)i,theta(i),dots
           endif
         endif
  100 continue
  110 format('  TH_0',i1,6x,1pe12.5,5x,1pe12.5)
  111 format('  TH_0',i1,6x,1pe12.5,5x,a12)
  112 format('  TH_',i2,6x,1pe12.5,5x,1pe12.5)
  113 format('  TH_',i2,6x,1pe12.5,5x,a12)
      do 200 i=1,neta
      do 200 j=1,i
         if(i.lt.10.and.j.lt.10)then
           if(varn(i,j).ne.SPEC)then
             if(sen(i,j).ne.SPEC)then
                write(unitp,210)i,j,varn(i,j),sen(i,j)
             else
                write(unitp,211)i,j,varn(i,j),dots
             endif
           else
                write(unitp,212)i,j,dots,dots
           endif
         elseif(i.lt.10.and.j.ge.10)then
           if(varn(i,j).ne.SPEC)then
             if(sen(i,j).ne.SPEC)then
                write(unitp,213)i,j,varn(i,j),sen(i,j)
             else
                write(unitp,214)i,j,varn(i,j),dots
             endif
           else
                write(unitp,215)i,j,dots,dots
           endif
         elseif(i.ge.10.and.j.lt.10)then
           if(varn(i,j).ne.SPEC)then
             if(sen(i,j).ne.SPEC)then
                write(unitp,216)i,j,varn(i,j),sen(i,j)
             else
                write(unitp,217)i,j,varn(i,j),dots
             endif
           else
                write(unitp,218)i,j,dots,dots
           endif
         else
           if(varn(i,j).ne.SPEC)then
             if(sen(i,j).ne.SPEC)then
                write(unitp,219)i,j,varn(i,j),sen(i,j)
             else
                write(unitp,220)i,j,varn(i,j),dots
             endif
           else
                write(unitp,221)i,j,dots,dots
           endif
         endif
  200 continue
  210 format('  OM_0',i1,'_0',i1,3x,1pe12.5,5x,1pe12.5)
  211 format('  OM_0',i1,'_0',i1,3x,1pe12.5,5x,a12)
  212 format('  OM_0',i1,'_0',i1,3x,a12,5x,a12)
  213 format('  OM_0',i1,'_',i2,3x,1pe12.5,5x,1pe12.5)
  214 format('  OM_0',i1,'_',i2,3x,1pe12.5,5x,a12)
  215 format('  OM_0',i1,'_',i2,3x,a12,5x,a12)
  216 format('  OM_',i2,'_0',i1,3x,1pe12.5,5x,1pe12.5)
  217 format('  OM_',i2,'_0',i1,3x,1pe12.5,5x,a12)
  218 format('  OM_',i2,'_0',i1,3x,a12,5x,a12)
  219 format('  OM_',i2,'_',i2,3x,1pe12.5,5x,1pe12.5)
  220 format('  OM_',i2,'_',i2,3x,1pe12.5,5x,a12)
  221 format('  OM_',i2,'_',i2,3x,a12,5x,a12)
      do 300 i=1,neps
      do 300 j=1,i
         if(i.lt.10.and.j.lt.10)then
           if(varn(neta+i,neta+j).ne.SPEC)then
             if(sen(neta+i,neta+j).ne.SPEC)then
               write(unitp,310),i,j,varn(neta+i,neta+j),
     1             sen(neta+i,neta+j)
             else
                write(unitp,311),i,j,varn(neta+i,neta+j),dots
             endif
           else
                write(unitp,312),i,j,dots,dots
           endif
         elseif(i.lt.10.and.j.ge.10)then
           if(varn(neta+i,neta+j).ne.SPEC)then
             if(sen(neta+i,neta+j).ne.SPEC)then
                write(unitp,313),i,j,varn(neta+i,neta+j),
     1             sen(neta+i,neta+j)
             else
                write(unitp,314),i,j,varn(neta+i,neta+j),dots
             endif
           else
                write(unitp,315),i,j,dots,dots
           endif
         elseif(i.ge.10.and.j.lt.10)then 
           if(varn(neta+i,neta+j).ne.SPEC)then
             if(sen(neta+i,neta+j).ne.SPEC)then
                write(unitp,316),i,j,varn(neta+i,neta+j),
     1             sen(neta+i,neta+j)
             else
                write(unitp,317),i,j,varn(neta+i,neta+j),dots
             endif
           else
                write(unitp,318),i,j,dots,dots
           endif
         else 
           if(varn(neta+i,neta+j).ne.SPEC)then
             if(sen(neta+i,neta+j).ne.SPEC)then
                write(unitp,319),i,j,varn(neta+i,neta+j),
     1             sen(neta+i,neta+j)
             else
                write(unitp,320),i,j,varn(neta+i,neta+j),dots
             endif
           else
                write(unitp,321),i,j,dots,dots
           endif
         endif 
  300 continue
  310 format('  SG_0',i1,'_0',i1,3x,1pe12.5,5x,1pe12.5)
  311 format('  SG_0',i1,'_0',i1,3x,1pe12.5,5x,a12)
  312 format('  SG_0',i1,'_0',i1,3x,a12,5x,a12)
  313 format('  SG_0',i1,'_',i2,3x,1pe12.5,5x,1pe12.5)
  314 format('  SG_0',i1,'_',i2,3x,1pe12.5,5x,a12)
  315 format('  SG_0',i1,'_',i2,3x,a12,5x,a12)
  316 format('  SG_',i2,'_0',i1,3x,1pe12.5,5x,1pe12.5)
  317 format('  SG_',i2,'_0',i1,3x,1pe12.5,5x,a12)
  318 format('  SG_',i2,'_0',i1,3x,a12,5x,a12)
  319 format('  SG_',i2,'_',i2,3x,1pe12.5,5x,1pe12.5)
  320 format('  SG_',i2,'_',i2,3x,1pe12.5,5x,a12)
  321 format('  SG_',i2,'_',i2,3x,a12,5x,a12)
C write number of estimated parameters to file
C      write(unitp,322) n99
C  322 format('  Number of estimated paramters:',i3)


C write correlations as vector
      if(avail(2).ne.0)then
      write(unitc,*)'    Parameters         Correlation'
      do 400 i=1,ntheta
        i1=i/10
        i2=MOD(i,10)
        S6(i)='TH'//char(i1+48)//char(i2+48)
  400 continue
      i=ntheta
      do 440 ii=1,neta
      do 440 ij=ii,neta
        i=i+1
        i1=ii/10
        i2=MOD(ii,10)
        j1=ij/10
        j2=MOD(ij,10)
        S6(i)='OM'//char(i1+48)//char(i2+48)//char(j1+48)//char(j2+48)
  440 continue
      do 460 ii=1,neps
      do 460 ij=ii,neps
        i=i+1
        i1=ii/10
        i2=MOD(ii,10)
        j1=ij/10
        j2=MOD(ij,10)
        S6(i)='SG'//char(i1+48)//char(i2+48)//char(j1+48)//char(j2+48)
  460 continue
      nsum=ntheta+neta*(neta+1)/2+neps*(neps+1)/2
      k=0
      do 500 i=1,nsum
      do 500 j=1,i
        k=k+1
        if(cor(k).ne.SPEC)then
          write(unitc,600)S6(i),S6(j),cor(k)
        else
          write(unitc,601)S6(i),S6(j),dots
        endif
  500 continue
  600 format('    'a6,':',a6,6x,1pe12.5)
  601 format('    'a6,':',a6,6x,a12)
      endif

C attempt to write covariance vector 5/10/01 wjb
C write correlations as vector
      if(avail(2).ne.0)then
      write(unitv,*)'    Parameters         Covariance'
      do 700 i=1,ntheta
        i1=i/10
        i2=MOD(i,10)
        S6(i)='TH'//char(i1+48)//char(i2+48)
  700 continue
      i=ntheta
      do 740 ii=1,neta
      do 740 ij=ii,neta
        i=i+1
        i1=ii/10
        i2=MOD(ii,10)
        j1=ij/10
        j2=MOD(ij,10)
        S6(i)='OM'//char(i1+48)//char(i2+48)//char(j1+48)//char(j2+48)
  740 continue
      do 760 ii=1,neps
      do 760 ij=ii,neps
        i=i+1
        i1=ii/10
        i2=MOD(ii,10)
        j1=ij/10
        j2=MOD(ij,10)
        S6(i)='SG'//char(i1+48)//char(i2+48)//char(j1+48)//char(j2+48)
  760 continue
      nsum=ntheta+neta*(neta+1)/2+neps*(neps+1)/2
      k=0
      do 800 i=1,nsum
      do 800 j=1,i
        k=k+1
        if(cov(k).ne.SPEC)then
          write(unitv,900)S6(i),S6(j),cov(k)
        else
          write(unitv,901)S6(i),S6(j),dots
        endif
  800 continue
  900 format('    'a6,':',a6,6x,1pe12.5)
  901 format('    'a6,':',a6,6x,a12)
      endif
      end
