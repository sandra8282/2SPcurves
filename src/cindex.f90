!------------------------------------------------------------
!     n             I  the # of participants
!     nevent        I  number of event types (M) not including censoring
!     time          I  event time
!     event         I  event type 1 to M
!     group         I  group - treatment or control
!     w             I  censoring weight
!     cind          O  naive estimator for concordance
!     cwind         O  weighted estimator for concordance
!------------------------------------------------------------
!  Author: Sandra Castro-Pearson
!  Division of Biostatistics, School of Public Health
!  University of Minnesota
!  Latest update: March 24, 2022
!------------------------------------------------------------
!     integer n, nevent, e, ew, i, j
!      double precision tau
!      double precision event(n)
!      double precision group(n)
!      double precision time(n)
!      double precision w(n)
!      double precision cind(nevent)
!      double precision cwind(nevent)
!      double precision numcount
!      double precision dnomcount

subroutine cindex(n,nevent,time,event,group, w, cind, cwind, tau)

   implicit none
   INTEGER :: n, nevent, i, j, e, ew
   INTEGER, DIMENSION(n) :: event, group
   DOUBLE PRECISION, DIMENSION(n) :: time, w
   DOUBLE PRECISION :: tau, cind, cwind, cnew, cwnew, dnomcount, numcount

   e=1

!  #get cindex
   dnomcount=0.
   numcount=0.
      do i=1,n
         do j=1,n
            if (i.eq.j) then
              dnomcount=dnomcount
              numcount=numcount
            else
              if (group(i).eq.1.and.time(i).le.tau.and.group(j).eq.0.and.time(j).le.tau) then
                  dnomcount = dnomcount + 1
                  if (time(j).le.time(i).and.event(i).eq.e.and.event(j).eq.e) then
                     numcount = numcount + 1
                  endif
              endif
            endif
         end do
      end do
   cind = numcount/dnomcount

!  #get weighted cindex
   dnomcount=0.
   numcount=0.
      do i=1,n
         do j=1,n
            if (i.eq.j) then
              dnomcount=dnomcount
              numcount=numcount
            else
              if (group(i).eq.1.and.time(i).le.tau.and.group(j).eq.0.and.time(j).le.tau) then
                  dnomcount = dnomcount + (1/(w(i)*w(j)))
                  if (time(j).le.time(i).and.event(i).eq.e.and.event(j).eq.e) then
                     numcount = numcount + (1/(w(i)*w(j)))
                  endif
              endif
            endif
         end do
      end do
   cind = numcount/dnomcount

end

