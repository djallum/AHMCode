module DOSsetup
  USE Inputs
  USE Tools
  USE PreAnalysis
  USE Diag, only: DiagCluster, PreSetUp
  implicit none
  SAVE




contains

  subroutine GetPotential( SitePotential, weakL, weakR, ClusterSize, SitesRemoved )
    implicit none
    real, dimension(:), intent(in) :: SitePotential
    real, dimension(ClusterSize) :: Sites
    integer, intent(in) :: ClusterSize, weakL, weakR, SitesRemoved
    integer EndSite, i
    if ( weakL .gt. weakR) then
       EndSite = (dim - sitesremoved) - weakL
       Sites(1:EndSite) = SitePotential((weakL+1):(dim - sitesremoved))
       Sites((EndSite+1):ClusterSize) = SitePotential(1:weakR)
    else if ( weakL .eq. weakR) then
       Sites = SitePotential
    else
       Sites = SitePotential((weakL+1):weakR)
    end if

       
    CALL Bin_Data(HistoData = Potential, Data = Sites, Max = POT_EMax, Min = POT_EMin)

    

  end subroutine GetPotential
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! GetDos routine. Inputs: Site Potentials, Cluster size, the labels for the left/right weak bonds of the cluster
  ! Outputs: The effecitve Hoppings between neighbour sites, contributions/weights to the density of states, SitesRemoved/Ignored
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  subroutine GetDos( SitePotential, Teff, ClusterSize, weakL, weakR, SitesRemoved, SitesIgnored )
    implicit none
    
                !---------------------------Inputs-------------------------------------------
    real, dimension(:), intent(in) :: SitePotential
    integer, intent(in) :: ClusterSize
    integer, intent(in) :: weakL, weakR

                !---------------------------Outputs------------------------------------------
    real, dimension(:), allocatable :: Cont
    real, dimension(:), allocatable :: Weight
    real, intent(out) :: Teff
    integer, intent(inout) :: SitesRemoved
    integer, intent(inout) :: SitesIgnored

                !---------------------------Programming Variables----------------------------
    real w1, w2, w3                                      ! The three possible grand potentials for a single site cluster
    integer Llabel                                       ! Must change weakL, relabel it so weakL doesn't have to have intent:inout
    real, dimension(ClusterSize) :: Sites
    integer :: EndSite
    integer :: i, j
    Llabel = weakL
    if ( ClusterSize .eq. dim+1 ) then
                !---------------------------Single Site Clusters-----------------------------
       SitesRemoved = SitesRemoved + 1                   ! Increment the number of sites removed
       

       if ( ( weakL .gt. weakR) ) then                      ! If weakL > weakR the cluster goes over the boundary, therefore weakL = weakR
          Llabel = weakR - ClusterSize                                      ! - ClusterSize
       end if
       if ( Size(SitePotential) .eq. 1 ) then
          Llabel = 0
       end if
       
       
       w1 = 0                                            ! Grand potential of unoccupied state
       w2 = SitePotential(Llabel+1) - ChemPot             !   "       "     of singly occupied state
       w3 = 2*SitePotential(Llabel+1) - 2*ChemPot + uSite !   "       "     of the doubly occupied state 
       
                 !---------------------------Determine the contributions/weights--------------
       ! Cont/Weight are allocatable because number of contributions depends on which ground state we are in.
       if ( (w1 .le. w2) .and. (w1 .le. w3) ) then       ! If w1 is the ground state
          allocate(Cont(1))
          allocate(Weight(1))

          ! Adding an electron when going from unoccupied state to singly occupied
          ! Convention for adding an electron is "Final grand potential" - "Initial Grand Potential"
          Cont = w2 - w1
          Weight = 1
          
       else if ( (w2 .le. w1) .and. (w2 .le. w3) ) then
          allocate(Cont(2))
          allocate(Weight(2))
  
          ! Removing an electron when going from singly occupied state to unoccupied
          ! Convention for removing an electron is "Initial Grand Potential" - "Final grand potential" 
          Cont(1) = w2 - w1
          Weight(1) = 0.5

          ! Adding an electron when going from singly occupied state to doubly occupied
          ! Convention for adding an electron is "Final grand potential" - "Initial Grand Potential"
          Cont(2) = w3 - w2
          Weight(2) = 0.5
          
       else if ( (w3 .le. w1) .and. (w3 .le. w2) ) then
          allocate(Cont(1))
          allocate(Weight(1))

          ! Removing an electron when going from doubly occupied state to singly occupied
          ! Convention for removing an electron is "Initial Grand Potential" - "Final grand potential"
          Cont = w3 - w2
          Weight = 1
          
       else
          print*, 'Error in DOS for cluster size 1'      ! Can occur if w1, w2 or w3 are NaN
       end if
       
                !---------------------------Assign Effective Hopping-------------------------
       Teff = 0
       
                !---------------------------Larger than single site clusters-----------------
    else if ( ( ClusterSize .ge. 1 ) .and. ( ClusterSize .le. ClusterMax ) ) then
       allocate(Cont(2*ClusterSize*(4**ClusterSize)))
       allocate(Weight(2*ClusterSize*(4**ClusterSize)))
	
       Cont = 0.0
       Weight = 0.0

       if ( (weakL .gt. weakR) .and. ( ClusterSize .gt. 1 ) ) then
          EndSite = (dim - sitesremoved) - weakL
          Sites(1:EndSite) = SitePotential((weakL+1):(dim - sitesremoved))
          Sites((EndSite+1):ClusterSize) = SitePotential(1:weakR)
       else if ( (weakL .gt. weakR) .and. ( ClusterSize .eq. 1 ) ) then
          Sites = SitePotential(weakR:weakR)
       else if ( size(SitePotential) .eq. ClusterSize ) then
          Sites = SitePotential
       else
          Sites = SitePotential((weakL+1):weakR)
       end if
       CALL DiagCluster(ClusterSize,4**ClusterSize,Sites,Cont,Weight,ChemPot,uSite,hop)

!       CALL Bin_Data(DOS, Cont, Weight, DroppedDos, DOS_EMax, DOS_EMin)

       if (DOS_MaxCluster .eq. ClusterMax) then
          do i=ClusterSize,DOS_MaxCluster
             if (i .eq. ClusterMax) then
                CALL Bin_Data(DOS(i,:), Cont, Weight, DroppedDos, DOS_EMax, DOS_EMin)
             else 
                CALL Bin_Data(HistoData = DOS(i,:), Data = Cont, Weights = Weight &
                     , Max = DOS_EMax, Min = DOS_EMin )
             end if
          end do
       else if (.not. (DOS_MaxCluster .eq. ClusterMax)) then
          CALL Bin_Data(DOS(1,:), Cont, Weight, DroppedDos, DOS_EMax, DOS_EMin)
       end if
       
       deallocate(Cont, Weight)
       Teff = 0
    else if ( ClusterSize .gt. clusterMax ) then                
       SitesIgnored = SitesIgnored + ClusterSize         ! Increment the number of sites ignored, happens when cluster is too large

       ! Effective Hopping still 0
       Teff = 0
       
    end if
  end subroutine GetDos

  subroutine MinCouplingLoc( WeakBonds, Bonds, ClusterStage, SitesRemoved, weakL, weakR)
    implicit none
   
    integer, dimension(:), intent(in) :: WeakBonds
    real, dimension(:), intent(in) :: Bonds
    integer, dimension(:), allocatable :: Temp
    real, dimension(:,:), allocatable :: Order
    integer, intent(out) :: weakL, weakR
    integer, intent(in) :: ClusterStage
    integer, intent(in) :: SitesRemoved
    

    integer numCluster
    integer nextWeak
    integer loop1, loop2, row
    integer ClusterSize

    allocate(Temp(Size(WeakBonds)))
    Temp = 0

    do loop1 = 1,size(WeakBonds)-1
       Temp(loop1 + 1) = WeakBonds(loop1)
    end do

    Temp(1) = WeakBonds(Size(WeakBonds)) - (dim - SitesRemoved)

    numCluster = Count(abs( Temp - WeakBonds ) .eq. ClusterStage)

    if (numCluster .gt. 0) then
       
       allocate(Order(numCluster,3))
       
       loop2 = 1
       do loop1 = 1,Size(WeakBonds) 
          weakL = WeakBonds(Loop1)                                  ! Currently cluster has weakL as its left neighbour
          
          nextWeak = Loop1 + 1                                      ! Next bond label in WeakBonds creates a cluster with weakL
          if ( (Loop1 + 1) .gt. size(WeakBonds) ) then              ! If weakL is the last weakbond in the system then search the beginning of the WeakBonds array
             nextWeak = 1                                           
             weakR = WeakBonds(nextWeak)                            ! weakR is the first weak bond in the system
             ClusterSize = weakR - ( weakL - (dim - SitesRemoved) ) ! Have to calculate the cluster size in this way due to boundary conditions
          else                                                      ! For all other cases, the weakR label is just the next element in the weakbonds array
             weakR = WeakBonds(nextWeak)
             ClusterSize = weakR - weakL
          end if
          
          
          if ( ClusterSize .eq. ClusterStage ) then
             Order(loop2, 1) = Bonds(weakL)**2 + Bonds(weakR)**2
             Order(loop2, 2) = weakL
             Order(loop2, 3) = weakR
             loop2 = loop2 + 1
          end if
          
          
       end do

       row = minloc( Order( 1:numCluster, 1 ), dim=1 ) 
       

       

       weakL = int(Order(row,2)); weakR = int(Order(row,3))

       deallocate(Order,Temp)
    else if ( numCluster .eq. 0 ) then
       weakL = 0; weakR = 0
    else
       print*, "Error in MinCouplingLoc", numCluster
       stop
    end if
    
    
  end subroutine MinCouplingLoc
  

  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! Full_DoS routine. Inputs: Site Potentials, Hopping Amplitudes, Bond Strengths, Weak bonds array, Density of states parameters
  ! Outputs: Density of states for the system, number of contributions dropped from the density of states, Number of sites not
  ! renormalized
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  subroutine System_DoS( SitePotential, Hopping, Bonds, WeakBonds, SitesIgnored )
    implicit none

                !---------------------------Inputs-------------------------------------------
    real, dimension(:), allocatable, intent(inout) :: SitePotential    
    real, dimension(:), allocatable, intent(inout) :: Hopping          
    real, dimension(:), allocatable, intent(inout) :: Bonds            
    integer, dimension(:), allocatable, intent(inout) :: WeakBonds
    
                !---------------------------Outputs------------------------------------------                                                    
    integer, intent(inout) :: SitesIgnored
    

                !---------------------------Programming Variables----------------------------
    integer SitesRemoved                                 ! Keeps track of the number of sites removed at each stage of the RG. Done to
                                                         ! change lengths of loops over the system size
    integer weakL                                        ! The label for the weak bond on the "left" of the current cluster
    integer weakR                                        ! The label for the weak bond on the "left" of the current cluster
    real Teff                                            ! Effecitve hopping amplitude between neighbouring sites after a cluster is
                                                         ! removed between them
    integer ClusterSize                                  ! Size of the current cluster
    integer IndexOfNeighbour, ClusterStage, Loop2, i, EndSite              ! Loop integers
    real :: start, finish

    real, allocatable :: sites(:)
    integer :: EndSites
    logical :: SkipDOS
                !---------------------------Initializing-------------------------------------
    
    SitesRemoved = 0
    weakL = 0
    weakR = 0
    

    ClusterStage = 1                                                          ! Current Cluster size, starting from 1 going to system size
    Loop2 = 0
    Teff = 0

    do while ( (size(WeakBonds) .gt. 1) .and. (WeakBonds(1) .ne. 0) .and. ( ClusterStage .le. ClusterMax ) ) ! Stop when WeakBonds is set to have 1 element or have the first
                                                                       ! element be 0. Either only one cluster left or no weak bonds.
 
                !---------------------------Determine Cluster Size and the bond labels-------
       CALL MinCouplingLoc( WeakBonds, Bonds, ClusterStage, SitesRemoved, weakL, weakR)
       
       if ( weakL .gt. weakR ) then                              ! If weakL is the last weakbond in the system then search the beginning of the WeakBonds array
          ClusterSize = weakR - ( weakL - (dim - SitesRemoved) ) ! Have to calculate the cluster size in this way due to boundary conditions
       else if ( weakL .lt. weakR) then                          ! For all other cases, the weakR label is just the next element in the weakbonds array
          ClusterSize = weakR - weakL
       else if ( (weakL .eq. 0) .and. (weakR .eq. 0) .and. ( (ClusterStage + 1) .le. ClusterMax) ) then
          ClusterSize = 0
          ClusterStage = ClusterStage + 1
          !Call PreSetUp(ClusterStage)
          CYCLE
       else if ( ( weakL .eq. 0 ) .and. (weakR .eq. 0) ) then
          ClusterStage = ClusterStage + 1
          exit
       else
          Print*, "weakL = weakR, but not 0 in Full_DoS"
          print*, WeakBonds
          print*, weakL, weakR
          stop
       end if
                !---------------------------Perform RG is cluster is of smallest size remaining-----------
       ! Search for smallest clusters first, labeled as ClusterStage, it increases when clusters of a certain size are eliminated
       if ( ClusterSize .eq. ClusterStage ) then                                 ! If the current cluster is small enough then remove it                      
          Loop2 = 0
          !---------------------------Calculate DoS for cluster and bin the contributions-----------
          SkipDOS = .false.
          allocate(Sites(ClusterSize))
          if ( (weakL .gt. weakR) .and. ( ClusterSize .gt. 1 ) ) then
             EndSite = (dim - sitesremoved) - weakL
             Sites(1:EndSite) = SitePotential((weakL+1):(dim - sitesremoved))
             Sites((EndSite+1):ClusterSize) = SitePotential(1:weakR)
          else if ( (weakL .gt. weakR) .and. ( ClusterSize .eq. 1 ) ) then
             Sites = SitePotential(weakR:weakR)
          else if ( size(SitePotential) .eq. ClusterSize ) then
             Sites = SitePotential
          else
             Sites = SitePotential((weakL+1):weakR)
          end if
          
          do i=1,ClusterSize-1
             IndexOfNeighbour = i+1

             if ( abs(Sites(i) - Sites(IndexOfNeighbour)-uSite) .lt. drop_cutoff ) then
                SkipDOS = .true.
             else if ( abs(Sites(i) - Sites(IndexOfNeighbour)+uSite) .lt. drop_cutoff ) then
                SkipDOS = .true.
             end if
          end do
          deallocate(Sites)
                
          If ( CalcDos .and. (.not. SkipDOS) ) CALL GetDoS( SitePotential &
               , Teff, ClusterSize, weakL, weakR, SitesRemoved, SitesIgnored )               ! Call the Get_DoS routine to calculate the contributions and their weights to DOS, weakL, weakR, SitesRemoved, SitesIgnored )
          If ( CalcPot ) &
               CALL GetPotential( SitePotential, weakL, weakR, ClusterSize, SitesRemoved )

          SitesRemoved = SitesRemoved + ClusterSize
                !---------------------------Removing the cluster from the relevant data-------------------

          CALL resize_array( SitePotential, ClusterSize, weakL + 1 )          ! Remove the sites in the cluster used above which has ClusterSize sites and is bounded on the left by weakL

          CALL resize_array( Hopping, ClusterSize, weakL + 1 )                ! Same process as above, but instead for the Hopping amplitude

          if ( weakL .gt. Size(Hopping) ) then                            ! weakL is used for the last cluster. If the cluster in question is around the boundary and weakL referenced
             Hopping(Size(Hopping)) = Teff                            ! a bond that has a label larger than the new system size, the effective bond between sites neighbouring the
          else                                                            ! removed cluster has a label of weakL-ClusterSize
             Hopping(weakL) = Teff                                        ! Otherwise the new bond has the same value as weakL
          end if

          CALL resize_array( Bonds, ClusterSize, weakL + 1 )                  ! Resize the Bonds array


          !do loop5 = 1,( dim - SitesRemoved )                             ! Use this loop to recalculate the Bonds array. This must be done before of Teff changing that single bond
             IndexOfNeighbour = weakL+1                                   ! In early stages of this code, teff will always be zero so could just set that same element to 0 in the bonds
             ! array. But in the future teff will be nonzero.
             if ( weakL .ge. ( dim - SitesRemoved) ) then                  ! Could also consider only recalculating the single bond but this current set up works and this is not the
                IndexOfNeighbour = 1                                      ! time consuming part of the code.
                weakL = dim - SitesRemoved
             end if                                                       ! Recalculating the whole thing as allowed me to catch an error in the past.
             
             CALL Strongest_Bond( Bonds, SitePotential, Hopping, weakL, & ! Recalculate the loop5'th strongest bond.
                  IndexOfNeighbour )
             
          !end do
          
          
          deallocate(WeakBonds)                                           ! deallocate the weakbonds array

          CALL CalcWeakBonds(Bonds,WeakBonds,SitesRemoved)                ! recalculate the weakbonds array, using the newly calculated Bonds array

          
          
       else if ( ClusterSize .gt. ClusterStage ) then
          print*, ClusterStage, weakL, weakR
          stop
       else if ( (ClusterSize .lt. ClusterStage) .and. (ClusterSize .ne. 0) ) then
          print*, "Small Clusters Found Late"                             ! This implicitly sets weakL to weakR and weakR to weakR+1
          print*, WeakBonds
          stop
       else if ( ClusterSize .eq. 0 ) then
          !If the cluster size is 0 then there are no more clusters of that size in the RG
       else
          print*, "Error in Full_Dos"
          stop
       end if

       
       
    end do
       
    ClusterSize = dim - SitesRemoved
    if ( (ClusterSize .le. ClusterMax) .and. (ClusterStage .ne. ClusterSize) ) then
       !CALL PreSetUp(ClusterSize)
       weakL = WeakBonds(1); weakR = WeakBonds(1)
       if ( CalcDos ) CALL GetDoS( SitePotential, Teff, ClusterSize &
            , weakL, weakR, SitesRemoved, SitesIgnored )
    
       If ( CalcPot ) &
                  CALL GetPotential( SitePotential, weakL, weakR, ClusterSize, SitesRemoved )
    end if
    !CALL PreSetUp(1)    
        
    
  end subroutine System_DoS

  subroutine Full_Dos()
    implicit none
    real, dimension(:), allocatable :: SitePotential, Hopping, Bonds
    integer, dimension(:), allocatable :: WeakBonds

    integer SitesRemoved
    integer SitesIgnored                                ! Keeps track of all the sites in clusters larger than this code can handle
    
    
    !----------------------------Coding Tools-------------------------------------
    integer i                                     ! Loop Integer
    real :: start, finish, TIME

 
    DOS = 0.0
    Potential = 0.0
    DroppedDos = 0
    PrunedBonds = 0
    

    call init_random_seed()
  
    SitesIgnored = 0
    SitesRemoved = 0

    CALL CPU_TIME(start)
    
    do i = 1,systemn
       if ( mod(real(i),0.1*systemn) == 0.0 ) then
          print*, int(i/real(systemn)*100), "%"
       end if
       
       print*, "hello"
       !---------------------------Create System------------------------------------
       allocate(SitePotential(dim),Hopping(dim),Bonds(dim))
       call create_AHM( SitePotential, Hopping, Bonds )
       call CalcWeakBonds( Bonds, WeakBonds, SitesRemoved )
       if ( WeakBonds(1) .eq. 0 ) then
          !print*, "All Bonds Strong"
          deallocate(SitePotential,Hopping,Bonds,WeakBonds)
          CYCLE
       end if
       CALL system_Dos( SitePotential, Hopping, Bonds, WeakBonds, SitesIgnored )
       deallocate(SitePotential, Hopping, Bonds, WeakBonds)
    end do

    CALL CPU_TIME(finish)
    TIME = finish - start
    !print*, "Sites in large clusters: ", SitesIgnored, " out of", n_d
    !print '("Time of New= ", I1.1," days,", I2.1," hours,", I2.1," minutes, and ", I2.1," seconds.")', &
    !int(TIME/real(86400)), mod(int(TIME/real(3600)),24), mod(int(TIME/real(60)),60), mod(int(TIME),60)

    
    If ( CalcDos ) then
       if ( .not. (DOS_MaxCluster .eq. ClusterMax) ) then
          CALL OpenFile(100, "DOS", "Density of States", "Energy", "Density of States" )
          write(100,*) "#Maximum cluster Size included: ", ClusterMax
          write(100,*) "#Time (s) = ", TIME 
          CALL PrintData(100, '(g12.5,g12.5)', DOS_EMin, DOS_EMax, bins, DOS(1,:), DroppedDos)
          Close(100)
       else
          do i=1,ClusterMax
             CALL OpenFile(100+i, "DOSInc"//trim(str(i))//"_", "Density of States", "Energy", "Density of States" )
             write(100+i,*) "#Maximum cluster Size included: ", i
             write(100+i,*) "#Time (s) = ", TIME
             !if (i .eq. 1) print*, DOS(i,1:100)
             DOS(i,:) = DOS(i,:)/(Sum(DOS(DOS_MaxCluster,:)) + DroppedDos)
             
             
             CALL PrintData(100+i, '(g12.5,g12.5)', DOS_EMin, DOS_EMax, bins, DOS(i,:))
             Close(100+i)
          end do
       end if      
    end If
    print*, "times", times
    If ( CalcPot ) then
       CALL OpenFile(200, "POT"//trim(str(Pot_ClSize))//"_", "Distribution of Site Potentials in counted clusters", &
            "Site Potentials", "Height of distribution" )
       write(200,*) "#Maximum cluster Size included: ", ClusterMax
       write(200,*) "#Time (s) = ", TIME
       Potential = Potential/Sum(Potential)
       print*, Potential
       CALL PrintData(200, '(g12.5,g12.5)', POT_EMin, POT_EMax, bins, Potential)
       Close(200)
    end If

  end subroutine Full_Dos

  end module DOSsetup
