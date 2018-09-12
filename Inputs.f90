! MODULE INPUTS
!
! This module contains all of the input variables. In principle, only this file should be adjusted for new runs, also occasionally the calling program as well. The other pieces of code should only be opened to add new features or to fix bugs.
!
! ---------------------------
! Inputs :
! ---------------------------
! PHYSICAL PARAMETERS:
!      systemn : Number of systems in the ensemble
!      dim : Number of sites in each system
!      DELTA : Disorder strength
!      hop : hopping amplitude
!      uSite: Interaction Strength
!      ChemPot : Chemical potential
! CODING PARAMETERS:
!      bond_cutoff : Defines the bond strength seperating weak bonds and strong bonds
!!      bins : number of elements in a number of arrays containing output data (aka number of bins in histograms)
!      ClusterMax : The maximum cluster size to be included in calculating the 'DOS' or 'Potential'
!      DOS : Density of States Array. Element i corresponds to the total contributions with energy
!            in the range [min+(max-min)*i/bins,min+(max-min)*(i+1)/bins]
!      DOS_EMax : Maximum energy of a DOS contribution allowed
!      DOS_EMin : Minimum energy of a DOS contribution allowed
!      


module inputs
  implicit none
  
  ! PHYSICAL PARAMETERS
  integer, parameter :: systemn = 1500000
  integer, parameter :: dim = 6
  real, parameter :: DELTA = 20
  real, parameter  :: hop = -1
  real, parameter :: uSite = 8
  real, parameter :: ChemPot = uSite/2
  
  ! CODING PARAMETERS
  real, parameter :: bond_cutoff = 0.5
  real, parameter :: prune_cutoff = 3
  integer, parameter :: bins = 1000
  integer, parameter :: n_d = systemn*dim
  integer, parameter :: ClusterMax = 6
  integer, parameter :: DOS_MaxCluster = 1!ClusterMax
  real DOS(DOS_MaxCluster,bins)
  real, parameter :: DOS_EMax = 2*uSite+1, DOS_EMin = -DOS_EMax
  !real :: DOS_EMax, DOS_EMin
  real :: DroppedDos !Weight of DOS contributions dropped, might be better to track the weights that have been dropped.
  integer :: PrunedBonds !# of bonds which are ignored due to the ignore_cutoff'
  integer :: NumStrongBonds !# of strong bonds to compare with the number of strong bonds that are ignored
  integer :: StrongestBondsPruned !# of bonds that were pruned where a strong bond was replaced with a weak bond
  integer, parameter :: Pot_ClSize = 8
  real Potential(bins)
  real, parameter :: POT_EMax = Delta/2, POT_EMin = -Pot_EMax
  real, parameter :: Bond_EMax = 20, Bond_EMin = -20
  real, parameter :: Diff_EMax = 3*DELTA, Diff_EMin = -3*DELTA
  real, parameter :: Avg_EMax = DELTA+Usite, Avg_EMin = -Avg_EMax
  real AvgEDist(bins)
  !real :: POT_EMax, POT_EMin, Bond_EMax, Bond_EMin, Diff_EMax, Diff_EMin, Avg_EMax, Avg_EMin
  ! RUN THIS PART?
  logical, private, parameter :: ff = .false., tt = .true.
  logical, parameter :: CalcFracSites = ff ! 
  logical, parameter :: CalcBondStrength = ff
  logical, parameter :: CalcEnergyDiff = ff
  logical, parameter :: CalcDos = tt
  logical, parameter :: CalcPot = ff
  logical, parameter :: CalcAvgEofB = ff

  real :: times


  

end module inputs

  
