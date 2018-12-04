module Tools
  USE Inputs
  implicit none
  SAVE

contains
  
  !     Bin_Data subroutine enables the binning of data into a histogram
  !
  !
  !     ------------------------------Arguments-------------------------------------
  !
  !     HistoData: A 1xN array, where N is the number of bins in the histogram. This contains the histogram
  !
  !     Data: A 1xM array, where M is the number of contributions to be added to HistoData. Each element is the value of a contribution
  !
  !     Weights: A 1xM array containing the weights of each contribution in Data
  !              This argument is optional
  !
  !     Dropped: A variable of type integer that keeps track of the number of contributions outside of a specified range
  !              This argument is optional
  !
  !     Max/Min: The maximum (Max) and minimum (Min) values to be counted. Contributions outside are not counted
  subroutine Bin_Data( HistoData, Data, Weights, Dropped, Max, Min )
    implicit none
    real, dimension(:), intent(in) :: Data
    real, dimension(:), Optional, intent(in) :: Weights
    real, intent(inout) :: HistoData(bins)
    real, intent(in) :: Max
    real, intent(in) :: Min
    real, optional, intent(inout) :: Dropped
    integer Loop1, Loop2
 
   
    If ( Present(Dropped) .and. Present(Weights) ) then
       
       do Loop1 = 1,size(Data)
          if ( Weights(Loop1) .eq. 0 ) CYCLE 
          if ( ( Data(Loop1) .lt. Min ) .or. ( Data(Loop1) .gt. Max ) ) then
             Dropped = Dropped + Weights(Loop1)
          else
             Loop2 = FLOOR(bins*(Data(Loop1) - Min)/(Max - Min) + 1)
             if ( (Loop2 .gt. bins) .or. (Loop2 .lt. 1) ) CYCLE             
             HistoData(Loop2) = HistoData(Loop2) + Weights(Loop1)  
             
          end if
          
       end do
    
    else if ( Present(Dropped) .and. (.not. Present(Weights) ) ) then

       do Loop1 = 1,size(Data) 
          if ( ( Data(Loop1) .lt. Min ) .or. ( Data(Loop1) .gt. Max ) ) then
             Dropped = Dropped + Weights(Loop1)
          else             
             Loop2 = FLOOR(bins*(Data(Loop1) - Min)/(Max - Min) + 1)
             if ( (Loop2 .gt. bins) .or. (Loop2 .lt. 1) ) CYCLE
             HistoData(Loop2) = HistoData(Loop2) + 1
             
          end if
          
       end do

    else if ( (.not. Present(Dropped)) .and. Present(Weights) ) then

       do Loop1 = 1,size(Data) 
          if ( ( Data(Loop1) .lt. Min ) .or. ( Data(Loop1) .gt. Max ) ) CYCLE
          Loop2 = FLOOR(bins*(Data(Loop1) - Min)/(Max - Min) + 1)
          if ( (Loop2 .gt. bins) .or. (Loop2 .lt. 1) ) CYCLE
          HistoData(Loop2) = HistoData(Loop2) + Weights(Loop1)
             
       end do
          
          

    else if ( .not. ( Present(Dropped) .or. Present(Weights) ) ) then

       do Loop1 = 1,size(Data) 
          if ( ( Data(Loop1) .lt. Min ) .or. ( Data(Loop1) .gt. Max ) ) CYCLE
          Loop2 = FLOOR(bins*(Data(Loop1) - Min)/(Max - Min) + 1)
          if ( (Loop2 .gt. bins) .or. (Loop2 .lt. 1) ) CYCLE   
          HistoData(Loop2) = HistoData(Loop2) + 1
          
       end do

    else
       print*, "Bin_Data options exceeded"
       STOP
    end if
    
     

  end subroutine Bin_Data

  subroutine OpenFile( Unum, Name, Contents, Column1, Column2 )
    implicit none
    integer, intent(in) :: Unum
    character(len=*), intent(in) :: Name, Contents, Column1, Column2
    character(len=25) :: filename
    filename = Name
    filename = FileNamer(filename)
    open (unit = Unum, file = "output/"//trim(filename)//"" ,status ='unknown' )
    
    
    write(Unum,*) "#This file contains the following data type for the ensemble: ", Trim(Contents)
    write(Unum,*) "#First column is: ", Trim(Column1)
    write(Unum,*) "#Second column is: ", Trim(Column2)
    write(Unum,*) "#Disorder Strength = ", DELTA
    write(Unum,*) "#Bond Cutoff = ", bond_cutoff
    write(Unum,*) "#Pruned bonds cutoff = ", prune_cutoff
    write(Unum,*) "#Fraction of n.n. sites in which at least one possible bond (of three) was pruned = ", &
         real(PrunedBonds)/real(n_d)
    write(Unum,*) "#Fraction of prunings that are useful = ", real(StrongestBondsPruned)/real(PrunedBonds)
    write(Unum,*) "#Hopping = ", hop
    write(Unum,*) "#Interaction strength = ", uSite
    write(Unum,*) "#Chemical Potential = ", ChemPot
    write(Unum,*) "#Dimensions = ", dim
    write(Unum,*) "#Number of systems = ", systemn

  end subroutine OpenFile
  
  !
  ! Prints data to file with unit number, Unum. Inputs: Unum = unit number for file, Form = format for the data, Bins array: xaxis, Data = Data.
  ! DoSMin + (DoSMax-DoSMin)*iii/real(bins)
  
  subroutine PrintData( Unum, Form, Min, Max, BinNum, Data_sp, Data_dp, Dropped )
    implicit none
    integer, intent(in) :: Unum, BinNum
    character(len=*), intent(in) :: Form
    real, intent(in), dimension(:), optional :: Data_sp
    real(dp), intent(in), dimension(:), optional :: Data_dp
    real, intent(in) :: Min, Max
    real, optional, intent(in) :: Dropped
    real :: BinWidth, sum1
    integer :: Loop1

    if ( Present(Data_dp) ) then
       BinWidth = (Max - Min)/real(BinNum)               ! Area of each bin (bins have equal width)
       do Loop1=1,BinNum
          write(Unum,(Form)) Min + ( Max - Min )*Loop1/real(BinNum) - BinWidth/2, Data_dp(Loop1)/( Sum(Data_dp)*BinWidth )
       end do
       sum1 = 0.0
       do Loop1=1,ClusterMax
          sum1 = sum1 + Data_dp(Loop1)/( Sum(Data_dp)*BinWidth )
       end do
       print*, "For W, U, B:", DELTA, uSite, bond_cutoff
       print*, "Total fraction of sites in 4-site clusters or less is: ", sum1
       print*, ""

       
    else if ( Present(Dropped) ) then
       BinWidth = (Max - Min)/real(BinNum)
       do Loop1=1,BinNum
          write(Unum,( Form)) Min + ( Max - Min )*Loop1/real(BinNum) - BinWidth/2, &
               Data_sp(Loop1)/( (Sum(Data_sp) + Dropped) * BinWidth )
       end do
    else
       
       BinWidth = (Max - Min)/real(BinNum)               ! Area of each bin (bins have equal width)
       print*, "BinWidth", BinWidth
       do Loop1=1,BinNum
          write(Unum,(Form)) Min + ( Max - Min )*Loop1/real(BinNum) - BinWidth/2, Data_sp(Loop1)/( Sum(Data_sp)*BinWidth )
       end do
       sum1 = 0.0
       do Loop1=1,ClusterMax
          sum1 = sum1 + Data_sp(Loop1)/( Sum(Data_sp)*BinWidth )
       end do
       print*, "For W, U, B:", DELTA, uSite, bond_cutoff
       print*, "Total fraction of sites in 4-site clusters or less is: ", sum1
       print*, ""
    end If
    
    Close(Unum)

  end subroutine PrintData

  subroutine resize_array(array, numRemove, Start)
    real, dimension(:), allocatable :: tmp_arr
    real, dimension(:), allocatable, intent(inout) :: array
    integer, intent(in) :: numRemove                           ! Number of elements to remove from 'array'
    integer, intent(in) :: Start                               ! The element in which to start removing the array
    integer i, j                                               ! Looping integer
    

    allocate(tmp_arr(size(array) - numRemove))

    tmp_arr(:) = 0

    j=1

    if ( (Start + numRemove) .le. size(array) ) then


       do i = 1, Start - 1

          tmp_arr(i) = array(i)
          j = j + 1

       end do


       do i = (Start + numRemove), Size(array)

          tmp_arr(j) = array(i)
          j = j + 1

       end do



    else if ( (Start + numRemove) .gt. size(array) ) then

       do i = (Start + numRemove - size(array)) , Start-1

          tmp_arr(j) = array(i)
          j = j + 1

       end do


    end if

    deallocate(array)
    allocate(array(size(tmp_arr)))
    
    array = tmp_arr
    
   
  end subroutine resize_array
  
  character(len=20) function str(k)                  ! function to change an integer to a string, used in FileNamer subroutine
     integer, intent(in) :: k

     write (str, *) k
     str = adjustl(str)
   end function str

   character(len=25) function FileNamer(Name)                                ! Increments filename labels Blah023.dat label:023
     character(len=25), intent(in) :: Name                                   ! Initial filename. Would be 'Blah' from above example
     integer num                                                             ! Keeps track of the label increment as this continues to check for name existence
     logical Exist                                                           ! Used to store if a filename exists

     Exist = .true.
     num = 0
     FileNamer = ""//trim(Name)//"000.text"                                  ! Takes Blah and puts it into the format 'Blah000.text'
  
     do while (Exist)                                                        ! Continue to change filename if the current one exists
        inquire( file = "output/"//trim(FileNamer)//"", exist = Exist )                           ! Checks the current folder (where this .f90 file is found) if the current filename exists
        if ( Exist ) then                                                    ! If the file name exists
           num = num + 1                                                     ! Increment the label

           if ( len(trim(str(num))) .eq. 1 ) then                            ! Labels can be from 000 to 999. To preserve filename length the leading zeros are required
              FileNamer = ""//trim(Name)//"00"//trim(str(num))//".text"      ! This if statement (and the one below) check the length of the label (1 has length 1, 23 has length 2, 450 has length 3)
           else if ( len(trim(str(num))) .eq. 2) then                        ! This determines the number of leading zeros to add to the filename label
              FileNamer = ""//trim(Name)//"0"//trim(str(num))//".text"
           else
              FileNamer = ""//trim(Name)//trim(str(num))//".text"
           end if
           
        else if ( .not. Exist ) then                                         ! If it does not exist. Note: Didn't increment 'num'
           if ( len(trim(str(num))) .eq. 1 ) then
              FileNamer = ""//trim(Name)//"00"//trim(str(num))//".text"
           else if ( len(trim(str(num))) .eq. 2) then
              FileNamer = ""//trim(Name)//"0"//trim(str(num))//".text"
           else
              FileNamer = ""//trim(Name)//trim(str(num))//".text"
           end if
           Exist = .false.
        else
           print*, "You belong in a museum. Error: FileNameErr"              ! Outputs error if something strange happens. Can't imagine what
        end if
     end do
   end function FileNamer

 

 end module Tools
 
 
 
