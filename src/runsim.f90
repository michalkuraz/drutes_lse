module runsim

  contains
  
    subroutine computehydro()
    
      do
        time = time + time_step
        
        call timestamps%fill(time)
    
    end subroutine computehydro

end module runsim
