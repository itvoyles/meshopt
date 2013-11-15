module select_precision 

! use selected_real_kind to select desired kinds of real variables in 
! a processor-independent manner 

implicit none

save

! Declare Parameters:
integer, parameter :: sngl  = selected_real_kind(p=6,  r=37)
integer, parameter :: dbl   = selected_real_kind(p=13, r=200)
integer, parameter :: extnd = selected_real_kind(p=17, r=2000)
integer, parameter :: quad  = selected_real_kind(p=26, r=200)
integer, parameter :: prec = dbl 

end module select_precision
