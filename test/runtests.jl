using Codes, Test

@test binarydot(0b01,0b01) == 1 
@test binarydot(0b11,0b11) == 0 
@test binarydot(0b10,0b01) == 0
