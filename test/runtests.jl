using Codes, Test

@testset "Test binary-based dot product" begin
    @test binarydot(0b01, 0b01) == 1 
    @test binarydot(0b11, 0b11) == 0 
    @test binarydot(0b10, 0b01) == 0
end;

@testset "Test biterror function" begin
    @test biterror(0b00, 0b00) == 0
end;

