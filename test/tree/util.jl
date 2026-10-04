module TestUtil
using CrystalShift
using CrystalShift: find_first_unassigned, cos_angle
using CrystalShift: precision, recall, get_phase_number
using CrystalShift: get_ground_truth, top_k_accuracy, in_top_k

using Test

t = Vector{CrystalPhase}(undef,10)

@testset "Helper functions" begin
    @test cos_angle([1,2,3], [1,2,3]) == 1
    @test find_first_unassigned(t) == 1
end

@testset "model evaluation utils" begin
    answer = Array{Int64}(undef, (7, 6))
    gt = Array{Int64}(undef, (7, 5))
    _answer = [[1,0,0,0,0,1],
            [0,1,0,0,0,0],
            [0,0,1,0,0,0],
            [0,0,0,1,0,1],
            [0,0,0,0,1,1],
            [0,1,0,1,0,0],
            [1,0,1,0,1,0]]

    _gt     = [[1,0,0,0,0],
            [0,0,1,0,0],
            [0,1,0,0,0],
            [0,0,0,1,0],
            [0,0,0,0,1],
            [0,1,0,1,0],
            [1,0,0,0,1]]

    for i in 1:size(answer, 1)
        answer[i,:] =_answer[i]
        gt[i,:] =_gt[i]
    end

    @test precision(answer=answer[1:1,:], ground_truth=gt[1:1,:]) == 1/2
    @test precision(answer=answer, ground_truth=gt) == 7/13
    @test recall(answer=answer[1:1,:], ground_truth=gt[1:1,:]) == 1
    @test recall(answer=answer, ground_truth=gt) == 7/9
end

# @test get_ground_truth(["1,2,3,4,5,6,7", "2,3,4,5,6,7,8"]) == [2 3 4 5 6 7 0; 3 4 5 6 7 8 0]


@testset "top k tests" begin
    _answer = [[[1,0,0,0,0],
               [0,1,0,0,0],
               [0,0,1,0,0],
               [0,0,0,1,0]],
               [[0,0,1,0,1],
               [0,1,0,1,0],
               [0,0,0,0,1],
               [1,1,1,1,1]] ]

    _gt = [[1,0,0,0,0],
           [0,0,0,0,1]]
    answer = Array{Int64}(undef, (2, 4, 5))
    gt = Array{Int64}(undef, (2 ,5))

    for i in 1:size(answer, 1)
        for j in 1:size(answer,2)
            answer[i,j,:] =_answer[i][j]
        end
        gt[i,:] =_gt[i]
    end
    
    @test in_top_k(answer[1,:,:], gt[1,:], 1) == true
    @test in_top_k(answer[2,:,:], gt[2,:], 1) == false
    @test in_top_k(answer[2,:,:], gt[2,:], 3) == true
    @test top_k_accuracy(answer, gt, 1) == 0.5
    @test top_k_accuracy(answer, gt, 3) == 1.
end
println("End of util.jl test")
end # module