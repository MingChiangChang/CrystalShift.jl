module Testnode
using CrystalShift: Node, add_child!, get_nodes_at_level
using CrystalShift: is_immidiate_child, is_child, get_level
using CrystalShift: get_phase_ids, check_same_phase, remove_child!,get_node_with_id
using CrystalShift: CrystalPhase, get_param_nums
using Test

# CrystalPhase object creation
path = "../data/"
phase_path = path * "tree_test_sticks.csv"
f = open(phase_path, "r")

if Sys.iswindows()
    s = split(read(f, String), "#\r\n") # Windows: #\r\n ...
else
    s = split(read(f, String), "#\n")
end

if s[end] == ""
    pop!(s)
end

cs = Vector{CrystalPhase}(undef, size(s))
@. cs = CrystalPhase(String(s))
println("$(size(cs, 1)) phase objects created!")
x = LinRange(8, 45, 512)
y = cs[1].(x)+cs[2].(x)
y /= max(y...)
# Creating objects for testing
root = Node()
node1 = Node(cs[1], 2)
node2 = Node([cs[1], cs[2]], 3)
node3 = Node([cs[1], cs[2]], 4)
node4 = Node([cs[1], cs[2], cs[3]], 5)
show(node1)
testnode = Node(node2, [CrystalPhase(cs[1], rand(get_param_nums(cs[1]))),
                        CrystalPhase(cs[2], rand(get_param_nums(cs[2])))],
                        x, y)
testnode = Node(testnode.phase_model, [node1, node2], testnode.id, x, y)

add_child!(root, node1)
add_child!(node1, node2)
fake_tree = [root, node1, node2, node3, node4]

@testset "Basic Node properties" begin
    @test is_immidiate_child(node1, node2)
    @test node2 == node3
    @test node1 != node3
    @test is_child(node1, node2)
    @test is_child(node1, node4) # two level down
    @test root.child_node == [node1]
    @test root.child_node[1].child_node == [node2]
    @test is_immidiate_child(node1, node2)
    @test is_immidiate_child(node2, node4)
    @test !is_immidiate_child(node1, node4)
    @test get_level(node4) == 3
    @test get_phase_ids(node4) == [0, 1, 2]
    @test get_nodes_at_level(fake_tree, 3) == [node4]
    @test check_same_phase(node2, [cs[1], cs[2]])
    @test check_same_phase(node2.phase_model, [cs[1], cs[2]])
    @test check_same_phase([cs[1], cs[2]], node2.phase_model)
    remove_child!(testnode, node1)
    @test testnode.child_node == [node2]
    @test get_node_with_id(fake_tree, 0)[1] == node1
    @test get_node_with_id(fake_tree, [0, 1]) == [node1, node2, node3, node4]
end
println("End of node.jl test")
end # module