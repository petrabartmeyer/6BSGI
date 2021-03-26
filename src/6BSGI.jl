using JuMP, Ipopt

include("read_data.jl")
include("models.jl")


hydrosys_folder = "p1"
hydrosys_instance = "i2"

hydro_data,hydro_instance = read_data(hydrosys_folder,hydrosys_instance)

model = Model(optimizer_with_attributes(
    Ipopt.Optimizer,
    "constr_viol_tol" => 1.0e-2,
    "acceptable_tol"  => 1.0e-5,
    "print_level"     => 4
))

create_model(model,hydro_data,hydro_instance)
optimize!(model)

println(value.(model[:pg]))
println(value(model[:mercado]))
