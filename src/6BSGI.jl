
include("read_data.jl")
include("models.jl")


hydrosys_folder = "p1"
hydrosys_instance = "i2"

hydro_data,hydro_instance = read_data(hydrosys_folder,hydrosys_instance)

model = Model(Gurobi.Optimizer)
create_model(model,hydro_data,hydro_instance)
optimize!(model)
