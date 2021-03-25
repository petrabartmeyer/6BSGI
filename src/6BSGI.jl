
include("read_data.jl")
include("models.jl")


hydrosys_folder = "p1"

hydro_data = read_data(hydrosys_folder)

model = Model(Gurobi.Optimizer)
create_model(model,hydro_data)

