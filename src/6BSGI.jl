
include("read_data.jl")
include("models.jl")


hydrosys_folder = "p1"

data = read_data(hydrosys_folder)

usinas = data["usinas"]
# model = Model(Gurobi.Optimizer)
#model = Model(Ipopt.Optimizer)

# create_model(model)

