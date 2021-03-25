using PyCall

"""
    Install Conda package `pandas` using
    ```
    Conda.add("pandas")
    ```
    in order to get read_data.py working.
"""

function read_data(hydrosys_folder::String,instance_folder::String)
    read_datapydir = joinpath(@__DIR__,"../data")
    problem_dir = joinpath(read_datapydir,hydrosys_folder)
    instance_dir = joinpath(problem_dir,instance_folder)

    pushfirst!(PyVector(pyimport("sys")."path"), read_datapydir)

    hydrosys_data = pyimport("read_data")[:load_data](problem_dir)
    hydrosys_instance = pyimport("read_data")[:load_instance_df](problem_dir,instance_dir,hydrosys_data)
    return hydrosys_data, hydrosys_instance
end