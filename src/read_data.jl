using PyCall

"""
    Install Conda package `pandas` using
    ```
    Conda.add("pandas")
    ```
    in order to get read_data.py working.
"""

function read_data(hydrosys_folder::String)
    read_datapydir = joinpath(@__DIR__,"../data")
    p1_dir = joinpath(read_datapydir,hydrosys_folder)

    pushfirst!(PyVector(pyimport("sys")."path"), read_datapydir)

    hydrosys_data = pyimport("read_data")[:load_data](p1_dir)
end