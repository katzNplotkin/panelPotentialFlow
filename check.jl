using FileIO

function getMesh(filename)
    mesh = load(filename)
    println(length(mesh))
    # println(Base.invokelatest(length, mesh))
    return mesh
end

m = getMesh("sphere.stl")
