using Pkg

# Get a list of all installed packages
installed_packages = keys(Pkg.project().dependencies)
installed_packages_vector = collect(installed_packages)
installed_packages_vector = push!(installed_packages_vector, "GLM")
println(installed_packages_vector)

# List of packages to exclude
packages_to_remove = ["DynamicGrids", "Dispersal", "CUDA", "EcologicalNetworksDynamics", "OrdinaryDiffEq", "Clustering"]

# Remove specified packages from the list
updated_packages = filter(x -> !(x in packages_to_remove), installed_packages_vector)
println(updated_packages)

for pkg_name in installed_packages_vector
    Pkg.rm(pkg_name)
end

for pkg_name in updated_packages
    Pkg.add(pkg_name)
end

