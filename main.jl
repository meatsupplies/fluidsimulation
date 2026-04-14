using Plots

@kwdef mutable struct Particle
    const position::Tuple{Int, Int}
    height::Float64 = 0.0
    velocity::Float64 = 1.0
    temperature::Float64 = 50.0
    viscosity::Float64 = 1.0
end

@kwdef mutable struct Liquid
    size::Float64
    resolution::Int
    Δx::Float64
    sea_level::Float64
    particles::Array{Particle, 2}

    function Liquid(size::Float64, resolution::Int, sea_level::Float64 = 0.0)::Liquid
        Δx = size / resolution
        particles = Array{Particle, 2}(undef, resolution, resolution)

        for i in 1:resolution
            for j in 1:resolution
                particles[i, j] = Particle((i, j), sea_level)
            end
        end

        return Liquid(size, resolution, Δx, sea_level, particles)
    end

    function get_cardinal_neighbors(liquid::Liquid, i::Int, j::Int)::Vector{Particle}
        cardinal_neighbors = Particle[]
        if i > 1
            push!(cardinal_neighbors, liquid.particles[i-1, j])
        end
        if i < resolution
            push!(cardinal_neighbors, liquid.particles[i+1, j])
        end
        if j > 1
            push!(cardinal_neighbors, liquid.particles[i, j-1])
        end
        if j < resolution
            push!(cardinal_neighbors, liquid.particles[i, j+1])
        end
        return cardinal_neighbors
    end

    function get_diagonal_neighbors(liquid::Liquid, i::Int, j::Int)::Vector{Particle}
        diagonal_neighbors = Particle[]
        if i > 1 && j > 1
            push!(diagonal_neighbors, liquid.particles[i-1, j-1])
        end
        if i > 1 && j < resolution
            push!(diagonal_neighbors, liquid.particles[i-1, j+1])
        end
        if i < resolution && j > 1
            push!(diagonal_neighbors, liquid.particles[i+1, j-1])
        end
        if i < resolution && j < resolution
            push!(diagonal_neighbors, liquid.particles[i+1, j+1])
        end
        return diagonal_neighbors
    end

    function acceleration(viscosity::Float64, Δh::Float64, Δx::Float64)::Float64
        return -viscosity * Δh * (1 - (Δx / √(Δx^2 + Δh^2)))
    end

    #TODO Using explicit euler, change to something better
    function update_positions!(liquid::Liquid, dt::Float64 = 0.1)::Nothing
        for i in 1:resolution
             for j in 1:resolution
                particle::Particle = liquid.particles[i, j]
                particle.height += particle.velocity * dt

                cardinal_neighbors::Vector{Particle} = liquid.get_cardinal_neighbors(liquid, i, j)
                diagonal_neighbors::Vector{Particle} = liquid.get_diagonal_neighbors(liquid, i, j)
                acceleration::Float64 = 0
                for neighbor in cardinal_neighbors
                    Δh = particle.height - neighbor.height
                    acceleration += acceleration(neighbor.viscosity, Δh, liquid.Δx)
                end
                for neighbor in diagonal_neighbors
                    Δh = particle.height - neighbor.height
                    acceleration += 0.5*acceleration(neighbor.viscosity, Δh, liquid.Δx)
                end
                acceleration *= (1 - particle.viscosity)
                particle.velocity += acceleration * dt
        end
    end
end
# Initialize water
water = Liquid(size = 10.0, resolution = 10, sea_level = 0.0)



# Plotting 
x = [water.particles[i].position[1] for i in 1:resolution^2]
y = [water.particles[i].position[2] for i in 1:resolution^2]
z = [water.particles[i].height for i in 1:resolution^2]
v = [water.particles[i].velocity for i in 1:resolution^2]
plt = plot(xlim=(0,10), ylim=(0,10), zlim=(-5,5))
scatter!(x, y, z, color=:blue)
quiver!(x, y, z, quiver = (0,0,v))



