# Black Hole

![Example Render](examples/thindisk.png)


Black Hole rendering work written in Julia.

Project is work in progress. I used SageMath to calculate the symbolic expressions for the Christoffel symbols of the Eddington-Finkelstein metric, and then hard-coded the non-zero components as Julia functions, which are then used for a geodesic ray-tracing approach to rendering accretion disks around a black hole.

### Currently working on:
- GPU accelerated plots
- Finishing implementation of optically thin disks
- Relativistic effects

## Example Usage
```julia
s = EddingtonFinkelstein(2.0)

geodesics = calcgeodesics(s)

disk = GeometricDisk()

image = renderdisk(g, geodesics, width=720, height=480)

# save image
using Images
save("render.png", image)
```

For an optically thin disk, may use something like
```julia
disk = GaussianThinDisk(
    α=π/50, 
    β=0.0, 
    rinner=12,
    router=44,
    s=s,
    σ=4
)
```
as a drop in replacement for `disk`. Note, there is currently no GPU support for the optically thin disks.