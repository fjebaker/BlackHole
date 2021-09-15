
""" Utility functions """


@doc raw"""

    dspher2xyz(r, θ, ϕ, dr, dθ, dϕ, disk::AccretionDisk)

Transform differential coordinates from spherical to cartesian, at a given point `r, θ, ϕ`,
also performing relevant transformation for plane of disk.

Returns `(dx, dy, dz)`.

Takes Jacobian into account. Symbolic expressions originally calculated with SageMaths.

Complement to `dxyz2spher`.
"""
function dspher2xyz(r, θ, ϕ, dr, dθ, dϕ, disk::AccretionDisk)::Tuple
    # trigonometric values
    ca = cos(disk.α)
    cb = cos(disk.β)
    sa = sin(disk.α)
    sb = sin(disk.β)

    ct = cos(θ)
    cp = cos(ϕ)
    st = sin(θ)
    sp = sin(ϕ)

    # component calculations
    dx = (
        (-cb * ct * sa - (sa * sb * sp - ca * cp) * st) * dr +
        (r * cb * sa * st - (sa * sb * sp - ca * cp) * r * ct) * dθ +
        (-(cp * sa * sb + ca * sp) * r * st) * dϕ
    )

    dy = (
        (cb * sp * st - sb * ct) * dr +
        (r * ct * cb * sp + r * sb * st) * dθ +
        (r * cp * cb * st) * dϕ
    )

    dz = (
        (ca * cb * ct + (ca * sb * sp + cp * sa) * st) * dr +
        (-r * ca * cb * st + (ca * sb * sp + cp * sa) * r * ct) * dθ +
        ((ca * cp * sb - sa * sp) * r * st) * dϕ
    )

    (dx, dy, dz)
end

@doc raw"""

    dxyz2spher(x, y, z, dx, dy, dz)

Transform differential coordinates from cartesian to spherical, at a given point `x, y, z`.

Returns `(dr, dθ, dϕ)`.

Takes Jacobian into account. Symbolic expressions originally calculated with SageMaths.

Complement to `dspher2xyz`.
"""
function dxyz2spher(x, y, z, dx, dy, dz)
    r = sqrt(x^2 + y^2 + z^2)
    dr = (x * dx + y * dy * z * dz) / r
    dθ = ((x * z * dx + y * z * dy) / r^2 / sqrt(x^2 + y^2) - sqrt(x^2 + y^2) * dz) / r^2
    dϕ = -y / (x^2 + y^2) * dx + x / (x^2 + y^2) * dy

    (dr, dθ, dϕ)
end


@doc raw"""

    spher2xyz(t, r, θ, ϕ)

Transform coordinates of a space time point in spherical to cartesian coordinates.

Returns (t, x, y, z).
"""
function spher2xyz(t, r, θ, ϕ)
    (t, r * sin(θ) * cos(ϕ), r * sin(θ) * sin(ϕ), r * cos(θ))
end

# tupled argument so it can accept arrays as input
spher2xyz((t, r, θ, ϕ, _...)) = spher2xyz(t, r, θ, ϕ)


@doc raw"""

    inclination_transform(x, y, z, α, β)

Transform a point on x,y,z to the inclination given by a rotation of α around x, and β around z.

Returns `(x, y, z)` in the rotated plane.
"""
function inclination_transform(x, y, z, α, β)
    y, z = (y * cos(β) - z * sin(β), z * cos(β) + y * sin(β))
    (x * cos(α) - z * sin(α), y, z * cos(α) + x * sin(α))
end


@doc raw"""

    truncator(px)::UInt8

Clamps value of `px` between 0 and 255.
"""
function truncator(px)::UInt8
    convert(UInt8, trunc(min(max(0, px), 255)))
end
