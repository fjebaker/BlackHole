
""" Coordinate Conversion Functions """

function spher2xyz(r, θ, ϕ, dr, dθ, dϕ, d::AccretionDisk)
    """ 
    Also performs relevant transformation for plane of disk.
    Takes Jacobian into account.

    Symbolic expressions calculated with SageMaths.
    """

    # trigonometric values 
    ca = cos(d.α)
    cb = cos(d.β)
    sa = sin(d.α)
    sb = sin(d.β)

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

function xyz2spher(x, y, z, dx, dy, dz)
    r = sqrt(x^2 + y^2 + z^2)
    dr = (x * dx + y * dy * z * dz) / r
    dθ = ((x * z * dx + y * z * dy) / r^2 / sqrt(x^2 + y^2) - sqrt(x^2 + y^2) * dz) / r^2
    dϕ = -y / (x^2 + y^2) * dx + x / (x^2 + y^2) * dy

    (dr, dθ, dϕ)
end
