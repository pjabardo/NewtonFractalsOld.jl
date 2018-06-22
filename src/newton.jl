

using Polynomials
using Images
using ImageView

function newton_iteration(x0::T, p, dp, errmax=10*eps(), itmax=1000) where T

    imax = 0
    x1 = zero(T)
    
    for i = 1:itmax
        y0 = p(x0)
        dy0 = dp(x0)
        dx = -y0 / dy0

        x1 = x0 + dx
        x0 = x1

        if abs(dx) < errmax
            return x1, i, true
        end
    end

    return x1, itmax, false
    

end

function setupnewton(p, n=512, fator=0.3333) 

    z = roots(p)
    nz = length(z)
    zr = real.(z)
    zi = imag.(z)

    xmax = maximum(zr)
    xmin = minimum(zr)
    Lx = xmax-xmin

    ymax = maximum(zi)
    ymin = minimum(zi)
    Ly = ymax - ymin

    Lmax = max(Lx, Ly)
    Lmin = min(Lx, Ly)
    Lmax1 = Lmax * fator
    Lmin1 = max(Lmin*fator, Lmax1)
    nmax = n
    nmin = ceil(Int, n * (Lmin + 2*Lmin1) / (Lmax + 2*Lmax1))
    
    if Lx > Ly
        xx = linspace(xmin-Lmax1, xmax+Lmax1, nmax)
        yy = linspace(ymin-Lmin1, ymax+Lmin1, nmin)
    else
        yy = linspace(ymin-Lmax1, ymax+Lmax1, nmax)
        xx = linspace(xmin-Lmin1, xmax+Lmin1, nmin)
    end
    

    dp = polyder(p)

    
    return z, p, dp, xx, yy
    
  
end

function newtonfractal(p, n=512; fmin=0.2, imax=20)

    z, p, dp, xx, yy = setupnewton(p, n)
    frac = newtonfracloop(z, p, dp, xx, yy)
    
    img = newton2img(length(z), frac, fmin=fmin, imax=imax)
    imshow(img)
end

    

function newtonfracloop(z, p, dp, xx, yy, errmult=20, itmax=1000) 
    nx = length(xx)
    ny = length(yy)
    nz = length(z)
    frac = Matrix{Tuple{Int,Int}}(ny, nx)

    for ix = 1:nx
        x = xx[ix]
        for iy = 1:ny
            y = yy[iy]
            errmax = errmult*eps(hypot(x,y))
            r, niter, converged = newton_iteration(x+y*im, p, dp, errmax)
            rmin = Inf
            imin = -1

            for k = 1:nz
                ee = abs(r-z[k])
                if ee < rmin
                    rmin = ee
                    imin = k
                end
            end
            frac[iy, ix] = (imin, niter)
        end
    end

    return frac
            
end

#=function cfdcolors(x, xmin=0.0, ymin=0.0;
                   clamp=TRUE, colmin=RGB{Float64}(0,0,0),
                   colmax=RGB{Float64)(1,1,1))
    #ξ = (x - xmin) / (xmax - xmin)
    if clamp
        if ξ < 0
            return RGB{Float64}(0,0,1)
        elseif ξ > 1.0
            return RGB{Float64}(1,0,0)
        end
    else
        if ξ < 0
            return colmin
        elseif ξ > 0
            return colmax
        end
    end

    if ξ <= 0.25
        return RGB{Float64}(0.0, 4*ξ, 1.0)
    elseif ξ > 0.25 && ξ < 0.5
        return RGB{Float64}(0.0, 1.0, 2.0 - 4*ξ)
    else
        return RGB{Float64}(0.0, 4.0-4*ξ, 1.0)
    end

=#
function cfdcolors(ξ)


    if ξ <= 0.25
        return RGB{Float64}(0.0, 4*ξ, 1.0)
    elseif ξ > 0.25 && ξ < 0.5
        return RGB{Float64}(0.0, 1.0, 2.0 - 4*ξ)
    elseif ξ >= 0.5 && ξ < 0.75
        return RGB{Float64}(-2.0 + 4*ξ, 1.0, 0.0)
    else
        return RGB{Float64}(1.0, 4.0-4*ξ, 0.0)
    end
    
    
end

                         
function newton2img(nz, frac; fmin=0.3, imax=20)
    basecol = cfdcolors.(linspace(0,1,nz))

    ny, nx = size(frac)

    img = zeros(RGB{Float64},ny,nx)
    for ix = 1:nx
        for iy = 1:ny
            icol, niter = frac[iy,ix]
            niter = min(niter,imax)
            r = red(basecol[icol])
            g = green(basecol[icol])
            b = blue(basecol[icol])
            fator = clamp((fmin - 1.0) * (niter - 4) / (imax - 4) + 1, fmin, 1.0)
            img[iy,ix] = RGB{Float64}(fator*r, fator*g, fator*b)
        end
    end
    return img
    
end


