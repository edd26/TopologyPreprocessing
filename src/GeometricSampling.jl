using LinearAlgebra

# integrate sinh^(n)(ax) from 0 to r
function integrate_sinh(n; r=1.0, a=1.0)
    if n==0
        return r
    elseif n==1
        return (cosh(a*r)-1)/a
    else
        return (sinh(a*r)^(n-1))*cosh(a*r)/(a*n) - (n-1)/n * integrate_sinh(n-2,r=r, a=a)
    end               
end    
    

function hyp_radial_density(r, d; curvature=-1.0, radius=1.0)
#    k = 1/(curvature^2)
    k = 1.0
    hyp_r = (sinh(r/k))^(d-1)/integrate_sinh(d-1,r=radius, a=radius/k)
    return hyp_r
end    


function euc_radial_density(r, d; radius=1.0)
    return d*(r^(d-1))/radius^d
end

# rejection sampling n=numofpts points from density dens, where the argument lies between 0 and maxval

function rejection_sampling(dens::Function, maxval,numofpts=1)
    max_val = dens(maxval)
    iter = 1
    rands = Array{Float64,1}(undef, numofpts)
    while iter <= numofpts
        x = rand()*maxval
        val = dens(x)
        u = rand()        
        if u*max_val < val
            rands[iter] = x
            iter+=1
        end
    end
    return rands    
end


function sample_hyp_rad(d, numofpts=1; curvature=-1.0, radius=1.0)
    rands = rejection_sampling(x->hyp_radial_density(x,d,curvature=curvature, radius=radius), radius,numofpts)
    
    # ...radius within the Poincare ball
    euc_rands = map(x->tanh(x/2.0),rands)
    return euc_rands
end

function sample_euc_rad(d, numofpts=1; radius=1.0)
    rands = rejection_sampling(x->euc_radial_density(x,d,radius=radius), radius, numofpts)
    return rands
end

function sample_sph(d, numofpts=1; curvature=1.0)
    rands = []
    i=0
    while i<=numofpts
        vec = randn(d+1)
        if vec[d+1]>0    
            push!(rands, normalize(vec))
            i+=1
        end
    end
    return rands
end

function sample_sphere(d, numofpts=1)
    vecs = randn(d, numofpts)
    rands = []
    for i=1:numofpts
        push!(rands, normalize(vecs[:,i]))
    end
    return rands
end

function sample_hyp(d, numofpts=1; radius=1.0, curvature=-1)
    sphere_pts = sample_sphere(d,numofpts)
    radii = sample_hyp_rad(d, numofpts, radius=radius, curvature=curvature)
    ball_pts = [radii[i]*sphere_pts[i] for i=1:numofpts]
    return ball_pts    
end

function sample_euc(d, numofpts=1; radius=1.0)
    sphere_pts = sample_sphere(d,numofpts)
    radii = sample_euc_rad(d, numofpts, radius=radius)
    ball_pts = [radii[i]*sphere_pts[i] for i=1:numofpts]
    return ball_pts    
end


function sample_ball(d, numofpts=1; radius=1.0, curvature=0.0)
    sphere_pts = sample_sphere(d,numofpts)
    if curvature < 0
        radii = sample_hyp_rad(d, numofpts, radius=radius, curvature=curvature)
        ball_pts = [radii[i]*sphere_pts[i] for i=1:numofpts]
    elseif curvature == 0.0
        radii = sample_euc_rad(d, numofpts, radius=radius)
        ball_pts = [radii[i]*sphere_pts[i] for i=1:numofpts]
    elseif curvature > 0
        ball_pts = sample_sph(d, numofpts, curvature=curvature)      
    end
end

function hyp_distance(pts; curvature=-1.0)
    distances = zeros(length(pts), length(pts))
    for i=1:length(pts)
        for j=1:i-1
            nx = 1-(norm(pts[i]))^2
            ny = 1-(norm(pts[j]))^2
            delta = 2 * norm(pts[i]-pts[j])^2/(nx*ny)
            distances[i,j] = acosh(1+delta)
            distances[j,i] = distances[i,j]
        end
    end
    return distances
end

function euc_distance(pts)
    distances = zeros(length(pts), length(pts))
    for i=1:length(pts)
        for j=1:i-1
            distances[i,j] = norm(pts[i]-pts[j])
            distances[j,i] = distances[i,j]
        end
    end
    return distances
end

function sph_distance(pts; curvature=1.0)
    distances = zeros(length(pts), length(pts))
    for i=1:length(pts)
        for j=1:i-1
            distances[i,j] = acos(dot(pts[i],pts[j]))
            distances[j,i] = distances[i,j]
        end
    end
    return distances
end

function distance_matrix(pts; curvature=0.0)
    if curvature < 0
        return hyp_distance(pts, curvature=curvature)
    elseif curvature == 0
        return euc_distance(pts)
    elseif curvature > 0
        return sph_distance(pts, curvature=curvature)
    end
end

function to_density(matr)
    dens_matr = zeros(size(matr,1), size(matr,2))
    n = size(matr)[1]
    all_entries = sort(setdiff(unique(matr), 0.0))
    total = binomial(n,2)
    for i=1:n
        for j=i+1:n
            dens_matr[i,j] = (findfirst(x->x==matr[i,j], all_entries))/total
            dens_matr[j,i] = dens_matr[i,j]
        end
    end
    return dens_matr           
end

function to_density!(matr)
    matr = to_density(matr)
    return matr           
end
 
