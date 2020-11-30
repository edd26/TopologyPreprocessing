# Module taken from:  https://github.com/alexyarosh/hyperbolic

using Plots, Eirene, Ripserer, Statistics

# compute the persistent betti numbers of the Vietoris-Rips complex given by the distance matrix `matr`
# if mintime, maxtime and numofsteps are specified -- returns a `numofsteps x maxdim` array
# if either of the keyword arguments is not specified or is set to Inf, returns `maxdim` arrays for method=:eirene, or error for :ripser

# function bettis(matr, maxdim; mintime=-Inf, maxtime=Inf, numofsteps=Inf, method=:ripser)
#     if (method == :ripser) || (method == :Ripser)
#         if VERSION < v"0.7.0"
#             error("Ripser requires at least Julia v 0.7.0")
#         end
#         return bettis_ripser(matr, maxdim, mintime=mintime, maxtime=maxtime, numofsteps=numofsteps)
#     elseif (method == :eirene) || (method == :Eirene)
#         return bettis_eirene(matr, maxdim, mintime=mintime, maxtime=maxtime, numofsteps=numofsteps)
#     else
#         error("Method $(method) is not supported. Supported methods are: method=:eirene, method=:ripser")
#     end
# end

function bettis_ripser(matr, maxdim; mintime=-Inf, maxtime=Inf, numofsteps=Inf)
    if (mintime == -Inf) || (maxtime == -Inf) || (numofsteps == -Inf)
        error("To use Ripser, specify parameters mintime, maxtime, numofsteps")
    end
    r = ripser(matr,  dim_max = maxdim, threshold = maxtime)

    int_length = maxtime-mintime
    step_length= int_length/numofsteps
    betts = zeros(numofsteps, maxdim)

    for dim=1:maxdim
        ints = r[dim+1]
        for intl in ints
            st = Int(ceil((intl[1]-mintime)/step_length))
            if intl[2] == Inf
                fin = numofsteps
            else
                fin = Int(ceil((intl[2]-mintime)/step_length))
            end
            betts[st:fin, dim] = map(x->x+1, betts[st:fin, dim])
         end
    end
    return betts
end
#
# # Original function returns 2 different types of betti curves. If no default
# # value parameters is given, it returns vector of matrices. If num of steps is
# # given, then it return matrix maxdim x numsteps.
# function bettis_eirene(matr, maxdim; mintime=-Inf, maxtime=Inf, numofsteps=Inf)
#     c = eirene(matr, minrad = mintime, maxrad= maxtime, numrad= numofsteps, maxdim=maxdim)
#
#     int_length = maxtime-mintime
#     step_length= int_length/numofsteps
#
#     if (mintime == -Inf) || (maxtime == Inf) || (numofsteps == Inf)
#         # return [betticurve(c, dim=maxdim) for d=1:maxdim]
#         return hcat([betticurve(c, dim=d)[:,2] for d=1:maxdim]...)
#     end
#
#     betts = zeros(numofsteps, maxdim)
#     # For every dimension compute betti curve
#     for dim=1:maxdim
#         bet = betticurve(c, dim=dim)
#
#         #for every element in betti curve return betti value if index is positive
#         for i=1:size(bet,1)
#             b = bet[i,:]
#             ind = Int(ceil((b[1]-mintime)/step_length))
#             if ind > 0
#                 betts[ind,dim]=b[2]
#             else
#                 betts[1,dim]=b[2]
#             end
#         end
#     end
#     return betts
# end

# average betti numbers over arrs
# assuming arrs is an array of arrays, where each arrs[j] is the same size

function average_bettis(arrs; maxdim=-1)
    if size(arrs,2) > 1
        return arrs
    end
    md = maxdim
    if maxdim == -1
        md = size(arrs[1],2)
    end
    numofints = size(arrs[1],1)
    av_bet = zeros(numofints,md)

    for i=1:numofints
        for d=1:md
            av_bet[i,d] = mean([arrs[j][i,d] for j=1:length(arrs)])
        end
    end
    return av_bet
end

# compute standard deviation of betti numbers over arrays in arrs
# assuming arrs is an array of arrays, where each arrs[j] is the same size

function std_bettis(arrs; maxdim=-1)
    md = maxdim
    if maxdim == -1
        md = size(arrs[1],2)
    end
    numofints = size(arrs[1],1)
    std_bet = zeros(numofints,md)

    if size(arrs,2) > 1
        return std_bet
    end

    for i=1:numofints
        for d=1:md
            std_bet[i,d] = std([arrs[j][i,d] for j=1:length(arrs)])
        end
    end
    return std_bet
end

# plot average curves at values `xval`, with averages given by `means` and standard deviations given by `std`

function plot_averages(xvals, means, stds; ribbon=true, label="", linestyle=:solid, color=:auto)
    if ribbon
         return plot(xvals, means, ribbon=stds,fillalpha=.3, labels=label, linestyle=linestyle, color=color)
    else
        return plot(xvals, means, labels=label, linestyle=linestyle, c=color)
    end
end

function plot_averages!(xvals, means, stds; ribbon=true, label="", linestyle=:solid, color=:auto)
    if ribbon
         return plot!(xvals, means, ribbon=stds,fillalpha=.3, labels=label, linestyle=linestyle, color=color)
    else
        return plot!(xvals, means, labels=label, linestyle=linestyle, c=color)
    end
end

function load_bettis(filename)
    dict = load(filename)
    for (matr_name, matr) in dict
        return matr
    end
end

# plot average curves at values `xval`, given that the bettis numbers are saved in `file`

function plot_averages(xvals, file::String; dim=1, ribbon=true, label="", linestyle=:solid, color=:auto)
    matr = load_bettis(file)

    av = average_bettis(matr)[:,dim]
    if ribbon
         st = std_bettis(matr)[:,dim]
         return plot(xvals, av, ribbon=st,fillalpha=.3, labels=label, linestyle=linestyle, c=color)
    else
        return plot(xvals, av, labels=label, linestyle=linestyle, c=color)
    end
end


function plot_averages!(xvals, file::String; dim=1, ribbon=true, label="", linestyle=:solid, color=:auto)
    matr = load_bettis(file)

    av = average_bettis(matr)[:,dim]
    if ribbon
         st = std_bettis(matr)[:,dim]
         return plot!(xvals, av, ribbon=st,fillalpha=.3, labels=label, linestyle=linestyle, c=color)
    else
        return plot!(xvals, av, labels=label, linestyle=linestyle, c=color)
    end
end
