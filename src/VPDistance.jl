# # ======================
# # Example usage
# using CSV
# using Plots
#
# file_name = "sts_for_VP_test.csv"
# csv_matrix = CSV.read(file_name)[:,2:end]
# almost_not_csv_matrix = Matrix(csv_matrix)
#
# file_name2 = "spikes.csv"
# spikes = load_csv_file_to_array(file_name2)
#
# file_name3 = "th_chunks.csv"
# th_chunks = load_csv_file_to_array(file_name3)
#
# sts = generate_spike_matrix(spikes; th_chunks=th_chunks)
# VPd = [get_selfdist(s; n_chan=32, cost=60., dt=0.01) for s in sts]
#
#
# plot_set = Any[]
#     for matrix in VPd
#         push!(plot_set, heatmap(matrix, color=:Reds_9, colorbar=false, yflip = true))
#     end
#
# plot(plot_set[1], plot_set[2], plot_set[3],
#     plot_set[4], plot_set[5], plot_set[6],
#     plot_set[7], plot_set[8], plot_set[9],
#     layout=(1,9), size=(9*1200,1100), legend=false, colorbar=false)



"""
    get_selfdist(st_inp; n_chan=32, cost=60., dt=0.01)

Method for computing pair-wise spike distances from a range of spike trains.
Function copied from Mikolaj SETCOmodel

Inputs:
    st_inp: [2 x N] array with spike times and indices of neurons.
    N - number of spikes generated, 1st row - index of neuron generating given spikes, 2nd row - spike time.
    n_chan - number of neurons (default: 32)
    cost - cost parameter for VP spike distance, in ms (default: 60 ms)
    dt - simulation timestep, in ms (default: 0.01 ms -> 100 kHz)
Output:
    pc - [n_chan x n_chan] matrix containing pairwise VP spikes distances for each pair of neurons.
"""
function get_selfdist(st_inp; n_chan=32, cost=60., dt=0.01)

    sts_new = Any[]
    for i in 1:n_chan
        push!(sts_new, st_inp[2,findall(x->x==i, st_inp[1,:])])
    end

    # sts = [st_inp[0,st_inp[1,:]==i] for i in 1:n_chan]
    pc = zeros(n_chan, n_chan)

    for i in 1:n_chan, j in 1:n_chan
        pc[i,j] = spkd(sts_new[i], sts_new[j], dt/(cost))
    end
    return pc
end


# TODO Add test with MATLAB code run for some spike train and compare with
# results from this
"""
    spkd(s1, s2, cost)

Fast implementation of victor-purpura spike distance (faster than neo & elephant python packages)
Direct Python port of http://www-users.med.cornell.edu/~jdvicto/pubalgor.htmlself.
The below code was tested against the original implementation and yielded exact results.
All credits go to the authors of the original code.

Code was translated from Frotran to Matlab, from Matlab to Python, from
Python to Julia. It was veryfied with MATLAB code.

Input:
    s1,2: pair of vectors of spike times
    cost: cost parameter for computing Victor-Purpura spike distance.
    (Note: the above need to have the same units!)
Output:
    d: VP spike distance.
"""
function spkd(s1, s2, cost)
    nspi=length(s1)
    nspj=length(s2)

    # Why not to use this?
    if cost==0
       return d=abs(nspi-nspj)
    elseif cost==Inf
       return d=nspi+nspj
   end

    scr=zeros(nspi+1, nspj+1)

    # initialize margins with cost of adding a spike
    scr[:,1]=0:nspi
    scr[1,:]=0:nspj

    for i = 2:nspi+1, j = 2:nspj+1
        component1 = scr[i-1,j]+1
        component2 = scr[i,j-1]+1
        component3 = scr[i-1,j-1]+cost*abs(s1[i-1]-s2[j-1])
        scr[i,j] = min(component1, component2, component3)
	end

    d=scr[end,end]

    return d
end


"""
	generate_spike_matrix(spikes; th_chunks=[[0,0]])

Generates matrix of the form [2xN] array with spike times and indices of
neurons. N - number of spikes generated, 1st row - index of neuron generating
given spikes, 2nd row - spike time.

Resulting matrix is time sorted. Resulting matrix may be splint into fragments
by setting 'th_chunks' to a list of fragments in a way such that
'th_chunks[k,1]' is a starting index and 'th_chunks[k,2]' is an ending index of
k'th fragment.
"""
function generate_spike_matrix(spikes; th_chunks=[[0,0]], val_range=20:52)

    if th_chunks == [[0,0]]
        th_chunks[1,1] = 1
        th_chunks[1,2] = size(spikes,1)
    end

    spike_train_simplified = Any[]
    for i = 1:size(th_chunks,1)
        #1 Get a chunk of spike trains of interest
        syllable_spikes = spikes[th_chunks[i, 1]:th_chunks[i, 2], val_range]
        # find all spikes
        all_spikes = findall(x->x==1,syllable_spikes)
        # convert to [2xN] matrix
        total_spikes = length(all_spikes)
        sorted_spikes = zeros(Int, 2,total_spikes)
        for k = 1:total_spikes
            sorted_spikes[1,k] = all_spikes[k][2]
            sorted_spikes[2,k] = all_spikes[k][1]
        end

        push!(spike_train_simplified, sorted_spikes[:,sortperm(sorted_spikes[2,:])])
    end

    return spike_train_simplified
end
