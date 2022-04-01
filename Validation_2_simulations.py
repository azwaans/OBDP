import numpy as np
import ete3
import obdp
import nt
import os


def simulateOneDataset(dirname):
    os.mkdir(dirname)
    os.chdir(dirname)
    ##### SIMULATION OF PARAMETERS ######

    # tor prior: uniform on (tmin, tmax)
    tmin = 1
    tmax = 5
    tor = tmin + np.random.random()*(tmax-tmin)
    
   # div_mean = np.log( np.log(5) / tor )
   # div_sd = 0.2
   # lambMinusMu = np.random.lognormal(mean=div_mean, sigma=div_sd)

    # mu prior: exponential distribution with mean 1
    mu = np.random.exponential(1.)

    # lambda prior: knowing mu, it's mu + exponential distribution with mean 0.01 (i.e. rate 100)
    lambMinusMu = np.random.exponential(0.01)
    lamb = mu + lambMinusMu

    # psi and omega priors are exponential distributions with mean 0.2 (i.e. rate 5)
    psi = np.random.exponential(0.2)
    omega = np.random.exponential(0.2)

    # r is uniformly chosen on (0,1)
    r = np.random.random()

    # rho is uniformly distributed on (0.5, 1)
    rho = 0.8 + np.random.random()*0.2

    # mutation rate prior: exponential distribution with mean 0.05 (i.e. rate 20)
    alpha = np.random.exponential(0.05)

    # we record the simulated parameters
    params_obdp = lamb, mu, rho, psi, r, omega
    print(params_obdp)
    params_nt = [alpha]
    allParams = [lamb, mu, rho, psi, r, omega, alpha]
    allParamsNames = ['lamb', 'mu', 'rho', 'psi', 'r', 'omega', 'alpha']
    obdp.exportParams(allParams, allParamsNames)
    #print(params_obdp)


    ##### SIMULATION OF THE TREE ######

    # We simulate trees conditioned on having one leaf reaching the present
    t, obs, ttruth = obdp.simTOconditionedOnSurvival(params_obdp, tor)
    #t.show()
    print(len(obs[0]), len(obs[1]), len(obs[2]), len(obs[3]), len(obs[4]), len(obs[5]))

    # record the reconstructed tree
    t.write(format=1, outfile="tree_reconstructed.nw")
    # record the full tree (with extinct branches)
    ttruth.write(format=1, outfile="tree_full.nw")
    # record the list of occurrences
    obdp.exportOccurrences( obs )
    # and finally, the list of psi and rho sampled individuals
    obdp.exportTaxa( t )


    ##### SIMULATION OF SEQUENCES ######

    # number of nucleotides in the sequence
    m_nt = 1000

    # simulation of genetic sequences along the reconstructed tree, under a strict clock, with model JC69
    tree_seq = nt.simSeqAlongTree(params_nt, t, m_nt)

    # we export sequences for psi and rho sampled individuals (i.e. "epidemiology-style")
    nt.writeNexus( tree_seq, "extant-extinct" )

    os.chdir("..")
    return 'done'


# uncomment the following line to get a deterministic output
np.random.seed(1)
for i in range(0, 1000):
    dirname = "dataset"+str(i)
    simulateOneDataset(dirname)


