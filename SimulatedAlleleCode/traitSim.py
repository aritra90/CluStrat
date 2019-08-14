import numpy as np
import math
from sklearn.cluster import KMeans

################################################################################
################################################################################
################################################################################

def simulate(normX,S,v,m,n,d):
    # Define the amount of variances for each effect
    v_gen = v[0]/100.0;
    v_env = v[1]/100.0;
    v_noise = v[2]/100.0;
    # set number of causal SNPs
    num_causal = 10;

    # get the SNP effects for the user defined number of SNPs
    effects = np.random.normal(0,0.5,num_causal)
    # rest of the SNPs will be set to 0 for the effects
    # Just appends zeros to the remainder of the effects array
    effects = np.append(effects,np.zeros((m-num_causal,1)))

    # genetic effect (calculating just for causal variants as rest are 0
    for i in range(0,num_causal):
       normX[i,:]  = normX[i,:]*effects[i]
    # sum of each column of normX
    gen_eff = normX.sum(axis = 0)

    # rescale genetic effects by the variance factor
    fact = math.sqrt(v_gen)/(np.std(gen_eff));
    # element-wise rescaling of the genetic effects
    gen_eff = fact*gen_eff;

    noise = np.zeros((n,1));
    # k-means clustering of the columns of the population admixture matrix
    # s.t. each individual falls in one of 'd' clusters
    theclustering = KMeans(n_clusters=d).fit(S.T)
    # print(theclustering.labels_)
    lambda_k = theclustering.labels_ + 1
    # some gamma distribution for the noise
    gamvec = 1/np.random.gamma(3,1,size=(d,1))

    for i in range(0,n):
        noise[i] = np.random.normal(0,gamvec[lambda_k[i]-1])

    # rescale lambda_k (environmental effects)
    mc_lambda_k = np.square((lambda_k - np.mean(lambda_k)));
    fact = math.sqrt(v_env)/math.sqrt(np.sum(mc_lambda_k)/(n-1));
    # element-wise rescaling of the environmental effects
    lambda_k = fact*lambda_k;

    # rescale noise term
    mc_noise = np.square((noise - np.mean(noise)))
    fact = math.sqrt(v_noise)/math.sqrt(np.sum(mc_noise)/(n-1));
    # element-wise rescaling of the noise effects
    noise = fact*noise;

    # Get the simulated quantitative trait for each individual
    # apostrophe --> transpose

    gen_eff = gen_eff.reshape((n,1))
    lambda_k = lambda_k.reshape((n,1))
    traits = gen_eff + lambda_k + noise
    print(' ')

    # Binary traits
    bin_trait = traits-noise;
    # probability that the individual has case status (vs control)
    prob_case = np.power(10,bin_trait)/(1+np.power(10,bin_trait))
    # print(np.where(prob_case>0.5)[0])
    status = np.zeros((n,1));
    status[np.where(prob_case > 0.5)[0]] = 1;

    return traits, status


if __name__ == '__main__':
    print("HEY HEY HEY")
