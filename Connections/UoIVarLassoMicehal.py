#!/usr/bin/env python
# coding: utf-8

# In[ ]:


# Misc
import time, os, sys, pdb

# Base
import numpy as np
import pandas as pd
import scipy.sparse as sparse
from scipy.sparse.linalg import inv
from sklearn.linear_model import Lasso, LogisticRegression #, LassoCV, LinearRegression
from functools import reduce


# Plot
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.colors import LinearSegmentedColormap



### Create function to pull top n most active neurons

# neuro_top_n() will take inputs:
# df: a pandas df where each column is a neuron and each row corresponds to a time step
# n: the desired number of "top n 'most active' neurons." Default is 100.
# method: the method by which "most active" is determined. Default is mean.
# It will return outputs:
# df_top_n: a pandas dataframe consisting of the top most active n neurons, based on the given method
# I plan to build this function out to use a variety of methods, such as pth percentile, but will start with just mean
def neuro_top_n(df, n = 100, method = "mean"):
    if method == "mean":
        # Find the mean value for each neuron
        df_means = np.mean(df, axis = 0)
        # Find the n neurons with highest mean values
        cols = df_means.nlargest(n)
        df_top_n = df.loc[:,cols.index]
        # This was the method used previously
        #df_top_n = df.iloc[:,cols.index]
        return df_top_n
    else:
        print("Please enter a valid method.")



### Create a function to split the entire trial into n second intervals, to be analyzed independently

# df_splitter will take inputs:
# df: a pandas df where each column is a neuron and each row corresponds to a time step
# Hz: the frequency of the time series contained in df, in Hz (observations per second). Default is 10 Hz.
# chunk_secs: the desired length (in seconds) of the new dataframes contained in the list. Default is 10 seconds.
# It will return outputs:
# df_list: a list consisting of contiguous pandas dataframes of a specified length, together composing the original dataframe
def df_splitter(df, Hz = 10, chunk_secs = 10):
    chunk_size = Hz * chunk_secs # sampling frequency in Hz times intended chunk size in seconds gives appropriate number of rows per chunk
    # compute the number of chunks 
    # we're floor dividing, so each resulting coefficient matrix is of the same size 
    # This means we can drop between 0 and (Hz*chunk size) - 1 rows from the original data
    num_chunks = len(df) // chunk_size
    # Create an empty list to hold all of the chunks of df
    df_list = []
    # loop over the dataframe to extract each chunk_secs second chunk into its own dataframe
    for i in range(num_chunks):
        # Compute index value to start chunk
        start = i * chunk_size
        # Compute index value to end chunk
        end = start + chunk_size
        # Take the ith consecutive 10 second chunk from the dataframe
        chunk_df = df.iloc[start:end]
        # Append this datafrane to our list
        df_list.append(chunk_df)
    # How many chunks did we get?
    # the f"<words {var} words>" form tells python to use f-string formatting, which is quick, and displays variables in {} in the output
    print(f"Created {len(df_list)} dataframes")
    return df_list



### Create a moving block bootstrap function

# This doesn't inherit the row index - I'm not sure if I prefer that or not

# block_length is the number of observations in a single block - not the number of seconds
def mbbs(data, n_blocks = 10, block_length = None):
    # number of datapoins
    n = len(data)
    
    # assign a default value for block_length
    # this ensures the resaulting sample will be as close to the size of the original sample as possible, 
    # while prioritizing all blocks being of the same length, and the result should not be longer than the original.
    if block_length is None:
        block_length = n // n_blocks
        
    # number of samples to take
    n_samples = n_blocks * block_length
    # initiate empty array to store bootstrap sample
    mbbs_sample = pd.DataFrame().reindex_like(data)
    
    # generate starting indices for each block
    block_indices = np.random.randint(low = 0, high = n - block_length, size = n_blocks)
    
    # Create the bootstrap sample
    for i in range(n_blocks):
        start_index = block_indices[i]
        end_index = start_index + block_length
        mbbs_sample.iloc[i * block_length : (i + 1) * block_length] = data.iloc[start_index:end_index]

    return mbbs_sample


# In[ ]:


### Create a function that will create the design matrrix for a VAR(D) model

def VAR_design(df, d = 1, n_obs = None, n_var = None, intercept = False):
    # Because we won't always feel like defining n_var and n_obs first
    if n_var is None:
        n_var = df.shape[1]
    if n_obs is None:
        n_obs = df.shape[0]
    #X = np.ones((n_obs-d, n_var*d))
    X = np.full([n_obs-d, n_var*d], np.nan)
    for i in range(d, n_obs):
        for j in range(d):
            X[i-d, n_var*j:n_var*(j+1)] = df.iloc[i-j-1, :]
    if intercept == True:
        X = np.insert(X, 0, np.ones((1, n_obs-d)), axis = 1)
    return X


# In[ ]:


def unique_generator(lst):
    seen = set()
    for item in lst:
        tupled = tuple(item)
        if tupled not in seen:
            seen.add(tupled)
            yield item


# In[ ]:


### Create a function the will combine the intersection and union steps into one UoI_VAR function

def UoI_VAR_LASSO(df, reg_path, B1 = 100, d = 1, intercept = False, method = 'OLS',
                  n_blocks = 10, block_length = None,
                  B2 = 100, gamma = 1):

    # Setting up parameters
    n_obs = df.shape[0]
    n_var = df.shape[1]
    I_M = np.identity(n_var) 
    K = len(reg_path)
    supp_list = [[] for i in range(K)]
    ### Intersection Step
    for i in range(B1):
        ## draw bootstrap sample
        bs = mbbs(df, n_blocks, block_length)
        ## Construct Y*, U* based on boostrap, according to equation (2.2) (Ruiz 2020)
        # Create design matrix for lag 1 that includes a column of 1s to explicitly solve for the intercept term
        
        X = VAR_design(bs, d, n_obs, n_var, intercept)
        
        # we aren't really vectorizing X, but we'll call it x_vec for consistency with y
        # we'll take the Kronecker of X and I_100, as per Ruiz (2020)
        X_vec = np.kron(X, I_M)
        
        # Convert X to a sparse matrix, to hopefully speed up computation time
        X_csr = sparse.csr_matrix(X_vec)
        
        # we need to lose the first d observations of each neuron to account for the lag order of d
        Y = bs.iloc[d:, :]
        # Y needs to be an ndarray to use reshape
        Y_ndarray = Y.to_numpy()
        # We're vectorizing Y via the column stacking operator
        # -1 in the reshape function tells python to infer the number of rows based on the specified number of columns, which is 1.
        Y_vec = Y_ndarray.reshape(-1, 1)
        
        if method == 'OLS': 
            for k in range(K):
                ## estimate B_bk
                mod = Lasso(alpha = reg_path[k], fit_intercept=False)
                mod.fit(X_csr, Y_vec)
                coef = mod.coef_
                ## compute support S_bk
                S_bk = np.nonzero(coef) # trying to make this a list, to see if this allows me to map the set function
                supp_list[k].append(S_bk) #indices related to non-zero support
                
        elif method == 'Logistic':
            for k in range(K):
                ## estimate B_bk
                mod = LogisticRegression(penalty = 'l1', C = 1/reg_path[k], solver = 'liblinear', fit_intercept = False)
                mod.fit(X_csr, Y_vec)
                coef = mod.coef_
                ## compute support S_bk
                # Below, we take the second element [1] because for some reason, 
                # this outputs an array of 0's as the first element, unlike in the OLS case
                S_bk = np.nonzero(coef)[1] # trying to make this a list, to see if this allows me to map the set function
                supp_list[k].append(S_bk) #indices related to non-zero support
                
        else:
            print('Please enter a valid method.\nValid methods include OLS and Logistic')
    
    # compute the aggregared support sets for each value of lambda by intersecting each candidate support for each value of lambda
    S_k = [reduce(np.intersect1d, tuple(supp_list[i])) for i in range(K)]
    
    # drop any elements of S_k which are empty
    S_k_star = [x for x in S_k if x.size != 0]
    sizes = np.array([S_k[i].size for i in range(len(S_k))])
    # index of empty elements of S_k which were dropped
    
    # For which values of k are the support sets empty?
    where = np.where(sizes == 0)[0]
    # If there's any values of k where the support set is empty, print
    if where.size != 0:
        print(f"empty candidate support sets for k = {where}\n or lambda = {[reg_path[i] for i in where]}\n these sets were dropped")
    
    # Remove duplicates from S_k
    S_k_star2 = list(unique_generator(S_k_star))
    
    if len(S_k_star2) != len(S_k_star):
        print('Some candidate support sets appeared more than once. The duplicates were dropped.')
    
    ### Union Step
    K2 = len(S_k_star2)   
    
    # this is where we'll store the one-step forecast error
    F_bk = np.full([B2, K2], np.nan)
    
    # create an empty list to store k* (minimum error lambda-index for each bootstrap resample)
    K_star_b = []
    
    # create a list of B2 many lists
    # each of the B2 lists will be a list of k many constrained OLS estimates
    BETA_bk = [[] for i in range(B2)]
    for i in range(B2):
        bs1 = mbbs(df, n_blocks, block_length)
        bs2 = mbbs(df, n_blocks, block_length)
        
        ## Construct Y*, U* based on boostrap, according to equation (2.2) (Ruiz 2020)
        # Create design matrix for lag d that excludes a column of 1s to explicitly solve for the intercept term
        X1 = VAR_design(bs1, d, n_obs, n_var, intercept)

        # we aren't really vectorizing X, but we'll call it x_vec for consistency with y
        # we'll take the Kronecker of X and I_M, as per Ruiz (2020)
        X1_vec = np.kron(X1, I_M)
    
        # Convert X to a sparse matrix, to hopefully speed up computation time
        X1_csr = sparse.csr_matrix(X1_vec)
    
        # we need to lose the first observation of each neuron to account for the 
        Y1 = bs1.iloc[d:, :]
        # Y needs to be an ndarray to use reshape
        Y1_ndarray = Y1.to_numpy()

        # We're vectorizing Y via the column stacking operator
        # -1 in the reshape function tells python to infer the number of rows based on the specified number of columns, which is 1.
        Y1_vec = Y1_ndarray.reshape(-1, 1)

        ## Now do the same for the second bootstrap sample
        # Create design matrix for lag d that excludes a column of 1s to explicitly solve for the intercept term
        X2 = VAR_design(bs2, d, n_obs, n_var, intercept)
        
        X2_vec = np.kron(X2, I_M)
    
        # Convert X to a sparse matrix, to hopefully speed up computation time
        X2_csr = sparse.csr_matrix(X2_vec)

        # we need to lose the first observation of each neuron to account for the 
        Y2 = bs2.iloc[d:, :]
        # Y needs to be an ndarray to use reshape
        Y2_ndarray = Y2.to_numpy()

        # We're vectorizing Y via the column stacking operator
        # -1 in the reshape function tells python to infer the number of rows based on the specified number of columns, which is 1.
        Y2_vec = Y2_ndarray.reshape(-1, 1)
        
        if method == 'OLS':
            for k in range(K2):
                # For k = 0 to K - 1 do
                X1_csr_con = X1_csr[:, S_k_star2[k]]
                beta_con_ols = inv(X1_csr_con.transpose() @ X1_csr_con) @ (X1_csr_con.transpose() @ Y1_vec)
                beta_bk = np.zeros((n_var*n_var, 1)) # ((n_obs*n_obs, 1))
                beta_bk[S_k_star2[k], :] = beta_con_ols
                BETA_bk[i].append(beta_bk)
                f_bk = np.sum(np.square(Y2_vec - (X2_csr @ beta_bk)))
                F_bk[i, k] = f_bk
            K_star_b.append(F_bk[i].argmin()) # In the case of ties, first element is returned
        
        elif method == 'Logistic':
            for k in range(K2):
                X1_csr_con = X1_csr[:, S_k_star2[k]]
                X1_csr_con

                mod_union = LogisticRegression(penalty = 'none', C = reg_path[k], solver = 'lbfgs', fit_intercept = False)
                mod_union.fit(X1_csr_con, Y1_vec)
                beta_con_logistic = mod_union.coef_
                beta_bk = np.zeros((n_var*n_var, 1))
                beta_bk[S_k_star2[k], :] = beta_con_logistic.T
                BETA_bk[i].append(beta_bk)
        
                Y2_sigmoid = 1/(1+np.exp(X2_csr @ beta_bk)) # np.exp(X2_csr @ beta_bk) / (1 + np.exp(X2_csr @ beta_bk))
                Y2_hat = np.rint(Y2_sigmoid)
                # np.sum(np.square(Y2_vec - Y2_hat)) # This is the sum of misclassifications
                f_bk = 1 - (np.sum(np.square(Y2_vec - Y2_hat)) / len(Y2_vec)) # This is classification accuracy
                F_bk[i, k] = f_bk
            K_star_b.append(F_bk[i].argmin()) # In the case of ties, first element is returned
        
        else:
            print('Please enter a valid method.\nValid methods include OLS and Logistic')

    BETA_b = [sublst[index] for sublst, index in zip(BETA_bk, K_star_b)]
    
    # I previously had the below line, but I think that was wrong, so I'm trying the new one below that
    # a = [np.count_nonzero(BETA_b) for i in BETA_b]
    a = [np.count_nonzero(i) for i in BETA_b]

    q = round(gamma * len(BETA_b))

    ind = np.argpartition(a, -q)[-q:]

    BETA_b_Q = [BETA_b[i] for i in ind]

    BETA_hat = np.mean(BETA_b_Q, axis = 0)
    
    ### Y-hat and redsiduals
    X = VAR_design(df, d, n_obs, n_var, intercept)
    X_vec = np.kron(X, I_M)
    X_csr = sparse.csr_matrix(X_vec)
        
    Y = df.iloc[d:, :]
    Y_ndarray = Y.to_numpy()
    Y_vec = Y_ndarray.reshape(-1, 1)
    
    if method == 'OLS':
        Y_hat_ = X_csr @ BETA_hat
        resid_ = Y_vec - Y_hat_
    
    elif method == 'Logistic':
        Z = X_csr @ BETA_hat
        Y_sigmoid = 1 / (1 + np.exp(-1 * Z)) #1 / (1 + np.exp(Z))
        Y_hat_ = np.rint(Y_sigmoid)
        resid_ = Y_vec - Y_hat_
        
    else: 
        print('Please enter a valid method.\nValid methods include OLS and Logistic')
        
    UoI_dict = {'S_k' : S_k_star2,
                'Beta_hat' : BETA_hat,
                'Y_hat' : Y_hat_,
                'resid' : resid_,
                'k-star' : K_star_b,
                'F_bk' : F_bk}
    
    return UoI_dict


# In[ ]:


lambda_grid = np.linspace(0.0001, 0.01, 50) # [0.0001, 0.01]

# apply UoI-VAR(1)-LASSO to this chunk of the trial for the VIS data
start = time.time()
UoI_VIS_0_dict = UoI_VAR_LASSO(df = df0, reg_path = lambda_grid, 
                               B1 = 5, d = 1, intercept = False, method = 'OLS',
                               n_blocks = 10, block_length = 10, B2 = 5, gamma = 1)
end = time.time()
run_time = end - start
print(f"time to run (secs): {run_time}")


# In[ ]:


list(UoI_VIS_0_dict.keys())


# In[ ]:


np.count_nonzero(UoI_VIS_0_dict['Beta_hat'])


# In[ ]:


plt.hist(UoI_VIS_0_dict['Beta_hat'])


# In[ ]:


n_var = df0.shape[1]
A_hat_0_0 = UoI_VIS_0_dict['Beta_hat'].reshape((n_var, n_var))
list(zip(*A_hat_0_0.nonzero()))


# In[ ]:


# Create an array whick contains the index location of each of the untransformed neurons involved in the network 
neuron_indx_loc = np.unique(A_hat_0_0.nonzero(), return_counts = True)
neuron_indx_loc


# In[ ]:


myColors = ((0.8, 0.0, 0.0, 1.0), (0.0, 0.0, 0.0, 0.0), (0.0, 0.0, 0.8, 1.0))
cmap = LinearSegmentedColormap.from_list('Custom', myColors, len(myColors))
A_round_0_0 = A_hat_0_0.copy()
A_round_0_0[A_round_0_0 > 0] = 1
A_round_0_0[A_round_0_0 < 0] = -1
sns.heatmap(A_round_0_0, cmap = cmap)


# In[ ]:


def one_step_holdout(df, A, idx_loc):
    df_pred = df.copy()
    df_pred.iloc[1:, idx_loc] = np.nan
    
    for i in range(1, df_pred.shape[0]):
        X_t = np.matrix(A) @ df_pred.iloc[i-1,:]
        df_pred.iloc[i, idx_loc] = X_t[idx_loc]
        
    plt.plot(df_pred.iloc[:, idx_loc])
    plt.plot(df.iloc[:, idx_loc])


# In[ ]:


df0_test = df0.iloc[69:100, :]


# In[ ]:


for idx in neuron_indx_loc[0]:
    plt.figure()
    one_step_holdout(df0_test, A_hat_0_0, idx)
plt.close()


# In[ ]:




