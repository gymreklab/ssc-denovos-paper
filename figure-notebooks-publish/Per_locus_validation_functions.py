# The functions in this file are used to generate the following from Mitra, et al.:
# Extended Data Figure 10a: Comparison of true versus inferred selection coefficients
# Extended Data Figure 10b: Power to detect negative selection as a function of s

# This file contains helper functions for SISTR per-locus validation.

########## Imports ##########

from LRT_functions import *
from scipy import stats

########## Per-locus Validation Functions ##########

"""Validate the per-locus method

Parameters
----------
per : int
    Period to validate
opt_allele: string
    Optimal alleles to validate (separated by commas)
s_vals: string
    S values to validate (separated by commas)
use_het: string
    Whether to use heterozygosity as a summary statistic
use_common : string
    Whether to use number of common alleles as a summary statistic
use_bins: string
    Whether to use bins of allele frequencies as a summary statistic
num_bins: int
    Number of bins to use when binning allele frequencies
abc_model: string
    Prior/demographic model to use for ABC
lrt_model : string
    Demographic model to use for LRT
ground_truth_model: string
    The ground truth model for simulations
all_pers: boolean 
    Whether to validate all periods (2,3,4) or just one
first_200: boolean
    Whether to infer s for first 200 allele frequencies in file
trials_2000: boolean
    Whether to infer s for 2000 allele frequencies (default: infer s for 200 allele frequencies)
LRT_200: boolean
    Whether to perform 200 simulations for the LRT (default: perform 2000 simulations for the LRT)
    
Returns
-------
opt_allele_list: list
    List of optimal alleles validated 
s_vals_dic: dictionary
    Contains information about inferred s values 
errors_s_dic: dictionary
    Contains information about error for inferred s values
s_vals: list
    List of s values validated
p_vals_dic: dictionary
    Contains information about p values obtained from LRT
errors_p_dic: dictionary
    Contains information about error for p values
eps_bins: float
    Error threshold for bins summary statistic
LogLR_vals_dic: dictionary
    Contains information about LogLR values
"""
def validate_per_locus(per, opt_allele, s_vals, use_het, use_common, use_bins, \
                       num_bins, abc_model, lrt_model, ground_truth_model, all_pers=False, \
                       first_200=True, trials_2000=False, LRT_200=False):
    
    # Process list of optimal alleles
    opt_allele_list = list(opt_allele.split(','))
    opt_allele_list = list(map(int, opt_allele_list)) 
    
    # Process list of s values
    s_vals = list(s_vals.split(','))
    s_vals = list(map(float, s_vals))
    
    # ABC parameters
    constant_het = 0.005
    denom_het = 3
    constant_common = 1
    denom_common = 4
    eps_bins = 0.3
    
    # LRT parameters
    # Default: 2000 sims for LRT
    LRT_num_sims = 2000 
    if LRT_200 == True:
        LRT_num_sims = 200
        
    # Each dictionary contains values for all optimal alleles
    # Key: optimal allele
    # Value: list of mean values for each s 
    s_vals_dic = {}
    errors_s_dic = {}
    p_vals_dic = {}
    errors_p_dic = {}
    LogLR_vals_dic = {}
    
    # Initialize dictionaries above
    for opt_allele in opt_allele_list:
        s_vals_dic[opt_allele] = []
        errors_s_dic[opt_allele] = []
        p_vals_dic[opt_allele] = []
        errors_p_dic[opt_allele] = []
        LogLR_vals_dic[opt_allele] = []
      
    # Run ABC and LRT for each opt_allele
    for opt_allele in opt_allele_list:
        # If plotting all periods on one graph, assign each optimal allele a period
        if all_pers == True:
            per = 2
            if opt_allele == 16 or opt_allele == 40:
                per = 1
            if opt_allele == 5 or opt_allele == 13 or opt_allele == 6 or opt_allele == 12:
                per = 3
            if opt_allele == 7 or opt_allele == 10 or opt_allele == 8 or opt_allele == 9:
                per = 4
                
        LogLR_list = []
        l0_list = []
        ls_list = []
        print('Running per: ' + str(per) + ' optimal allele: ' + str(opt_allele))
        no_ABC_accept = 0 # Number of s values with <10 ABC acceptances
        
        ### Get list of s values in LRT simulations file ###
        s_list_available = []
        lrtFile = '/gymreklab-tscc/bonnieh/lrt/results/' + lrt_model + '/' + str(per) + '_' + str(opt_allele) + '_freqs.txt'
        lrt_file = open(lrtFile, 'r')
        header = lrt_file.readline().strip()
    
        for line in lrt_file:
            info = line.strip().split('\t')
            s = float(info[0])
            if s not in s_list_available:
                s_list_available.append(s)
            
        est_s_dic = {} # Dictionary of estimated s values; Key: True value of s, Value: list of estimated s values
        
        # Put summary statistics in dictionaries
        obs_het_dic = {} # Key - s; value - list of het
        obs_comm_dic = {} # Key - s; value - list of number of common alleles (frequency >= 5%)
        obs_bins_dic = {} # Key - s; value - list of bins (bins are given as lists)
        
        # Same as dictionaries above except without the simulations 
        # where s could not be estimated using ABC (<10 ABC acceptances)
        # These are the simulations used for LRT
        obs_het_dic_lrt = {}
        obs_comm_dic_lrt = {}
        obs_bins_dic_lrt = {}
    
        lrtFile_ground_truth = '/gymreklab-tscc/bonnieh/lrt/results/' + ground_truth_model + '/' + str(per) + '_' + str(opt_allele) + '_freqs.txt'
        
        # Fill in dictionaries of summary statistics
        for s in s_vals:
            
            obs_het_dic[s] = []
            obs_comm_dic[s] = []
            obs_bins_dic[s] = []
            
            # Default: Get first 200 allele frequencies for given s 
            # Later, s will be inferred for these 200 allele frequencies
            freqs_list_raw = GetLRTListFreq200(lrtFile_ground_truth, s) # Get list of allele frequencies for this s
            if first_200 == False:
                freqs_list_raw = GetLRTListFreq200(lrtFile_ground_truth, s, False)
            if trials_2000 == True:
                freqs_list_raw = GetLRTListFreq(lrtFile_ground_truth, s)
                
            # Get summary statistics from frequencies
            for freq_string in freqs_list_raw:
                
                obs_het, obs_common, obs_bins = GetSummStats(freq_string, num_bins)
                obs_het_dic[s].append(obs_het)
                obs_comm_dic[s].append(obs_common)
                obs_bins_dic[s].append(obs_bins)
            
        abcFile = '/gymreklab-tscc/bonnieh/abc/results/' + abc_model +'/' + str(per) + '_' + str(opt_allele) + '.txt' 

        # Read abcFile line by line and place in lookup table in the form of a list
        abc_list = GetABCList(abcFile, num_bins)
        
        # Perform ABC
        for s in s_vals:
            
            list_est_s = [] # List of posterior estimates of s
            
            # Lists of summary statistics of s
            list_het = []
            list_common = []
            list_bins = []
            
            for i in range(0, len(obs_het_dic[s])):
                obs_het = float(obs_het_dic[s][i])
                obs_common = int(obs_comm_dic[s][i])
                obs_bins = obs_bins_dic[s][i]
                
                s_ABC, lower_bound, upper_bound, num_accepted, s_accepted = Get_S_ABC(abc_list, \
                                       obs_het, obs_common, obs_bins, constant_het, \
                                       denom_het, constant_common, denom_common, eps_bins, use_het, \
                                       use_common, use_bins)
                if s_ABC != -1:
        
                    s_ABC_round = get_LRT_bin(s_ABC)
            
                    # Get nearest s in LRT file to perform LRT
                    if s_ABC_round not in s_list_available:
                        s_ABC_round = getNearestS(s_ABC_round, s_list_available)
                        
                    list_est_s.append(s_ABC_round) 
                    list_het.append(obs_het)
                    list_common.append(obs_common)
                    list_bins.append(obs_bins)
                    
                else:
                    no_ABC_accept = no_ABC_accept + 1
                    
            est_s_dic[s] = list_est_s
            obs_het_dic_lrt[s] = list_het
            obs_comm_dic_lrt[s] = list_common
            obs_bins_dic_lrt[s] = list_bins
            
            # Put mean of esimated s in s_vals_dic and calculate standard deviation
            s_vals_dic[opt_allele].append(np.mean(list_est_s)) 
            std_err = np.std(list_est_s)
            errors_s_dic[opt_allele].append(std_err)
            
        # Get LRT summary statistic tables for s = 0
        lrtFile_for_s_0 = '/gymreklab-tscc/bonnieh/lrt/results/' + lrt_model + '/' + str(per) + '_' + str(opt_allele) + '_15_freqs.txt' 
        freqs_list_raw_0 = GetLRTListByRow(lrtFile_for_s_0, 0)
    
        LRT_table_0_het = []
        LRT_table_0_common = []
        LRT_table_0_bins = []
        
        # Get summary statistics from allele frequencies
        for freq_string in freqs_list_raw_0:
                
            obs_het, obs_common, obs_bins = GetSummStats(freq_string, num_bins)
            LRT_table_0_het.append(obs_het) 
            LRT_table_0_common.append(obs_common) 
            LRT_table_0_bins.append(obs_bins)
                
        # Perform LRT
        for s in est_s_dic:
            p_vals_list = []
            LRlog_list = []
            obs_het_list = obs_het_dic_lrt[s]
            obs_common_list = obs_comm_dic_lrt[s]
            obs_bins_list = obs_bins_dic_lrt[s]
            s_ABC_list = est_s_dic[s]
            
            for i in range(0, len(obs_het_list)):
                obs_het = obs_het_list[i]
                obs_common = obs_common_list[i]
                obs_bins = obs_bins_list[i]
                s_ABC_round = s_ABC_list[i]
                
                # Get LRT summary statistic tables for s = s_ABC_round
                
                # If estimated s from ABC is 0, LRT summary statistic tables from file with simulations for s = 0
                if s_ABC_round == 0:
                    freqs_list_raw_s = GetLRTListByRow(lrtFile_for_s_0, 1)
                else:
                    freqs_list_raw_s = GetLRTListFreq(lrtFile, s_ABC_round)
            
                LRT_table_s_het = []
                LRT_table_s_common = []
                LRT_table_s_bins = []
                
                # Get summary statistics from allele frequencies
                for freq_string in freqs_list_raw_s:
                    
                    obs_het_s_ABC, obs_common_s_ABC, obs_bins_s_ABC = GetSummStats(freq_string, num_bins)
                    LRT_table_s_het.append(obs_het_s_ABC) 
                    LRT_table_s_common.append(obs_common_s_ABC) 
                    LRT_table_s_bins.append(obs_bins_s_ABC)
                
                LRT_table_0_het = LRT_table_0_het[0:LRT_num_sims]
                LRT_table_s_het = LRT_table_s_het[0:LRT_num_sims]
                LRT_table_0_common = LRT_table_0_common[0:LRT_num_sims]
                LRT_table_s_common = LRT_table_s_common[0:LRT_num_sims]
                LRT_table_0_bins = LRT_table_0_bins[0:LRT_num_sims]
                LRT_table_s_bins = LRT_table_s_bins[0:LRT_num_sims]
                
                # Get p value
                if len(LRT_table_s_het) != 0:
                    likelihood_0, likelihood_s_ABC, LR, LogLR, pval = LikelihoodRatioTest(LRT_table_0_het, LRT_table_s_het, \
                                LRT_table_0_common, LRT_table_s_common, LRT_table_0_bins, LRT_table_s_bins, LRT_num_sims, \
                                obs_het, obs_common, obs_bins, constant_het, denom_het, \
                                constant_common, denom_common, eps_bins, use_het, use_common, use_bins)
                
                    p_vals_list.append(pval) 
                    LRlog_list.append(LogLR)
                    
                    if s == 0:
                        LogLR_list.append(LogLR)
                        l0_list.append(likelihood_0)
                        ls_list.append(likelihood_s_ABC)
               
            # Put power of p values in p_vals_dic 
            LogLR_vals_dic[opt_allele] = LogLR_list
            p_vals_dic[opt_allele].append(len([i for i in p_vals_list if i < 0.05])/len(p_vals_list)*100) 
            std_err = stats.sem(p_vals_list, ddof = 0)
            errors_p_dic[opt_allele].append(std_err)
            
    return opt_allele_list, s_vals_dic, errors_s_dic, s_vals, p_vals_dic, errors_p_dic, eps_bins, LogLR_vals_dic