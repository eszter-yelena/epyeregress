import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import lognorm

#input
#w: serial interval distribution
#X: the raw co variate matrix starting from at least |w| days before the start time
#I: incidence curve for local case counts
#I_imp: incidence curve for imported case counts
#start, end: start(>=2) and end time for Rt's to be estimated, the last day ("end" day)'s Rt estimate is obtained through extrapolate rather than MCMC inference
#select: list of variables selected to include in the model
#draws: number of draws used for point estimates and credible intervals
#thin: thin in MCMC iterations

#output
#Rt: Rt point estimates (posterior mean and median) and credible intervals (50% and 95%, in the form of quantiles) for each day from the start day to the end day
#beta_estimates: point estimates (posterior mean and median) and credible intervals (50% and 95%, in the form of quantiles)
#dic: DIC value of the model, used to assess the fit

from create_matrix import accumulator_normalized, accumulator_matrix
from read_cases import local_cases, imported_cases
from EpiRegress_Functions import EpiRegress, summary_stat, case_interval_t, case_prob, plot_Rt, plot_variances

start = 30  
end = accumulator_matrix.shape[0]  
select = list(range(accumulator_matrix.shape[1]))  # select all columns
# select = list(range(7))   # select all columns


print("getting serial distribution...")
# Serial interval distribution
w = lognorm.rvs(1.132, 0.742, size=50)
w = w / np.sum(w)

print("getting output ...")

# w = np.exp(-(np.arange(1, 51) - 1 - 1.132)**2 / (2 * 0.742**2))
# w = w / np.sum(w)


# Step 3: Call the EpiRegress function
# Assuming EpiRegress is defined and imported correctly
out = EpiRegress(accumulator_matrix, w, start, end, select)

# Step 4: Plot the Rt estimates
# Assuming 'plot' function is defined to plot Rt estimates
plot_Rt(out['Rt'])
# plot_variances()
# plt.xlabel('Time')
# plt.ylabel('Rt')
# # plt.title('Rt Estimates')
# plt.savefig('Rt_estimates.png', dpi=300)
# plt.show()