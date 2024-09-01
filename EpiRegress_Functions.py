import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scipy.stats as stats
from scipy.stats import nbinom
from scipy.ndimage import gaussian_filter1d
import pyjags
from create_matrix import accumulator_normalized, accumulator_matrix
from read_cases import local_cases, imported_cases

#serial_distribution: serial interval distribution
#accumulator_matrix: the raw covariate matrix starting from at least |w| days before the start time
#local_cases: incidence curve for local case counts
#imported_cases: incidence curve for imported case counts
#start, end: start and end time for Rt's to be estimated, the last day ("end" day)'s Rt estimate is obtained through extrapolate rather than MCMC inference
#select: list of variables selected to include in the model

start = accumulator_matrix.iloc[0, 0]
end = accumulator_matrix.iloc[-1, -1]
select = accumulator_matrix.columns

#serial interval distributions
sequence = np.arange(1, 51)
mean = 1.132 
std_dev = 0.742
serial_distribution = np.random.normal(mean, std_dev, size=sequence.size)
serial_distribution = serial_distribution / np.sum(serial_distribution)

def EpiRegress(accumulator_matrix, serial_distribution, start, end, select=None, burn_in=30, draw=30, thin=10):
    if select is None:
        # select = list(range(accumulator_matrix.shape[1]))\
        select = list(range(7)) 

    estim_period = np.arange(max(start - len(serial_distribution), 0), end)
    factor = np.std(accumulator_matrix.iloc[estim_period, :], axis=0)
    selected_columns = accumulator_matrix.columns[select]
    estim_period = np.arange(max(start - len(serial_distribution), 0), end)
    mean_acc = np.mean(accumulator_matrix.iloc[estim_period, :], axis=0)
    std_acc = np.std(accumulator_matrix.iloc[estim_period, :], axis=0)
    factor = std_acc.iloc[select]

    accumulator_normalized = (accumulator_matrix.loc[:, selected_columns] - mean_acc) / factor
    # print ("accumulator_normalized", accumulator_normalized)

    init = {
        'a': 0,
        'b': np.zeros(len(selected_columns))
    }
    data = {
        'accumulator_matrix': accumulator_normalized,
        'local_cases': local_cases,
        'imported_cases': imported_cases,
        'inv_serial_distribution': serial_distribution[::-1],
        'n': len(select),
        'm': len(serial_distribution),
        'start': start,
        'end': end,
        'lambda_': 55
    }

    # Define model and sampling
    with open("epiregress_model.txt", "r") as f:
        model_code = f.read()

    variables = ['a', 'b', 'phi', 'tau', 'R', 'mu']
    jagmod = pyjags.Model(model_code, data=data, chains=1, init = init)
    jagmod.sample(burn_in, vars=['a', 'b', 'phi', 'tau', 'R', 'mu'])
    samples = jagmod.sample(draw * thin, vars=['a', 'b', 'phi', 'tau', 'R', 'mu'], thin=thin)

    variable_matrices = {}

    for var in variables:
        sample_shape = samples[var][0].shape  # Get the shape of the first sample set
        var_matrix = np.zeros((len(samples[var]), *sample_shape))  # Initialize with the correct shape
        
        for i, sample_set in enumerate(samples[var]):
            var_matrix[i] = sample_set  # Assign the sample set directly
        variable_matrices[var] = var_matrix  # Store the matrix for the current variable in the dictionary
    

    print(len(variable_matrices['a']))
    
    Rt = np.zeros((variable_matrices['R'].shape[0], 6))

    for i, row in enumerate(variable_matrices['R']):
        Rt[i] = summary_stat(row)

    beta_summary = np.zeros((variable_matrices['b'].shape[0], 6))
    for i, row in enumerate(variable_matrices['b']):
        beta_summary[i] = summary_stat(row)

    # Dynamically generate covariate names
    print("Shape of beta summary: ", beta_summary.shape)
    num_covariates, num_timepoints, num_samples = variable_matrices['b'].shape
    # covariates = [f'covariate_{i+1}' for i in range(num_covariates)]
    covariate_names = accumulator_normalized.columns.tolist()

    # Create a DataFrame for easier handling
    df = pd.DataFrame(beta_summary, columns=['mean', 'ci_95_low', 'ci_50_low', 'median', 'ci_50_high', 'ci_95_high'])
    df['factor'] = covariate_names

    plot_variances(df)

    mean_case_summary = np.zeros((variable_matrices['mu'].shape[0], 6))
    for i, row in enumerate(variable_matrices['mu']):
        mean_case_summary[i] = summary_stat(row)

    sigma = 2  # Adjust the sigma value as needed
    Rt_smoothed = np.zeros_like(Rt)

    for i in range(Rt.shape[1]):
        Rt_smoothed[:, i] = gaussian_filter1d(Rt[:, i], sigma=sigma)
        np.savetxt('Rt_smoothed_mobility&vax_select not not.csv', Rt_smoothed, delimiter=',')

    # case_interval = pd.DataFrame([case_interval_t(t) for t in range(start, end + 1)])
    mu_matrix = variable_matrices['mu']
    tau_matrix = variable_matrices['tau']
    
    # case_interval = np.vstack([case_interval_t(t, mu_matrix, tau_matrix) for t in range(start, end + 1)])

    # loglikelihood_raw = np.vstack([case_prob(t, local_cases, mu_matrix, tau_matrix, start) for t in range(start, end + 1)])
    # loglikelihood_sum_raw = np.sum(loglikelihood_raw, axis=1)
    # deviance_i_c = -2 * np.mean(loglikelihood_sum_raw) + 0.5 * np.mean(np.var(-2 * loglikelihood_sum_raw))


    # print("Beta Summary:")
    # print(beta_summary)

    print("Rt Summary:")
    print(Rt)

    print("Rt Smoothed:")
    print(Rt_smoothed)

    # print("Mean Case Summary:")
    # print(mean_case_summary)

    # print("DIC:")
    # print(deviance_i_c)

    # print("case interval matrix:")
    # print(case_interval)

    return {
        'beta_summary': beta_summary,
        'Rt': Rt_smoothed,
        'mean_case_summary': mean_case_summary
        # 'deviance_i_c': deviance_i_c,
        # 'case_interval':case_interval
    }

def summary_stat(data):
    return np.concatenate(([np.mean(data)], np.quantile(data, [0.025, 0.25, 0.5, 0.75, 0.975])))


def case_interval_t(t, mu_matrix, tau_matrix):
    # Repeat mu for each row of tau or vice versa
    mu_t = np.repeat(mu_matrix, len(tau_matrix), axis=0)
    tau_t = np.tile(tau_matrix, (len(mu_matrix), 1))
    
    # Generate simulated cases for the given day
    sim = np.hstack([
        nbinom.rvs(mu_t[i] / (tau_t[i] - 1), 1 / tau_t[i], size=20)
        for i in range(len(tau_t))
    ])
    
    # Calculate quantiles for the simulated cases
    return np.quantile(sim, [0.025, 0.97])

def case_prob(t, local_cases, mu_matrix, tau_matrix, start):
    log_probs = []
    print("t: ", t)
    for i in range(mu_matrix.shape[1]):  
        mu_t = mu_matrix[t, i]
        tau = tau_matrix[t, i]
        log_prob = nbinom.logpmf(local_cases[t + start - 1], mu_t / (tau - 1), 1 / tau)
        log_probs.append(log_prob)
    return np.array(log_probs)


def plot_variances(covariates):
    num_factors = len(covariates)
    fig, ax = plt.subplots(figsize=(10,6))
    y_pos = np.arange(num_factors)

    #50% intervals
    for i, (ci_low, ci_high) in enumerate(zip(covariates['ci_50_low'], covariates['ci_50_high'])):
        ax.plot([ci_low, ci_high], [y_pos[i], y_pos[i]], color="#718AB1", linewidth=5, alpha=1)

    # 95% intervals
    for i, (ci_low, ci_high) in enumerate(zip(covariates['ci_95_low'], covariates['ci_95_high'])):
        ax.plot([ci_low, ci_high], [y_pos[i], y_pos[i]], color="#718AB1", linewidth=5, alpha=0.3)
        
    ax.set_yticks(y_pos)
    ax.set_yticklabels(covariates['factor'])
    ax.invert_yaxis()
    ax.set_xlabel("impact on $R_t$")
    ax.set_title("impacts of covariates")
    plt.tight_layout()
    plt.savefig('covariate_value(30,mobility&vax select not not).svg', format ='svg')
    plt.show()


def plot_Rt(R2):
    n = R2.shape[0]  # Number of time steps

    # Create the plot using matplotlib
    fig, ax = plt.subplots(figsize=(8, 8))

    ax.fill_between(range(n), R2[:, 1], R2[:, 5], color="#718AB1", alpha=0.3, label="Credible Interval")

    # Plot point estimate (if available)
    if R2.shape[1] > 4:
        ax.plot(R2[:, 3], marker='o', color='#718AB1', linestyle='', label='Point Estimate')

    # Plot threshold line at Rt=1
    ax.axhline(y=1, color='palevioletred', linestyle='--', label='Rt=1')

    # Set labels and title
    ax.set_xlabel('Day')
    ax.set_ylabel('R(t)')
    ax.set_title('Effective Reproductive Number (Rt)')

    plt.xticks(rotation=45)
    plt.legend()
    plt.tight_layout()
    plt.savefig('epiregress_RT.svg', format ='svg')
    plt.show()
