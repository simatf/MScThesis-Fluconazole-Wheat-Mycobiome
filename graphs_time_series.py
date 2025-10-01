import matplotlib as mpl
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from scipy.interpolate import interp1d
import statsmodels.api as sm
import seaborn as sns
import scipy.optimize as opt
from scipy import stats

#------------------------------------------------------------------

data = pd.read_csv("plate readerSF43\Kinetic_OD600_23_SF_250611_new_173126_mod.txt", sep= '\t', header= 0)

# normalize
for letter in ['A', 'B', 'C', 'D', 'E', 'F']:
    for i in range(1, 13):
        col = f'{letter}{i}'
        h_col = f'G{i}'
        if col in data.columns and h_col in data.columns:
            data[col] = data[col] - data[h_col]

# set time
time = pd.to_timedelta(data['Time'])

strain_dict = {'A': 'SF43_rep1', 'B': 'SF43_rep2', 'C': 'SF43_rep3'}

row = 'A'
row2 = 'B'
row3 = 'C'

#concentrations_ = ['0.42 mg/mL', '0.52 mg/mL', '0.625 mg/mL', '0.75 mg/mL', '0.9 mg/mL']
concentrations = ['900 mg/L','750 mg/L' , '625 mg/L', '520 mg/L', '430 mg/L', '200 mg/L', '10 mg/L', '3 mg/L', '0 mg/L']
concRPMI = [1, 2, 3, 4, 5, 6, 7, 8, 10]


for i, conc in enumerate(concRPMI):
    if i == 8: # control
        reps_c = np.column_stack((data.loc[:,f'{row}{conc}'], data.loc[:,f'{row2}{conc}'], data.loc[:,f'{row3}{conc}']))
        mean_c = np.mean(reps_c, axis= 1)
        plt.plot(time/pd.Timedelta(days= 1), mean_c, label= concentrations[i], color= 'dimgray', ls= '--') 
    else: 
        reps = np.column_stack((data.loc[:,f'{row}{conc}'], data.loc[:,f'{row2}{conc}'], data.loc[:,f'{row3}{conc}']))
        mean = np.mean(reps, axis= 1)
        plt.plot(time/pd.Timedelta(days= 1), mean, label= concentrations[i])

#plt.title(strain_dict[row])
plt.ylabel('OD540', fontsize= 16)
plt.xlabel('Time [days]', fontsize= 16)
plt.xticks(fontsize= 14)
plt.yticks(fontsize= 14)
plt.xlim(0, 5)

plt.legend(fontsize= 14)
plt.tight_layout()
plt.show()


### Kill percent after 6 days

# time
continues_time = time/pd.Timedelta(days= 1)
data.iloc[:, 0] = continues_time
print(data.loc[240, 'Time'])
# get values after 5 days for each concentration (stacked 3 reps)
vals5 = np.column_stack((data.loc[240, f'{row}{concRPMI[0]}':f'{row}{concRPMI[-2]}'], data.loc[240, f'{row2}{concRPMI[0]}':f'{row2}{concRPMI[-2]}'], data.loc[240, f'{row3}{concRPMI[0]}':f'{row3}{concRPMI[-2]}']))
#vals5 = data.loc[240, f'{row}{concRPMI[0]}':f'{row}{concRPMI[-2]}']
# get control value
#control_val5 = data.loc[240, f'{row}{concRPMI[-1]}']
control_val5 = np.column_stack((data.loc[240, f'{row}{concRPMI[-1]}'], data.loc[240, f'{row2}{concRPMI[-1]}'], data.loc[240, f'{row3}{concRPMI[-1]}']))
c_val5_mean = np.mean(control_val5, axis= 1)

# calculate kill percentages
kill_perc = (1 - vals5 / c_val5_mean) * 100
#for val in vals5:
#    killp = (1-val/c_val5_mean)*100
#    kill_perc.append(killp)

#kill_perc.append(0.00)
#kill_perc.reverse()
kill_perc = np.flip(kill_perc, axis=0)  

# remove 10
kill_perc = np.delete(kill_perc, 1, axis=0)
# remove 200
kill_perc = np.delete(kill_perc, 1, axis=0)

kill_perc = np.array(kill_perc, dtype= np.float64)
# calc mean, std, sterr
means = np.mean(kill_perc, axis= 1)
std_devs = np.std(kill_perc, axis= 1, ddof= 1)
std_err = std_devs/np.sqrt(kill_perc.shape[1])

# concentrations
conc = [3, 430, 520, 625, 750, 900] # 10 and 200 removed

#define log function
def log_model(x, a, b):
    return a + b * np.log(x)

# fit
params, covariance = opt.curve_fit(log_model, conc, means)
a, b = params

# interpolate
#interp_func = interp1d(perc_kill_list, concentrations, kind="linear", fill_value="extrapolate")
# EC50 estimation
#EC50_interp = interp_func(50)
#print(f"Interpolated EC50: {EC50_interp:.2f} mg/L")

# EC50
EC50_interp = np.exp((50 - params[0]) / params[1])
print(f"Interpolated EC50 (50% kill) = {EC50_interp:.2f} mg/L")

# predictions
predicted = log_model(conc, a, b)

# calc R2
ss_res = np.sum((means - predicted) ** 2)
ss_tot = np.sum((means - np.mean(means)) ** 2)
r2 = 1 - (ss_res / ss_tot)

# statistical significance?
slope, intercept, r_value, p_value, std_err = stats.linregress(np.log(conc), means)

# # interpolate
# interp_func = interp1d(means, conc, kind="linear", fill_value="extrapolate")
# # EC50 estimation
# EC50_interp = interp_func(50)
# print(f"Interpolated EC50: {EC50_interp:.2f} mg/L")

# # linear regression
# X = np.array(conc).reshape(-1, 1)
# X = sm.add_constant(X)
# Y = np.array(means)

# model = sm.OLS(Y, X).fit()
# predictions = model.predict(X)

# print(model.summary())

# intercept, slope = model.params
# print(f"Regression Equation: Kill % = {intercept:.2f} + {slope:.4f} * Concentration")
# p_value = model.pvalues[1]  # p-value for the concentration coefficient
# r_squared = model.rsquared   # R² value

# # Generate smooth fit line for plotting
# x_fit = np.linspace(0, max(conc), 100)  # Smooth X values
# y_fit = intercept + slope * x_fit  # Corresponding Y values

### plot 
plt.figure(figsize=(8, 6))
ax = sns.scatterplot(x= conc, y= means, label= "Data", color= "tab:red", s= 80)
ax2 = sns.lineplot(x= conc, y= predicted, color='tab:blue', label="Log Fit", linewidth=2)
plt.errorbar(conc, means, yerr= std_err, fmt= 'none', color= 'tab:red', capsize= 5)

ax.text(550, 13, f"EC50 = {EC50_interp:.2e} mg/L", fontsize=14, color='black')
ax.text(550, 14.5, f"$p$-value = {p_value:.4f}", fontsize=14, color='black')
ax.text(550, 16, f"$R^2$ = {r2:.4f}", fontsize=14, color='black')

ax.set_xlabel('Fluconazole concentration [mg/L]', fontsize = 16)
ax.set_ylabel('Inhibition percentage [%]', fontsize = 16)
ax.tick_params(axis= 'x', labelsize = 14)
ax.tick_params(axis= 'y', labelsize = 14)

plt.legend(fontsize= 14)
plt.tight_layout()
plt.show()