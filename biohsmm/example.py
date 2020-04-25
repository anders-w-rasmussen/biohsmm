import numpy as np
from matplotlib import pyplot as plt
from model import model_1 as model
import states as states
import emissions as emits

# Let's simulate some observations
region_length = 5000

# Fix random state
np.random.seed(19680801)

# Allow normal obs to be observations that are normally distributed perhaps
# some measurement of regulatory function (see CRISPR-decryptr)

normal_observations = np.random.normal(0, 0.5, region_length)
normal_observations[1000:1100] = np.random.normal(0.5, 0.5, 100)
normal_observations[4400:4600] = np.random.normal(-0.5, 0.5, 200)

# Allow binary observations to be observations of some active
# enhancer / silencer mark

binary_observations = np.zeros(region_length)
binary_observations[np.random.randint(1000, 1100, 10)] = 1
binary_observations[np.random.randint(4400, 4600, 20)] = 1

# Add some background frequency of seeing our histone mark
binary_observations[np.random.randint(0, region_length, 10)] = 1

# Assume some data is missing (replae with NANs)
normal_observations[np.random.randint(0, region_length, 2500)] = np.nan
binary_observations[np.random.randint(0, region_length, 2500)] = np.nan
#
# # Plot these observations
x = np.arange(0, region_length)
plt.subplot(211)
plt.title('Normal Observations')
plt.plot(np.argwhere(np.isfinite(normal_observations == True)), normal_observations[np.argwhere(np.isfinite(normal_observations == True))])
plt.subplot(212)
plt.title('Bernoulli Observations')
plt.plot(np.argwhere(np.isfinite(binary_observations == True)), binary_observations[np.argwhere(np.isfinite(binary_observations == True))], c='r')
plt.show()

# State list
state_list = []

# Emission list for state 1
emit_1 = []
emit_1.append(emits.normal_dist_emit([0, 100, 5, 1], emit_name='norm_emit'))
emit_1.append(emits.binary([5000, 10], emit_name='bern_emit'))
pseudo_1 = np.array([100, 5000])
state_list.append(states.Negative_Binomial('bkg', 3, pseudo_1, emit_1))

# Emission list for state 2 (enhancer)
emit_2 = []
emit_2.append(emits.normal_dist_emit([0.5, 100, 5, 1], emit_name='norm_emit'))
emit_2.append(emits.binary([5000, 500], emit_name='bern_emit'))
pseudo_2 = np.array([100, 2000])
state_list.append(states.Negative_Binomial('enhancer', 3, pseudo_2, emit_2))

# Emission list for state 3 (silencer)
emit_3 = []
emit_3.append(emits.normal_dist_emit([-0.5, 100, 5, 1], emit_name='norm_emit'))
emit_3.append(emits.binary([5000, 500], emit_name='bern_emit'))
pseudo_3 = np.array([100, 2000])
state_list.append(states.Negative_Binomial('silencer', 3, pseudo_3, emit_3))

# Define the pi prior and tmat prior (initial probability and transition matrix)

pi_prior = [1, 1, 1]
tmat_prior = np.ones((3, 3)) - np.identity(3)
hsmm_model = model(pi_prior, tmat_prior, state_list)

# Train the model
start_base = 0
obs_list = [normal_observations, binary_observations]
hsmm_model.train(obs_list, start_base, 12)

# Return marginal probabilities of each state (decode)

enhancer_marg  = hsmm_model.give_gammas(obs_list, start_base, 1)
silencer_marg = hsmm_model.give_gammas(obs_list, start_base, 2)


print(hsmm_model.states[1].emits[0].mu)
print(hsmm_model.states[2].emits[0].mu)
print(np.exp(hsmm_model.states[1].emits[1].p1))
print(np.exp(hsmm_model.states[2].emits[1].p1))


# Plot these observations
plt.subplot(211)
plt.title('Enhancer MargProb')
plt.plot(enhancer_marg)
plt.subplot(212)
plt.title('Silencer MargProb')
plt.plot(silencer_marg)
plt.show()





















