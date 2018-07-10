import os
import netCDF4 as nc
import numpy as np
import configuration_file as config
import matplotlib.pyplot as plt

i = 0
taus = {}
final_pressures = {}
folder = r'../NCDF'
quantityOfInterest = 'pressure'

for file in os.listdir(folder):
    filepath = os.path.join(folder,file)
    data = nc.Dataset(filepath, 'r')

    setup = data.timestepping + " - " + data.coupling

    if not taus.has_key(setup):
        taus[setup] = []
    if not final_pressures.has_key(setup):
        final_pressures[setup] = []

    if \
            (hasattr(data, 'elasticity_module') and abs(data.elasticity_module - config.E) < 10**-10) \
            and \
            (hasattr(data, 'length') and abs(data.length - config.L) < 10**-5) \
            and \
            (hasattr(data, 'n_elem') and abs(data.n_elem - config.n_elem) < 10**-5):
        final_pressures[setup].append(data.variables[quantityOfInterest][:,-1])
        taus[setup].append(data.tau)
    else:
        print abs(data.length - config.L) < 10**-10
        print abs(data.elasticity_module == config.E) < 10**-10
        print "invalid!"
    data.close()
    i+=1

handles = []
labels = []
markers = ['*','|','_','.','1','2','3','4']
ref_algorithm = config.TimeStepping.TrapezoidalRule.name + " - " + config.CouplingAlgorithm.Monolitic.name

do_not_plot_reference_algorithm_results = False
do_not_plot_most_accurate_solution = True

for setup in taus.keys():
    if setup == ref_algorithm and do_not_plot_reference_algorithm_results:
        continue  # do not plot reference algorithm

    taus_for_setup = np.array(taus[setup])
    i_ref = np.argmin(taus[ref_algorithm])
    p_ref = final_pressures[ref_algorithm][i_ref]

    i = 0
    error_dict = {}

    for p in final_pressures[setup]:
        if abs(taus[setup][i] - taus[ref_algorithm][i_ref]) < 10**-10 and do_not_plot_most_accurate_solution:
            i += 1
            continue

        error_dict[taus[setup][i]] = np.linalg.norm(final_pressures[setup][i] - p_ref)/final_pressures[setup][i].shape[0]
        i += 1


    experiment_taus = error_dict.keys()
    errors = error_dict.values()

    sort_ids = np.argsort(experiment_taus)

    sorted_taus = [experiment_taus[sort_ids[j]] for j in range(len(sort_ids))]
    sorted_errors = [errors[sort_ids[j]] for j in range(len(sort_ids))]

    h = plt.loglog(sorted_taus[:], sorted_errors[:], markers.pop())[0]
    handles.append(h)
    labels.append(setup)

plt.legend(handles, labels)
plt.show()
