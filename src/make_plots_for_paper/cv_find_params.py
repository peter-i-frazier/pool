__author__ = 'jialeiwang'
import subprocess

root_path = '/fs/home/jw865/remote_deployment/ucsd_reversible_labeling/make_plots_for_nature_paper/'
alpha_0_list = [10, 50, 100, 300, 500, 800, 1000, 5000, 10000]
alpha_1_list = [0.01, 0.1, 1, 2, 5, 10, 20, 50]
p_1_list = [0.5, 0.3, 0.1, 1e-2, 1e-3, 1e-4, 1e-5, 1e-6]

for alpha_0 in alpha_0_list:
    for alpha_1 in alpha_1_list:
        for p_1 in p_1_list:
            subprocess.call(['jsub', '"' + root_path + 'cv_find_params.nbs', root_path + 'gen_roc_data.R',
                             str(alpha_0), str(alpha_1), str(p_1), '0"', '-mfail', '-email', 'charleswang304@gmail.com',
                             '-xhost', 'sys_pf'])
