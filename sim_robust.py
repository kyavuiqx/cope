import numpy as np
from scipy.special import expit
from opeuc import OPEUC
from opedr import OPEDR
from simulator_save import Simulator
from policy import target_policy
from utilize import true_q_function, true_ratio_function, true_pm_function, true_pa_function, false_q_function, false_ratio_function, false_pm_function, false_pa_function, true_q_function_drl
from utilize_prototypemodel import QModel, RatioModel, PMModel, PAModel, QModel2

gamma = 0.9
state_dim = 1
simulator = Simulator(dim_state=state_dim, model_type='toy')

q_model0 = QModel(true_q_function)
r_model0 = RatioModel(true_ratio_function)
pm_model0 = PMModel(true_pm_function)
pa_model0 = PAModel(true_pa_function)

drl_q_model = QModel2(true_q_function_drl)
drl_r_model = RatioModel(true_ratio_function)
drl_pa_model = PAModel(true_pa_function)

truth = 55.0559
num_time = 200

nrep = 100
num_trajectory_list = np.arange(1, 16) * 60

error0 = np.zeros((num_trajectory_list.shape[0], nrep))
error1 = np.zeros((num_trajectory_list.shape[0], nrep))
error2 = np.zeros((num_trajectory_list.shape[0], nrep))
error4 = np.zeros((num_trajectory_list.shape[0], nrep))
for i, num_trajectory in enumerate(num_trajectory_list):
    print("Trajectory:", num_trajectory, "start.")
    for r in range(nrep):
        if r % 10 == 0:
            print("Remain: ", 100 - r, "%.")
        pass
        simulator.sample_trajectory(num_trajectory, num_time, seed=r)
        simulator.trajectory2iid()
        sim_iid_dataset = simulator.iid_dataset
        target_action_s0 = np.apply_along_axis(target_policy, 1, sim_iid_dataset['s0']).flatten()
        target_action = np.apply_along_axis(target_policy, 1, sim_iid_dataset['state']).flatten()
        target_action_next = np.apply_along_axis(target_policy, 1, sim_iid_dataset['next_state']).flatten()
        sim_iid_dataset['policy_action_s0'] = target_action_s0
        sim_iid_dataset['policy_action'] = target_action
        sim_iid_dataset['policy_action_next'] = target_action_next

        opeuc_obj0 = OPEUC(sim_iid_dataset, q_model0, r_model0, pm_model0,
                           pa_model0, gamma=gamma, matrix_based_learning=True, policy=target_policy)
        opeuc_obj0.compute_opeuc()
        trl_value0 = opeuc_obj0.opeuc
        error0[i, r] = trl_value0 - truth

        q_model1 = QModel(true_q_function, seed=r)
        r_model1 = RatioModel(false_ratio_function, seed=r)
        pm_model1 = PMModel(true_pm_function, seed=r)
        pa_model1 = PAModel(true_pa_function, seed=r)
        opeuc_obj1 = OPEUC(sim_iid_dataset, q_model1, r_model1, pm_model1,
                           pa_model1, gamma=gamma, matrix_based_learning=True, policy=target_policy)
        opeuc_obj1.compute_opeuc()
        trl_value1 = opeuc_obj1.opeuc
        error1[i, r] = trl_value1 - truth

        q_model2 = QModel(false_q_function, seed=r)
        r_model2 = RatioModel(true_ratio_function, seed=r)
        pm_model2 = PMModel(true_pm_function, seed=r)
        pa_model2 = PAModel(false_pa_function, seed=r)
        opeuc_obj2 = OPEUC(sim_iid_dataset, q_model2, r_model2, pm_model2,
                           pa_model2, gamma=gamma, matrix_based_learning=True, policy=target_policy)
        opeuc_obj2.compute_opeuc()
        trl_value2 = opeuc_obj2.opeuc
        error2[i, r] = trl_value2 - truth

        opedr_obj0 = OPEDR(sim_iid_dataset, drl_q_model, drl_r_model, drl_pa_model, 
                           gamma=gamma, target_policy=target_policy, matrix_based_learning=True)
        opedr_obj0.compute_opedr()
        drl_value0 = opedr_obj0.opedr
        error4[i, r] = drl_value0 - truth
        pass
    print("Trajectory: ", num_trajectory, "done.")
    pass

np.savetxt('triple-robust-trl-truth.csv', error0, delimiter=',')
np.savetxt('triple-robust-trl-correct-1.csv', error1, delimiter=',')
np.savetxt('triple-robust-trl-correct-2.csv', error2, delimiter=',')
np.savetxt('triple-robust-drl-truth.csv', error4, delimiter=',')
