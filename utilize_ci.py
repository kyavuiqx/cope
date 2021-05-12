import pandas as pd
import numpy as np

def compute_std(num, time, data_point_estimation):
    trajectory_id = np.repeat(np.arange(num), time)
    df = pd.DataFrame(data={'id': trajectory_id, 'trajectory_estimation': data_point_estimation})
    df = df.groupby('id').mean()
    std_value = df['trajectory_estimation'].to_numpy().std()
    return std_value

def compute_ci(num, time, data_point_estimation):
    std_value = compute_std(num, time, data_point_estimation)
    se_value = std_value / np.sqrt(np.array(num*1.0))
    est = np.mean(data_point_estimation)
    ci = np.array([est - 1.96*se_value, est + 1.96*se_value])
    return ci

def cover_truth(ci, truth):
    cover = truth >= ci[0] and truth <= ci[1]
    cover = 1.0 * cover
    return cover
