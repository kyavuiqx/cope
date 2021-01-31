import numpy as np
import pandas as pd

def rmse(x):
    return np.sqrt(np.mean(np.square(x)))

def mse(x):
    return np.mean(np.square(x))

def result_display(x, verbose=True, trim=None, file_name=None):
    def bias_est(x):
        if trim is not None:
            x = x[np.abs(x) < trim]
        return np.mean(x)

    def std_est(x):
        if trim is not None:
            x = x[np.abs(x) < trim]
        return np.std(x)

    def rmse_est(x):
        if trim is not None:
            x = x[np.abs(x) < trim]
        return np.sqrt(np.mean(np.square(x)))

    bias = np.apply_along_axis(bias_est, 2, x)
    std = np.apply_along_axis(std_est, 2, x)
    mse = np.apply_along_axis(rmse_est, 2, x)

    if verbose:
        print(bias)
        print(std)
        print(mse)
    else:
        result = np.hstack([bias, std, mse])
        # np.savetxt("result.xlsx", result, delimiter=',')
        df = pd.DataFrame(result, columns=['Bias', 'Std', 'RMSE'])
        if file_name is None:
            file_name = 'result.xlsx'
        else:
            file_name = file_name + '.xlsx'
        df.to_excel(file_name, index=False)
