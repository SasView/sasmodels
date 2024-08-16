try:
    import sasmodels
except ModuleNotFoundError:
    raise ModuleNotFoundError("sasmodels is not installed. Please install it using 'pip install sasmodels'.")

from numpy import logspace, sqrt
from matplotlib import pyplot as plt
from sasmodels.core import load_model
from sasmodels.direct_model import call_kernel, call_Fq
import plotly.express as px
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from plotly.offline import plot

import pandas
import numpy as np
from itertools import product
from timeit import timeit

model_name = 'cylinder2'
try:
    model_cylinder_accurate = load_model(model_name)
except Exception as e:
    print(f"Error: {e}")
    print(f"Model {model_name} not found. Please check if the model exists in sasmodels.")
    raise e

q = np.linspace(0.001, 0.1, 200)
kernel_cylinder = model_cylinder_accurate.make_kernel([q])

n = 30
lengths = np.geomspace(10, 1000000, n)
# lengths = lengths[::4]
# rtols = [0.00001]

Iq_data = {length: {} for length in lengths}

radii = np.geomspace(10, 1000000, n)

res = dict()
n_runs = 10

def time_func(lr):
    length, radius = lr
    pars_cylinder = {'radius': radius, 'length': length, 'scale': 1, 'background': 0.001}
    Iq_data[length][radius] = call_kernel(kernel_cylinder, pars_cylinder)
    res[(length, radius)] = timeit(lambda: call_kernel(kernel_cylinder, pars_cylinder), number=n_runs) / n_runs

for lr in product(lengths, radii):
    time_func(lr)


index = pandas.MultiIndex.from_tuples(res.keys(), names=['length', 'radius'])

df = pandas.DataFrame(res.values(), index=index, columns=['time'])
df.to_excel(f'{model_name}_benchmark.xlsx')
