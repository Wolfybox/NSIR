import numpy as np
import torch

# print(np.random.randint(low=0, high=100, size=1))
# print(np.log(np.e))
# print(np.round(1.6))
# a = np.array([1, 2, 3])
# b = np.array([2, 3, 4])
# a[np.array([0, 2])] = b[:2]
# print(a)
# print(np.intersect1d(a, b))
# state_dict = {
#     'a': 1,
#     'b': 2,
# }
#
# print(list(state_dict.keys())[0])
# print(int(1.999999999))

a = torch.rand(size=(5, 5))
print(a)
indices = torch.where(a > 0.5)
indices = torch.nonzero(a > 0.5, as_tuple=True)
print(indices)
b = torch.ones(size=(5, 5))
a[indices] = b[indices]
print(a)
