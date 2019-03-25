import pickle

with open("pickle_example_mat", "rb") as fr:
    data = pickle.load(fr)

print(data)