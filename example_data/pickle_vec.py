import pickle

with open("pickle_example", "rb") as fr:
    data = pickle.load(fr)

print(data)