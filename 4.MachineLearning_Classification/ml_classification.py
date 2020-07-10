import csv
import numpy as np
from sklearn.linear_model import LogisticRegression
from sklearn.svm import SVC

train_data = []
train_label = []
with open("training.csv", "r") as f:
    reader = csv.reader(f)
    for i, row in enumerate(reader):
        if i == 0:
            continue
        tmp = row[1:]
        tmp = [int(i) for i in tmp]
        train_data.append(tmp[:-1])
        train_label.append(tmp[-1])

test_data = []
test_label = []
with open("test.csv", "r") as f:
    reader = csv.reader(f)
    for i, row in enumerate(reader):
        if i == 0:
            continue
        tmp = row[1:]
        tmp = [int(i) for i in tmp]
        test_data.append(tmp[:-1])
        test_label.append(tmp[-1])


svc = SVC()
svc.fit(train_data, train_label)

print(svc.score(test_data, test_label))