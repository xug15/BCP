import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import train_test_split,GridSearchCV
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_curve,auc

data = pd.read_csv('BreastCancer.csv',sep=',',index_col=0)
label_lut = pd.Series({"benign":0,"malignant":1})
y = label_lut.loc[data["Class"]].values　# 将样本标签由字符串转化为整数表示
data = data.iloc[:,:-1] #选取feature对应的列
data = data.fillna(data.mean(axis=0)) #用均值填充缺失值
X = data.values
X = StandardScaler().fit_transform(X) # Z score scaling

data = pd.read_csv('BreastCancer.csv',sep=',',index_col=0)
label_lut = pd.Series({"benign":0,"malignant":1})
y = label_lut.loc[data["Class"]].values　# 将样本标签由字符串转化为整数表示
data = data.iloc[:,:-1] #选取feature对应的列
data = data.fillna(data.mean(axis=0)) #用均值填充缺失值
X = data.values
X = StandardScaler().fit_transform(X) # Z score scaling

clf = GridSearchCV(LogisticRegression(penalty="l2"),
                   param_grid={"C":[0.01,0.1,1,10]}, cv=5,
                   scoring="roc_auc",refit=True,verbose=4)
selector = RFECV(clf, step=1, 
                 min_features_to_select=3,cv=5,
                 importance_getter=lambda clf:clf.best_estimator_.coef_,
                 scoring="roc_auc",verbose=4)
selector = selector.fit(X_discovery, y_discovery)
print(selector.support_)
# array([ True,  True,  True,  True,  True,  True,  True,  True,  True])
print(selector.ranking_)
# [1 1 1 1 1 1 1 1 1]

clf = GridSearchCV(LogisticRegression(penalty="l2"),
                   param_grid={"C":[100,1000]}, cv=5,
                   scoring="roc_auc",refit=True,verbose=4)
selector = RFECV(clf, step=1, 
                 min_features_to_select=3,cv=5,
                 importance_getter=lambda clf:clf.best_estimator_.coef_,
                 scoring="roc_auc",verbose=4)
selector = selector.fit(X_discovery, y_discovery)
print(selector.support_)
# [ True False  True  True False  True  True False  True]
print(selector.ranking_)
# [1 2 1 1 4 1 1 3 1]

from itertools import combinations
feature_combinations=[]
for i in combinations(list(range(9)), 3):
    feature_combinations.append(list(i))
print(len(feature_combinations))
# 84

from sklearn.base import ClassifierMixin, BaseEstimator
class MaskedLogisticRegression(BaseEstimator, ClassifierMixin):
    def __init__(self,feature_indices=None,**params):
        self.feature_indices = feature_indices
        self.estimator = LogisticRegression(**params)
    def mask(self,X):
        if self.feature_indices is None:
            return X
        else:
            return X[:,self.feature_indices]
    def fit(self, X, y=None):
        self.classes_ = np.unique(y)
        return self.estimator.fit(self.mask(X),y)
    def predict(self, X):
        return self.estimator.predict(self.mask(X))
    def predict_proba(self, X):
        return self.estimator.predict_proba(self.mask(X))

clf = GridSearchCV(MaskedLogisticRegression(),
                   param_grid={"feature_indices":feature_combinations}, cv=5,
                   scoring="roc_auc",refit=True,verbose=4)
clf = clf.fit(X_discovery, y_discovery)
print(list(data.columns[clf.best_params_['feature_indices']]))
# ['Cl.thickness', 'Cell.shape', 'Bare.nuclei']

y_pred_proba = clf.predict_proba(X_validation)[:,1]
fpr, tpr,_ = roc_curve(y_validation,y_pred_proba)
AUROC = auc(fpr, tpr)


plt.figure(figsize=(4,4))
plt.plot(fpr, tpr, '-', color='b', label='Validation AUC of {:.4f}'.format(AUROC), lw=2)
plt.plot([0, 1], [0, 1], '--', color=(0.6, 0.6, 0.6), label='Random Chance')
plt.xlim([-0.01, 1.01])
plt.ylim([-0.01, 1.01])
plt.title('ROC curve of test data')
plt.xlabel('FPR')
plt.ylabel('TPR')
plt.legend(loc='best',fontsize='small')
plt.tight_layout()
plt.show()
plt.close()
