import os,glob
import geopandas as gpd
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
# import rasterio
import elapid
from sklearnex import patch_sklearn, unpatch_sklearn
import re
import difflib

# import machine learning classifiers
from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import LogisticRegression
# from sklearn.ensemble import ExtraTreesClassifier
# from sklearn.linear_model import TweedieRegressor
from xgboost import XGBClassifier
from lightgbm import LGBMClassifier
from pyimpute import impute,load_training_vector,load_targets
from sklearn import model_selection
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import train_test_split
####import models#############
CLASS_MAP = {
    'rf': (RandomForestClassifier()),
    # 'et': (ExtraTreesClassifier()),
    'xgb': (XGBClassifier()),
    'maxent': (LogisticRegression()),
    # 'lgbm': (LGBMClassifier())
    }


###read monthly GBIF data#######
allPath = glob.glob('/root/autodl-tmp/GBIF_Data/monthData/*.csv')
auc_scoreLsit = []
####multi cores speed calculation###
patch_sklearn()
###climate data#######
raster_features = glob.glob('/root/autodl-tmp/clim_resampleData/*.tif')
for path in allPath:
    try:
        baseName = os.path.basename(path)
        spName = re.findall('[A-Z].*(?=_)',baseName)[0]
        month = re.findall('[0-9]{1,2}',baseName)[0]
        
        ##########read raster data###########
        train_data = pd.read_csv(path).dropna().to_numpy()
        target_xs, raster_info = load_targets(raster_features)
        # ########divide train and test datasets#########
        X_train, X_test, y_train, y_test = train_test_split(train_data[:,0:6],train_data[:,6],shuffle=True,random_state=42)
        ##create species folder##
        spPath = '/root/autodl-tmp/result/model/'+spName+'/'
        if not os.path.exists(spPath):
                os.mkdir(spPath)
        for name, (model) in CLASS_MAP.items():
            # cross validation for accuracy scores (displayed as a percentage)
            k = 5 # k-fold
            kf = model_selection.KFold(n_splits=k)
            accuracy_scores = model_selection.cross_val_score(model, X_train, y_train, cv=kf, scoring='accuracy',n_jobs=-1)
            # print(name + " %d-fold Cross Validation Accuracy: %0.2f (+/- %0.2f)"
            #     % (k, accuracy_scores.mean() * 100, accuracy_scores.std() * 200))

            # spatial prediction
            model.fit(X_train, y_train)
            auc_score = roc_auc_score(y_test, model.predict(X_test))
            # print('auc score: %0.2f', auc_score)
            auc = {spName+name+month:auc_score,'number':train_data.shape[1]}
            auc_scoreLsit.append(auc)
            filePath = '/root/autodl-tmp/result/model/'+spName+'/'+name+'_'+ month
            if not os.path.exists(filePath):
                os.mkdir(filePath)
                impute(target_xs, model, raster_info, outdir=filePath,
                    class_prob=True, certainty=True)
    except Exception as e:
        print(e)
with open('/root/autodl-tmp/result/auc_scoreLsit.txt', 'w') as f:
    for item in auc_scoreLsit:
        f.write("%s\n" % item)







