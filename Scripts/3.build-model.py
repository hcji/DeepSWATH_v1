# -*- coding: utf-8 -*-
"""
Created on Wed Jan  1 06:55:41 2020

@author: yn
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import normalize, MinMaxScaler
from sklearn.metrics import roc_curve, auc, accuracy_score, confusion_matrix, precision_score, recall_score
from tensorflow.keras.models import Model, load_model
from tensorflow.keras.layers import Dense, Input, Dropout, Conv1D, Flatten, concatenate
from tensorflow.keras.optimizers import Adam
from tensorflow.keras.callbacks import ModelCheckpoint, EarlyStopping, ReduceLROnPlateau


class DIANet:
    def __init__(self, X1, X2, Y):
        self.X1 = np.expand_dims(X1,-1)
        self.X2 = np.expand_dims(X2,-1)
        self.Y = Y
        
        dim_x = X1.shape[1]     
        
        inp_x1 = Input(shape=(dim_x,1))
        inp_x2 = Input(shape=(dim_x,1))
        
        # conv chrom x1
        hid_x1 = Conv1D(32, 2, activation='relu', input_shape=(dim_x,1))(inp_x1)
        hid_x1 = Conv1D(32, 2, activation='relu')(hid_x1)
        hid_x1 = Flatten()(hid_x1)
        
        # conv chrom x2
        hid_x2 = Conv1D(32, 2, activation='relu', input_shape=(dim_x,1))(inp_x2)
        hid_x2 = Conv1D(32, 2, activation='relu')(hid_x2)
        hid_x2 = Flatten()(hid_x2)
        
        # concat
        concat = concatenate([hid_x1, hid_x2])
        hid = Dense(16, activation='relu')(concat)
        hid = Dense(16, activation='relu')(hid)
        
        # output
        pred = Dense(2, activation='softmax')(hid)
        
        # optimizer
        opt = Adam(lr=0.001)
        
        # build model
        model = Model(inputs=[inp_x1, inp_x2], outputs=pred)
        model.compile(optimizer=opt, loss='categorical_crossentropy', metrics=['accuracy'])
        self.model = model
        
    def train(self, epochs=50, save_path = 'Model/DeepDIA_Model.h5'):
        X1 = self.X1
        X2 = self.X2
        Y = self.Y
        
        # call back
        earlyStopping = EarlyStopping(monitor='val_loss', patience=5, verbose=0, mode='min')
        mcp_save = ModelCheckpoint(save_path, save_best_only=True, monitor='val_loss', mode='min')
        reduce_lr_loss = ReduceLROnPlateau(monitor='val_loss', factor=0.1, patience=3, verbose=1, epsilon=1e-4, mode='min')
        
        history = self.model.fit([X1, X2], [Y], epochs=epochs, callbacks=[earlyStopping, mcp_save, reduce_lr_loss], validation_split=0.1)
        return history


def plot_roc(y_pred, y_real, classes):
    fpr = dict()
    tpr = dict()
    roc_auc = dict()
    for i in range(2):
        fpr[i], tpr[i], _ = roc_curve(y_real[:, i], y_pred[:, i], drop_intermediate=False)
        roc_auc[i] = auc(fpr[i], tpr[i])
    '''
    fpr["micro"], tpr["micro"], _ = roc_curve(y_pred.ravel(), y_real.ravel())
    roc_auc["micro"] = auc(fpr["micro"], tpr["micro"])
    '''
    plt.figure()
    lw = 2
    for i in range(len(fpr)):
        plt.plot(fpr[i], tpr[i], lw=lw, label=classes[i] + ': (area = {})'.format(round(roc_auc[0], 4)))
    plt.plot([0, 1], [0, 1], color='navy', lw=lw, linestyle='--')
    plt.xlim([0.0, 1.0])
    plt.ylim([0.0, 1.05])
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    plt.legend(loc="lower right")
    plt.show()    


if __name__ == '__main__':
    
    simu_precursor_eics = np.load('Data/simu_precursor_eics.npy')
    simu_fragment_eics = np.load('Data/simu_fragment_eics.npy')
    simu_decoy_eics = np.load('Data/simu_decoy_eics.npy')
    
    exp_precursor_eics = np.load('Data/all_precursor_eics.npy')
    exp_fragment_eics = np.load('Data/all_fragment_eics.npy')
    exp_decoy_eics = np.load('Data/all_decoy_eics.npy')
    
    precursor_eics = np.vstack((simu_precursor_eics, exp_precursor_eics))
    fragment_eics = np.vstack((simu_fragment_eics, exp_fragment_eics))
    decoy_eics = np.vstack((simu_decoy_eics, exp_decoy_eics))
    
    del(simu_precursor_eics)
    del(simu_fragment_eics)
    del(simu_decoy_eics)
    del(exp_precursor_eics)
    del(exp_fragment_eics)
    del(exp_decoy_eics)

    precursor_eics = normalize(precursor_eics, axis=1, norm='max')
    fragment_eics = normalize(fragment_eics, axis=1, norm='max')
    decoy_eics = normalize(decoy_eics, axis=1, norm='max')
    
    X1 = np.vstack((precursor_eics, precursor_eics))
    X2 = np.vstack((fragment_eics, decoy_eics))
    Y = np.append( np.ones(len(fragment_eics)), np.zeros(len(decoy_eics)))
    Y = np.array(pd.get_dummies(Y))
    
    del(precursor_eics)
    del(fragment_eics)
    del(decoy_eics)
    
    # train test split
    inds = np.arange(len(Y))
    np.random.shuffle(inds)
    n = int(len(Y) * 0.9)
    tr, ts = inds[:n], inds[n:]

    X1_tr, X1_ts = X1[tr], X1[ts]
    X2_tr, X2_ts = X2[tr], X2[ts]
    Y_tr, Y_ts = Y[tr], Y[ts]
    
    del(X1)
    del(X2)
    del(Y)
    
    mod = DIANet(X1_tr, X2_tr, Y_tr)
    history = mod.train(epochs=50, save_path='Model/DeepDIA_Model.h5')
    
    mod = load_model('Model/DeepDIA_Model.h5')
    X1_ts = np.expand_dims(X1_ts,-1)
    X2_ts = np.expand_dims(X2_ts,-1)
    Y_pred = mod.predict([X1_ts, X2_ts])
    roc = roc_curve(Y_ts[:, 0], Y_pred[:, 0], drop_intermediate=False)
    roc_auc = auc(roc[0], roc[1])
    
    Y_pred = np.round(Y_pred)
    acc = accuracy_score(Y_ts[:,0], Y_pred[:,0])
    confusion = confusion_matrix(Y_ts[:,0], Y_pred[:,0])
    precision = precision_score(Y_ts[:,0], Y_pred[:,0])
    recall = recall_score(Y_ts[:,0], Y_pred[:,0])
    
    fig, axes = plt.subplots(nrows=2, ncols=2, figsize=(12, 8))

    axes[0,0].plot(np.arange(1, len(history.epoch)+1), history.history['loss'], alpha= 0.8)
    axes[0,0].plot(np.arange(1, len(history.epoch)+1), history.history['val_loss'], alpha= 0.8)
    axes[0,0].set_ylabel('loss')
    axes[0,0].set_xlabel('epoch')
    axes[0,0].legend(['train', 'valid'], loc="lower left")
    
    axes[0,1].plot(np.arange(1, len(history.epoch)+1), history.history['acc'], alpha= 0.8)
    axes[0,1].plot(np.arange(1, len(history.epoch)+1), history.history['val_acc'], alpha= 0.8)
    axes[0,1].set_ylabel('accuracy')
    axes[0,1].set_xlabel('epoch')
    axes[0,1].legend(['train', 'valid'], loc="lower right")
    
    axes[1,0].plot(np.arange(1, len(history.epoch)+1), history.history['lr'], color = 'green', alpha= 0.8)
    axes[1,0].set_ylabel('learning rate')
    axes[1,0].set_xlabel('epoch')
    axes[1,0].legend(['learning rate'], loc="lower left")

    axes[1,1].plot(roc[0], roc[1], color = 'red', alpha= 0.8)
    axes[1,1].plot([0, 1], [0, 1], color = 'black', alpha= 0.8, linestyle='--')
    axes[1,1].set_ylabel('true positive rate')
    axes[1,1].set_xlabel('false negative rate')
    axes[1,1].legend(['ROC curve (area = %0.2f)' % roc_auc], loc="lower right")
    plt.savefig("Figure/model_eval.png")
    