import os
import copy

from mogonet.feat_importance import cal_feat_imp, summarize_imp_feat
from mogonet.train_test import train_test

data_folder = 'ICGC'
view_list = [1, 2]


def run_mogonet(num_class):
    os.chdir('./mogonet')
    num_epoch_pretrain = 500
    num_epoch = 2500
    lr_e_pretrain = 1e-3
    lr_e = 5e-4
    lr_c = 1e-3
    train_test(data_folder, view_list, num_class,
               lr_e_pretrain, lr_e, lr_c,
               num_epoch_pretrain, num_epoch)


def run_biomarker(num_class):
    os.chdir('./mogonet')
    model_folder = os.path.join(data_folder, 'models')
    featimp_list_list = []
    for rep in range(5):
        featimp_list = cal_feat_imp(data_folder, os.path.join(model_folder, str(rep + 1)),
                                    view_list, num_class)
        featimp_list_list.append(copy.deepcopy(featimp_list))
    summarize_imp_feat(featimp_list_list)
