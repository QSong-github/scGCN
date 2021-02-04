from __future__ import division
from __future__ import print_function

import datetime
import json
import logging
import os
import pickle
import time

import numpy as np
import optimizers
import torch
from config import parser
from models.base_models import NCModel, LPModel
from data_utils import *
from utils.train_utils import get_dir_name, format_metrics
from utils.data_utils import load_data

args = parser.parse_args()
print (args.dataset)
args.manifold='PoincareBall'
args.model='HGCN'

np.random.seed(args.seed)
torch.manual_seed(args.seed)
if int(args.double_precision):
    torch.set_default_dtype(torch.float64)

args.device ='cpu' #'cuda:' + str(args.cuda) if int(args.cuda) >= 0 else 'cpu'
args.patience = args.epochs if not args.patience else  int(args.patience)

np.random.seed(args.seed)
torch.manual_seed(args.seed)
if int(args.double_precision):
    torch.set_default_dtype(torch.float64)
if int(args.cuda) >= 0:
    torch.cuda.manual_seed(args.seed)
args.device = 'cpu'#'cuda:' + str(args.cuda) if int(args.cuda) >= 0 else 'cpu'
args.patience = args.epochs if not args.patience else  int(args.patience)
logging.getLogger().setLevel(logging.INFO)
if args.save:
    if not args.save_dir:
        dt = datetime.datetime.now()
        date = f"{dt.year}_{dt.month}_{dt.day}"
        models_dir = os.path.join(os.environ['LOG_DIR'], args.task, date)
        save_dir = get_dir_name(models_dir)
    else:
        save_dir = args.save_dir
    logging.basicConfig(level=logging.INFO,
                        handlers=[
                                  logging.FileHandler(os.path.join(save_dir, 'log.txt')),
                                  logging.StreamHandler()
                                ])
    
logging.info(f'Using: {args.device}')
logging.info("Using seed {}.".format(args.seed))

data = load_data(args, os.path.join(os.environ['DATAPATH'],args.dataset))
args.n_nodes, args.feat_dim = data['features'].shape

Model = NCModel
args.n_classes = int(data['labels'].max() + 1)
logging.info(f'Num classes: {args.n_classes}')

args.lr_reduce_freq = args.epochs

# Model and optimizer
model = Model(args)
logging.info(str(model))
optimizer = getattr(optimizers, args.optimizer)(params=model.parameters(), lr=args.lr,
                                                weight_decay=args.weight_decay)
lr_scheduler = torch.optim.lr_scheduler.StepLR(
    optimizer,
    step_size=int(args.lr_reduce_freq),
    gamma=float(args.gamma)
    )
tot_params = sum([np.prod(p.size()) for p in model.parameters()])
logging.info(f"Total number of parameters: {tot_params}")

# Train model
t_total = time.time()
counter = 0
best_val_metrics = model.init_metric_dict()
best_test_metrics = None
best_pred_metrics = None
best_emb = None

for epoch in range(args.epochs):
    t = time.time()
    model.train()
    optimizer.zero_grad()
    embeddings = model.encode(data['features'], data['adj_train_norm'])
    train_metrics = model.compute_metrics(embeddings, data, 'train')
    train_metrics['loss'].backward()
    if args.grad_clip is not None:
        max_norm = float(args.grad_clip)
        all_params = list(model.parameters())
        for param in all_params:
            torch.nn.utils.clip_grad_norm_(param, max_norm)
    optimizer.step()
    lr_scheduler.step()
    if (epoch + 1) % args.eval_freq == 0:
        model.eval()
        embeddings = model.encode(data['features'], data['adj_train_norm'])
        val_metrics = model.compute_metrics(embeddings, data, 'val')
        if model.has_improved(best_val_metrics, val_metrics):
            best_test_metrics = model.compute_metrics(embeddings, data, 'test')
            best_pred_metrics = model.compute_metrics(embeddings, data, 'pred')
            best_emb = embeddings.cpu()
            best_val_metrics = val_metrics
            counter = 0
        else:
            counter += 1
            if counter == args.patience and epoch > args.min_epochs:
                print("Early stopping")
                break


print ("Optimization Finished!")
print ("Total time elapsed: {:.4f}s".format(time.time() - t_total))
if not best_test_metrics:
    model.eval()
    best_emb = model.encode(data['features'], data['adj_train_norm'])
    best_test_metrics = model.compute_metrics(best_emb, data, 'test')
    best_pred_metrics = model.compute_metrics(best_emb, data, 'pred')
print (" ".join(["Val set results:", format_metrics(best_val_metrics, 'val')]))
print (" ".join(["Test set results:", format_metrics(best_test_metrics, 'test')]))
print (" ".join(["Pred set results:", format_metrics(best_pred_metrics, 'pred')]))
