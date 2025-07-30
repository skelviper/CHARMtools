import torch
from torch.utils.data import DataLoader, TensorDataset
from torch import nn, optim
from sklearn.model_selection import train_test_split
from scipy import stats
from sklearn import linear_model
import shap
import numpy as np
import pandas as pd
import tqdm
from joblib import Parallel, delayed

class MultiFeatureRegression(nn.Module):
    def __init__(self, input_dim, alpha,norm_type='None'):
        super().__init__()
        self.linear = nn.Linear(input_dim, 1)
        self.alpha = alpha  
        self.norm_type = norm_type
        
        nn.init.constant_(self.linear.weight, 0.0)
        nn.init.constant_(self.linear.bias, 0.0)
        
    def forward(self, x):
        return torch.exp(self.linear(x))
    
    def regularization(self):
        if self.norm_type == 'l1':
            return self.alpha * torch.sum(torch.abs(self.linear.weight))
        elif self.norm_type == 'l2':
            return self.alpha * torch.sum(self.linear.weight**2)
        else:
            raise ValueError(f"Unknown norm type: {self.norm_type}")

class GeneRegressor:
    """
    Use epigenomic feature to predict gene expression in single cell
    """
    def __init__(self, alpha=0.01, patience=20, norm_type='l1'):
        self.alpha = alpha
        self.patience = patience # for early stopping
        self.norm_type = norm_type
        
    def build_model(self, input_dim, alpha):
        return MultiFeatureRegression(input_dim, alpha, self.norm_type)
    
    def train_test_model(self, X,y,epochs,test_size=0.2,lr=1e-3,batch_size=16,patience=None,celltype:list=None):
        if patience is None:
            patience = self.patience
        
        #X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=test_size, random_state=42)
        if celltype is not None:
            celltype_np = np.array(celltype)
            X_train, X_test, y_train, y_test, celltype_train, celltype_test = train_test_split(
                X, y, celltype_np, test_size=test_size, random_state=42
            )
        else:
            X_train, X_test, y_train, y_test = train_test_split(
                X, y, test_size=test_size, random_state=42
            )
            celltype_test = None # 确保在没有提供celltype时，该变量存在

        model = self.build_model(X_train.shape[1], self.alpha)
        optimizer = optim.Adam(model.parameters(), lr=lr)

        train_dataset = TensorDataset(torch.FloatTensor(X_train), torch.FloatTensor(y_train))
        train_loader = DataLoader(train_dataset, batch_size=batch_size, shuffle=True)
        best_corr = -1
        best_corr_shuffle = -1
        best_celltype_corr = -1
        epochs_no_improve = 0
        best_weights = None
        for epoch in range(epochs):
            model.train()
            epoch_loss = 0
            for batch_X, batch_y in train_loader:
                optimizer.zero_grad()
                pred = model(batch_X).squeeze()

                loss = nn.PoissonNLLLoss(log_input=False)(pred, batch_y)
                # loss = nn.MSELoss()(pred,batch_y)

                if self.norm_type != "None":
                    loss += model.regularization()

                loss.backward()
                optimizer.step()

                with torch.no_grad():
                    model.linear.weight.data = torch.clamp(model.linear.weight.data, min=0)

                epoch_loss += loss.item()

            model.eval()
            with torch.no_grad():
                val_pred = model(torch.FloatTensor(X_test)).squeeze()
                val_pred_np = val_pred.numpy()

                current_corr,current_pv = stats.pearsonr(val_pred.numpy(), y_test)
                shuffle_corr,shuffle_pv = stats.pearsonr(np.random.permutation(val_pred.numpy()), y_test)

                current_celltype_corr = np.nan 
                current_celltype_pv = np.nan

                if celltype_test is not None:
                    df = pd.DataFrame({
                        'true': y_test,
                        'pred': val_pred_np,
                        'celltype': celltype_test
                    })
                    pseudo_bulk = df.groupby('celltype').mean()
                    if len(pseudo_bulk) > 1:
                        current_celltype_corr, current_celltype_pv = stats.pearsonr(pseudo_bulk['pred'], pseudo_bulk['true'])

            if current_corr > best_corr:
                best_corr = current_corr
                best_corr_pv = current_pv
                best_corr_shuffle = shuffle_corr
                best_corr_shuffle_pv = shuffle_pv

                if celltype_test is not None:
                    best_celltype_corr = current_celltype_corr
                    best_celltype_pv = current_celltype_pv

                best_weights = {
                    'weights': model.linear.weight.data.clone(),
                    'bias': model.linear.bias.data.clone()
                }
                epochs_no_improve = 0
            else:
                epochs_no_improve += 1
                if epochs_no_improve >= patience:
                    break
        if celltype is not None:
            result = {
                'test_corr': best_corr,
                'test_pv': best_corr_pv,
                'celltype_corr': best_celltype_corr,
                'celltype_pv': best_celltype_pv,
                'shuffle_corr': best_corr_shuffle,
                'shuffle_pv': best_corr_shuffle_pv,
                'weights': best_weights
            }

        else:
            result = {
                'test_corr': best_corr,
                'test_pv': best_corr_pv,
                'shuffle_corr': best_corr_shuffle,
                'shuffle_pv': best_corr_shuffle_pv,
                'weights': best_weights
            }
            
        return model,result


def set_gene_tile_significance_bootstrapped(x, y, w_mat, e, cell_info, celltype_col, clusters):
    """Compute standardized Shapley values for each group of cell clusters.

    adapt from SCARlink: https://www.nature.com/articles/s41588-024-01689-8
    
    Parameters
    ----------
    x : [[float]]
        Tile matrix for a gene.
    y : [float]
        Gene expression vector for same gene.
    w_mat : [float]
        Learned regression coefficients.
    e : float
        Learned regression bias.
    cell_info : data frame
        Cell metadata found in key cell_info in coassay_matrix.h5
    celltype_col : str
        Column in cell_info containing cell clusters.
    clusters : [str]
        List of clusters

    Returns
    -------
    zscore_d
        Dictionary of standardized z-scores for each tile and cell type
        for a given gene.
    """

    n = 5000
    clf = linear_model.PoissonRegressor()
    clf.fit(x, y)
    clf.coef_ = np.ravel(w_mat)
    clf.intercept_ = e
    ypred = clf.predict(x)
    explainer_linear = shap.LinearExplainer(clf, x) 
    pvals_d = {}
    shap_lst = []
    prod_lst = []
    b_clust = []
    w_mat = np.ravel(w_mat)
    shap_lst_c = []
    pv_lst = []
    for clust in clusters:
        clust_idx = (cell_info[celltype_col] == clust).values
        if np.sum(clust_idx) < 40: continue
        b_clust.append(clust)
        x_clust = x[clust_idx].copy()
        y_clust = y[clust_idx].copy()
        shap_lst_c = []
        ypred_subset = np.zeros((n, w_mat.shape[0]))
        ypred_minus_f = np.zeros((n, w_mat.shape[0]))
        for i in range(n):
            idx = np.random.choice(x_clust.shape[0], 30, replace=True)
            x_test_subset = x_clust[idx]
            shap_values_linear = explainer_linear.shap_values(x_test_subset)
            mean_shap = np.mean(shap_values_linear, axis=0)
            shap_lst_c.append(mean_shap)
        shap_lst.append(shap_lst_c)
    
    shap_lst = np.array(shap_lst)
    shap_lst = np.mean(shap_lst, axis=1)
    z_score = stats.zscore(np.ravel(shap_lst)).reshape(shap_lst.shape)
    zscore_d = dict(zip(b_clust, z_score))
    for c in clusters:
        if c not in zscore_d:
            zscore_d[c] = np.zeros(w_mat.shape[0])

    shap_lst_df = pd.DataFrame(shap_lst.T, columns=b_clust)
    return zscore_d, shap_lst_df