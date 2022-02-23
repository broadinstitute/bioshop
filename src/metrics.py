import numpy as np
from scipy.special import softmax
from sklearn.metrics import (
    auc,
    precision_recall_curve,
    roc_auc_score,
    f1_score,
    confusion_matrix,
    matthews_corrcoef,
)

def compute_two_class_metrics(name, pred_labels, pred_scores, labels):
    pred_scores = pred_scores[:, 1]
    roc_auc_pred_score = roc_auc_score(labels, pred_scores)
    precisions, recalls, thresholds = precision_recall_curve(labels, pred_scores)
    fscore = (2 * precisions * recalls) / (precisions + recalls)
    fscore[np.isnan(fscore)] = 0
    ix = np.argmax(fscore)
    threshold = thresholds[ix].item()
    pr_auc = auc(recalls, precisions)
    tn, fp, fn, tp = confusion_matrix(labels, pred_labels, labels=[0, 1]).ravel()
    result = {
        f'{name}_roc_auc': roc_auc_pred_score,
        f'{name}_threshold': threshold,
        f'{name}_pr_auc': pr_auc,
        f'{name}_recall': recalls[ix].item(),
        f'{name}_precision': precisions[ix].item(),
        f'{name}_f1': fscore[ix].item(),
        f'{name}_tn': tn.item(),
        f'{name}_fp': fp.item(),
        f'{name}_fn': fn.item(), 
        f'{name}_tp': tp.item()
    }
    return result

def compute_metrics_by_class(logits, labels):
    logits = logits.reshape(logits.shape[0], 2, 2)
    pred_labels = np.argmax(logits, axis=-1)
    pred_scores = softmax(logits, axis=-1)
    snp_mask = (labels < 2)
    # snp
    pred_scores_snp = pred_scores[:, 0][snp_mask]
    pred_labels_snp = pred_labels[:, 0][snp_mask]
    labels_snp = labels[snp_mask]
    snp_metrics = compute_two_class_metrics('snp', pred_labels_snp, pred_scores_snp, labels_snp)
    # indel
    pred_scores_indel = pred_scores[:, 1][~snp_mask]
    pred_labels_indel = pred_labels[:, 1][~snp_mask]
    labels_indel = labels[~snp_mask] - 2
    indel_metrics = compute_two_class_metrics('indel', pred_labels_indel, pred_scores_indel, labels_indel)
    metrics = {}
    metrics.update(snp_metrics)
    metrics.update(indel_metrics)
    return metrics

def compute_multi_metrics(logits, labels):
    pred_labels = np.argmax(logits, axis=-1)
    result = {
        'multi_acc': (pred_labels == labels).mean(),
        'multi_f1': f1_score(y_true=labels, y_pred=pred_labels, average='micro')
    }
    return result

def compute_metrics(pred):
    logits = pred.predictions[0]
    labels = np.squeeze(pred.label_ids)
    metrics = {}
    metrics.update(compute_multi_metrics(logits, labels))
    metrics.update(compute_metrics_by_class(logits, labels))
    return metrics
