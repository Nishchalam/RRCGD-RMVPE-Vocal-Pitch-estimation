import numpy as np
import torch
from collections import defaultdict

from src import to_local_average_cents, bce, PITCH_TH
from mir_eval.melody import raw_pitch_accuracy, to_cent_voicing, raw_chroma_accuracy, overall_accuracy
from mir_eval.melody import voicing_recall, voicing_false_alarm


def evaluate_batch(pitch_pred, pitch_label, hop_length):
    # Model outputs logits — convert to probabilities before metric computation.
    pitch_pred  = torch.sigmoid(pitch_pred).detach().cpu().numpy()
    pitch_label = pitch_label.detach().cpu().numpy()

    rpa_list, oa_list, vr_list = [], [], []

    for pred, label in zip(pitch_pred, pitch_label):
        cents_pred  = to_local_average_cents(pred,  None, 0.0)
        cents_label = to_local_average_cents(label, None, 0.0)

        freq_pred = np.array([10 * (2 ** (c / 1200)) if c else 0 for c in cents_pred])
        freq      = np.array([10 * (2 ** (c / 1200)) if c else 0 for c in cents_label])

        time_slice = np.array([i * hop_length / 1000 for i in range(len(cents_label))])

        ref_v, ref_c, est_v, est_c = to_cent_voicing(time_slice, freq, time_slice, freq_pred)

        rpa_list.append(raw_pitch_accuracy(ref_v, ref_c, est_v, est_c))
        oa_list.append(overall_accuracy(ref_v, ref_c, est_v, est_c))
        vr_list.append(voicing_recall(ref_v, est_v))

    return {
        "RPA": np.mean(rpa_list),
        "OA":  np.mean(oa_list),
        "VR":  np.mean(vr_list),
    }


def evaluate(dataset, model, hop_length, device, pitch_th=PITCH_TH):
    metrics = defaultdict(list)

    for data in dataset:
        roots       = data["roots"].to(device)
        pitch_label = data["pitch"].to(device)

        if roots.dim() == 2:
            roots = roots.unsqueeze(0)

        with torch.no_grad():
            logits     = model(roots).squeeze(0)           # raw logits
            pitch_pred = torch.sigmoid(logits)             # probabilities for metrics

        T          = pitch_label.shape[0]
        pitch_pred = pitch_pred[:T]
        logits     = logits[:T]

        # bce() now takes logits internally (BCEWithLogitsLoss)
        loss = bce(logits, pitch_label)
        metrics["loss"].append(loss.item())

        cents_pred  = to_local_average_cents(pitch_pred.detach().cpu().numpy(), None, pitch_th)
        cents_label = to_local_average_cents(pitch_label.detach().cpu().numpy(), None, pitch_th)

        freq_pred = np.array([10 * (2 ** (c / 1200)) if c else 0 for c in cents_pred])
        freq      = np.array([10 * (2 ** (c / 1200)) if c else 0 for c in cents_label])

        time_slice = np.array([i * hop_length / 1000 for i in range(len(cents_label))])

        ref_v, ref_c, est_v, est_c = to_cent_voicing(time_slice, freq, time_slice, freq_pred)

        metrics["RPA"].append(raw_pitch_accuracy(ref_v, ref_c, est_v, est_c))
        metrics["RCA"].append(raw_chroma_accuracy(ref_v, ref_c, est_v, est_c))
        metrics["OA"].append(overall_accuracy(ref_v, ref_c, est_v, est_c))
        metrics["VFA"].append(voicing_false_alarm(ref_v, est_v))
        metrics["VR"].append(voicing_recall(ref_v, est_v))

    return metrics
