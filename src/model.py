import torch
from torch import nn
from .deepunet import DeepUnet
from .constants import *
from .spec import RRCGDSpectrogram
from .seq import BiGRU
import torch.nn.functional as F


class E2E(nn.Module):
    def __init__(
        self,
        hop_length,
        n_blocks,
        n_gru,
        kernel_size,
        en_de_layers=5,
        inter_layers=4,
        in_channels=1,
        en_out_channels=16,
        device="cuda",
        writer=None
    ):
        super().__init__()

        # ------------------------------------------------------------
        # RRCGD feature extractor
        # ------------------------------------------------------------
        self.rrcgd = RRCGDSpectrogram(
            n_bins=N_BINS,
            r0=R0,
            hop_length=hop_length,
            device=device,
            save_dir=None,
            bin_method=BIN_METHOD,
            f_max=F_MAX,
        )

        # Effective frequency dimension after:
        #   - full representation
        #   - frequency truncation
        #   - bin truncation
        #   - model-compatible zero padding
        n_bins_eff = self.rrcgd.n_bins

        # ------------------------------------------------------------
        # DeepUNet
        # ------------------------------------------------------------
        self.unet = DeepUnet(
            kernel_size,
            n_blocks,
            en_de_layers,
            inter_layers,
            in_channels,
            en_out_channels,
            n_bins=n_bins_eff,
        )

        # ------------------------------------------------------------
        # CNN projection
        # ------------------------------------------------------------
        self.cnn = nn.Conv2d(
            en_out_channels,
            3,
            (3, 3),
            padding=(1, 1)
        )

        # ------------------------------------------------------------
        # Sequence / classification head
        # ------------------------------------------------------------
        if n_gru:
            self.fc = nn.Sequential(
                BiGRU(
                    3 * n_bins_eff,
                    256,
                    n_gru
                ),
                nn.Linear(
                    512,
                    N_CLASS
                ),
                nn.Dropout(0.25),

                # Model outputs raw logits.
                # Use BCEWithLogitsLoss / binary_cross_entropy_with_logits
                # during training and apply sigmoid explicitly when
                # probabilities are required.
            )

        else:
            self.fc = nn.Sequential(
                nn.Linear(
                    3 * n_bins_eff,
                    N_CLASS
                ),
                nn.Dropout(0.25),

                # Model outputs raw logits.
            )

    def forward(self, x):

        # ------------------------------------------------------------
        # RRCGD
        #
        # (B, N_BINS_eff, T)
        # ------------------------------------------------------------
        rrcgd = self.rrcgd(x).unsqueeze(1)

        # ------------------------------------------------------------
        # Temporal padding for DeepUNet
        # ------------------------------------------------------------
        T_orig = rrcgd.shape[-1]

        pad = (
            32 - (T_orig % 32)
        ) % 32

        if pad > 0:
            rrcgd = F.pad(
                rrcgd,
                (0, pad)
            )

        # ------------------------------------------------------------
        # DeepUNet + CNN
        # ------------------------------------------------------------
        x = self.cnn(
            self.unet(rrcgd)
        )

        # Remove only the temporal padding.
        x = x[:, :, :, :T_orig]

        # ------------------------------------------------------------
        # Rearrange:
        #
        # (B, 3, N_BINS_eff, T)
        #
        # -> (B, T, 3, N_BINS_eff)
        #
        # -> (B, T, 3*N_BINS_eff)
        # ------------------------------------------------------------
        x = x.permute(
            0,
            3,
            1,
            2
        ).flatten(-2)

        # ------------------------------------------------------------
        # Classification head
        # ------------------------------------------------------------
        for layer in self.fc:
            x = layer(x)

        return x