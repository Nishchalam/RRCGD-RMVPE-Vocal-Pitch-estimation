import numpy as np
import torch.nn.functional as F
import torch.nn as nn
import librosa
import torch
import os

from .constants import N_BINS, SAMPLE_RATE, F_MAX, BIN_METHOD, N_DENSE

import matplotlib.pyplot as plt


class RRCGDSpectrogram(nn.Module):
    """
    ACS of RRCGD spectrogram — optimised variant.

    Frequency representation modes:

    Base methods:
        linear_fmax
            Generate n_bins linearly spaced frequencies from 0 -> f_max.

        mel_aligned
            Generate n_bins Mel-aligned center frequencies from 0 -> f_max.

    Truncation methods:
        linear_fmax_upper_frequency_truncated
        mel_aligned_upper_frequency_truncated

            First generate the ORIGINAL n_bins frequency grid over
            0 -> fs/2. Then retain only bins whose center frequency
            is <= f_max.

            The retained bins are zero-padded to the smallest multiple
            of 2^5 = 32 required by the current DeepUNet architecture.

        linear_fmax_upper_bins_truncated
        mel_aligned_upper_bins_truncated

            First generate the ORIGINAL n_bins frequency grid over
            0 -> fs/2. Then retain the first n_bins // 2 bins.

            This is an INDEX-BASED truncation and does not use the
            physical frequency value for selecting bins.

    Important:
        f_max has two meanings depending on the mode:

        Full modes:
            f_max = actual upper frequency of the generated grid.

        Frequency-truncation modes:
            f_max = physical frequency cutoff.

        Upper-bin-truncation modes:
            f_max identifies the experiment/cutoff but does not
            determine bin selection; selection is by index.

    The original frequency grid for truncation modes ALWAYS spans
    0 -> fs/2 before truncation.
    """

    def __init__(
        self,
        n_bins=N_BINS,
        r0=1.008,
        hop_length=320,
        device="cuda",
        save_dir=None,
        bin_method=BIN_METHOD,
        f_max=None,
        chunk_size=64,
    ):
        super().__init__()

        self.r0 = r0
        self.r0_sq = r0 ** 2

        self.hop_length = hop_length
        self.step = 0
        self.save_dir = save_dir

        self.bin_method = bin_method
        self.chunk_size = chunk_size
        self.pool_k = None

        sr = SAMPLE_RATE

        # Preserve existing behaviour when f_max is not explicitly supplied.
        f_max = F_MAX if f_max is None else float(f_max)

        # Store the ACTUAL experiment f_max.
        # This is important because truncation experiments may use fs/4
        # while the original baseline may use fs/2.
        self.f_max = f_max

        # ------------------------------------------------------------
        # Identify base frequency representation and truncation mode
        # ------------------------------------------------------------

        truncation_suffixes = {
            "_upper_frequency_truncated": "freq_truncate",
            "_upper_bins_truncated": "bin_truncate",
        }

        base_method = bin_method
        mode = "full"

        for suffix, trunc_mode in truncation_suffixes.items():
            if bin_method.endswith(suffix):
                base_method = bin_method[: -len(suffix)]
                mode = trunc_mode
                break

        # ------------------------------------------------------------
        # Linear / Mel frequency representations
        # ------------------------------------------------------------

        if base_method in ("linear_fmax", "mel_aligned"):

            # --------------------------------------------------------
            # FULL REPRESENTATIONS
            #
            # E1:
            #     f_max = fs/2
            #
            # E2:
            #     f_max = fs/4
            #
            # The grid is directly generated over 0 -> f_max.
            #
            # --------------------------------------------------------
            if mode == "full":
                grid_fmax = f_max

            # --------------------------------------------------------
            # TRUNCATION REPRESENTATIONS
            #
            # E3/E4:
            #
            # First generate the ORIGINAL grid over 0 -> fs/2.
            # Only afterwards perform frequency/index truncation.
            #
            # --------------------------------------------------------
            else:
                grid_fmax = sr / 2.0

            # --------------------------------------------------------
            # Construct the frequency grid
            # --------------------------------------------------------

            if base_method == "mel_aligned":

                mel_hz = librosa.mel_frequencies(
                    n_mels=n_bins + 2,
                    fmin=0.0,
                    fmax=grid_fmax,
                )

                # Preserve the existing centre-frequency convention.
                centers_hz = mel_hz[1:-1]

            else:
                # linear_fmax

                max_angle = (
                    grid_fmax / (sr / 2.0)
                ) * np.pi

                w_full = np.linspace(
                    0.0,
                    max_angle,
                    n_bins,
                )

                centers_hz = (
                    w_full / (2.0 * np.pi)
                ) * sr

            # --------------------------------------------------------
            # Select bins
            # --------------------------------------------------------

            if mode == "full":

                # E1 / E2:
                # Keep the entire generated grid.
                keep_idx = np.arange(n_bins)

            elif mode == "freq_truncate":

                # E3:
                #
                # The source grid was generated over fs/2.
                # f_max is now the physical frequency cutoff.
                #
                # IMPORTANT:
                # Do not regenerate a new grid over f_max.
                keep_idx = np.where(
                    centers_hz <= self.f_max
                )[0]

            else:
                # E4:
                #
                # Index-based truncation.
                # Ignore the physical frequency spacing.
                keep_idx = np.arange(
                    n_bins // 2
                )

            # --------------------------------------------------------
            # Convert selected frequencies to angular frequencies
            # --------------------------------------------------------

            centers_kept = centers_hz[keep_idx]

            w_np = (
                2.0
                * np.pi
                * centers_kept
                / sr
            )

            # Number of actual frequency bins before padding.
            raw_n = len(w_np)

            self.raw_n_bins = raw_n

            # --------------------------------------------------------
            # DeepUNet compatibility
            #
            # Current architecture:
            #
            #     en_de_layers = 5
            #     2^5 = 32
            #
            # Therefore the frequency dimension must be divisible by 32.
            #
            # Only truncation modes need padding for the current
            # experiments. Padding contains ZERO values and therefore
            # does not introduce additional frequency content.
            # --------------------------------------------------------

            if mode != "full":

                model_divisor = 2 ** 5  # current DeepUNet: en_de_layers=5

                self.pad_amount = (
                    int(
                        np.ceil(
                            raw_n / model_divisor
                        )
                        * model_divisor
                    )
                    - raw_n
                )

            else:

                self.pad_amount = 0

            self.n_bins = raw_n + self.pad_amount

        # ------------------------------------------------------------
        # Existing maxpool mode
        # ------------------------------------------------------------

        elif bin_method == "maxpool":

            assert N_DENSE % n_bins == 0, (
                f"N_DENSE ({N_DENSE}) must be divisible "
                f"by N_BINS ({n_bins})"
            )

            max_angle = (
                f_max / (sr / 2.0)
            ) * np.pi

            w_np = np.linspace(
                0.0,
                max_angle,
                N_DENSE,
            )

            self.pool_k = N_DENSE // n_bins

            self.n_bins = n_bins
            self.pad_amount = 0
            self.raw_n_bins = n_bins

        # ------------------------------------------------------------
        # Existing Mel maxpool mode
        # ------------------------------------------------------------

        elif bin_method == "mel_maxpool":

            assert N_DENSE % n_bins == 0, (
                f"N_DENSE ({N_DENSE}) must be divisible "
                f"by N_BINS ({n_bins})"
            )

            mel_hz_fine = librosa.mel_frequencies(
                n_mels=N_DENSE + 2,
                fmin=0.0,
                fmax=f_max,
            )

            center_freqs_hz = mel_hz_fine[1:-1]

            w_np = (
                2.0
                * np.pi
                * center_freqs_hz
                / sr
            )

            self.pool_k = N_DENSE // n_bins

            self.n_bins = n_bins
            self.pad_amount = 0
            self.raw_n_bins = n_bins

        # ------------------------------------------------------------
        # Existing adaptive Mel mode
        # ------------------------------------------------------------

        elif bin_method == "mel_adaptive":

            mel_hz_fine = librosa.mel_frequencies(
                n_mels=N_DENSE + 2,
                fmin=0.0,
                fmax=f_max,
            )

            center_freqs_hz = mel_hz_fine[1:-1]

            w_np = (
                2.0
                * np.pi
                * center_freqs_hz
                / sr
            )

            self.n_bins = n_bins
            self.pad_amount = 0
            self.raw_n_bins = n_bins

        else:

            raise ValueError(
                f"Unknown bin_method '{bin_method}'. "
                "Choose 'mel_aligned', 'linear_fmax', "
                "'maxpool', 'mel_maxpool', 'mel_adaptive', "
                "or one of those with "
                "'_upper_frequency_truncated' / "
                "'_upper_bins_truncated' appended."
            )

        # ------------------------------------------------------------
        # Store frequency grid
        # ------------------------------------------------------------

        w = torch.tensor(
            w_np,
            dtype=torch.float32,
            device=device,
        )

        self.register_buffer(
            "w_grid",
            w.view(1, -1),
        )

    # ========================================================================
    # Internal helpers
    # ========================================================================

    def _chunked_rrcgd_sum(
        self,
        m: torch.Tensor,
        theta: torch.Tensor,
        BT: int,
        K: int,
    ) -> torch.Tensor:

        """
        Accumulate RRCGD over roots in chunks.

        m, theta:
            (BT, K)

        Returns:
            (BT, number_of_frequency_samples)
        """

        w_size = self.w_grid.shape[-1]

        out = torch.zeros(
            BT,
            w_size,
            device=self.w_grid.device,
            dtype=torch.float32,
        )

        for k0 in range(
            0,
            K,
            self.chunk_size,
        ):

            k1 = min(
                k0 + self.chunk_size,
                K,
            )

            m_c = (
                m[:, k0:k1]
                .unsqueeze(-1)
            )

            th_c = (
                theta[:, k0:k1]
                .unsqueeze(-1)
            )

            cos_t = torch.cos(
                self.w_grid.unsqueeze(0)
                - th_c
            )

            num = (
                m_c ** 2
                - m_c * self.r0 * cos_t
            )

            den = (
                self.r0_sq
                + m_c ** 2
                - 2.0 * m_c * self.r0 * cos_t
            )

            # den is bounded below by (r0 - 1)^2 for m in (0, 1], but that
            # floor (~6.4e-5 for r0=1.008) is still small enough that a
            # root near-resonant with a grid frequency produces a large
            # num/den spike. This feature feeds straight into the DeepUnet's
            # Conv2d/BatchNorm2d stack, which autocast runs in fp16 — an
            # unclamped spike there overflows fp16 and permanently corrupts
            # BatchNorm running stats (nan loss that never recovers, since
            # those stat updates aren't protected by GradScaler). Clamp each
            # chunk's contribution to keep the feature in a safe range.
            contrib = num / (den + 1e-9)
            if os.environ.get("EXP_RRCGD_CLAMP", "1") != "0":
                contrib = torch.clamp(contrib, -50.0, 50.0)

            out += contrib.sum(dim=1)

        return torch.nan_to_num(out, nan=0.0, posinf=50.0, neginf=-50.0)

    def _roots_to_m_theta(self, roots):

        """
        Unpack complex roots into clipped magnitude m and angle theta.
        """

        if (
            roots.dim() == 4
            and roots.shape[-1] == 2
        ):
            roots = torch.view_as_complex(
                roots.contiguous()
            )

        if roots.dim() == 2:
            roots = roots.unsqueeze(0)

        roots = (
            roots
            .to(self.w_grid.device)
            .to(torch.complex64)
        )

        B, T, K = roots.shape

        mag = torch.abs(roots)

        m = torch.where(
            mag > 1.0,
            1.0 / (mag + 1e-9),
            mag,
        )

        return (
            m.view(B * T, K),
            torch.angle(roots).view(B * T, K),
            B,
            T,
            K,
        )

    # ========================================================================
    # Forward
    # ========================================================================

    def forward(self, roots):

        """
        roots:
            (B,T,K)
            (T,K)
            (B,T,K,2)

        Returns:
            (B, effective_n_bins, T)
        """

        m, theta, B, T, K = (
            self._roots_to_m_theta(roots)
        )

        # ------------------------------------------------------------
        # Step 1:
        # Compute RRCGD over the selected frequency grid.
        #
        # For truncation modes this is already the truncated grid.
        # ------------------------------------------------------------

        rrcgd = self._chunked_rrcgd_sum(
            m,
            theta,
            B * T,
            K,
        )

        # ------------------------------------------------------------
        # Step 2:
        # Existing dense-grid pooling modes
        # ------------------------------------------------------------

        if self.bin_method in (
            "maxpool",
            "mel_maxpool",
        ):

            rrcgd = F.max_pool1d(
                rrcgd.unsqueeze(1),
                kernel_size=self.pool_k,
                stride=self.pool_k,
            ).squeeze(1)

        elif self.bin_method == "mel_adaptive":

            rrcgd = F.adaptive_max_pool1d(
                rrcgd.unsqueeze(1),
                output_size=self.n_bins,
            ).squeeze(1)

        # ------------------------------------------------------------
        # Step 3:
        # Zero-pad truncated representations.
        #
        # IMPORTANT:
        # Padding happens AFTER frequency selection.
        #
        # Therefore these zeros do NOT represent additional
        # frequency bins/content.
        # ------------------------------------------------------------

        elif getattr(self, "pad_amount", 0) > 0:
            # Checkpoints pickled before pad_amount existed (all "full"-mode
            # E1 runs) never set this attribute; "full" mode always resolves
            # it to 0 in current __init__, so defaulting missing to 0 here
            # is exactly consistent with what a fresh instance would have.

            rrcgd = F.pad(
                rrcgd,
                (0, self.pad_amount),
            )

        # ------------------------------------------------------------
        # Step 4:
        # Reshape:
        #
        #     (BT, effective_n_bins)
        #
        # ->  (B, T, effective_n_bins)
        #
        # ->  (B, effective_n_bins, T)
        # ------------------------------------------------------------

        out = (
            rrcgd
            .view(
                B,
                T,
                self.n_bins,
            )
            .transpose(1, 2)
        )

        # ------------------------------------------------------------
        # Optional visualization
        # ------------------------------------------------------------

        if (
            self.save_dir is not None
            and self.step % 1000 == 0
            and torch.cuda.current_device() == 0
        ):

            spec = (
                out[0]
                .detach()
                .cpu()
                .numpy()
            )

            os.makedirs(
                self.save_dir,
                exist_ok=True,
            )

            fig, ax = plt.subplots(
                figsize=(14, 4)
            )

            im = ax.imshow(
                spec,
                aspect="auto",
                origin="lower",
                cmap="magma",
            )

            fig.colorbar(
                im,
                ax=ax,
            )

            ax.set_title(
                f"ACS-RRCGD | "
                f"{self.bin_method} | "
                f"F_MAX={self.f_max} Hz | "
                f"step {self.step}"
            )

            fig.savefig(
                os.path.join(
                    self.save_dir,
                    f"acs_rrcgd_"
                    f"{self.bin_method}_"
                    f"step_{self.step}.png",
                ),
                dpi=150,
                bbox_inches="tight",
            )

            plt.close(fig)

        self.step += 1

        return out.float()