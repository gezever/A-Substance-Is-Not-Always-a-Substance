# GPU & CUDA Setup for sentence-transformers

## Hardware

| | |
|---|---|
| GPU | NVIDIA GeForce RTX 3060 (GA106 Lite Hash Rate) |
| VRAM | 12 288 MiB |
| PCI slot | 05:00.0 |
| OS | Linux Mint 22.2 |
| Kernel | 6.17.0-20-generic |

## Status (April 2026) — fully operational

| Component | Version / state |
|---|---|
| NVIDIA driver | 580.126.09 |
| CUDA (driver level) | 13.0 |
| PyTorch | 2.11.0+cu130 |
| sentence-transformers | 5.4.0 |
| `torch.cuda.is_available()` | `True` |
| `/dev/nvidia*` devices | present |
| `nvidia-smi` | available |

The GPU is fully usable. The R scripts (`08`, `09`, `09b`, `09c`) call
`model$encode()` via sentence-transformers in the reticulate uv environment;
CUDA is used automatically without any code changes.

---

## Reticulate uv environment

reticulate manages its own Python environment via `uv`. The active venv with
CUDA support is located at:

```
~/.cache/R/reticulate/uv/cache/archive-v0/Ki1QH51N4jeG6c5pJ9hTD/
```

Installing or upgrading packages in this environment:

```bash
~/.cache/R/reticulate/uv/bin/uv pip install <package> \
  --python ~/.cache/R/reticulate/uv/cache/archive-v0/Ki1QH51N4jeG6c5pJ9hTD/bin/python
```

Or activate the venv temporarily in a shell:

```bash
source ~/.cache/R/reticulate/uv/cache/archive-v0/Ki1QH51N4jeG6c5pJ9hTD/bin/activate
pip show torch          # shows Version: 2.11.0, Location: ...Ki1QH51N4j.../site-packages
deactivate
```

---

## Verification in R

```r
library(reticulate)
torch <- reticulate::import("torch")
torch$cuda$is_available()        # TRUE
torch$cuda$get_device_name(0L)   # "NVIDIA GeForce RTX 3060"
torch$version$cuda               # "13.0"
```

Expected `nvidia-smi` output (abbreviated):

```
+-----------------------------------------------------------------------------------------+
| NVIDIA-SMI 580.126.09  Driver Version: 580.126.09   CUDA Version: 13.0               |
+---------------------------+------------------------+-----------------------------------+
|   0  NVIDIA GeForce RTX 3060  Off | 00000000:05:00.0 On |                         N/A |
|  0%  44C  P8  16W / 170W |  617MiB / 12288MiB |    0%    Default                     |
+---------------------------+------------------------+-----------------------------------+
```

---

## Effect on the scripts

The embedding step in `08_embedding_clustering.R`, `09_embedding_chemont.R`,
`09b_chemont_model_comparison.R`, and `09c_chemont_model_validation.R` calls
`model$encode()` via sentence-transformers. When CUDA is available,
sentence-transformers automatically selects the GPU — no code changes required.

The scibert and biobert models in `09b` / `09c` are larger than miniLM and
benefit most from GPU acceleration. Cached embeddings
(`.rds` files in `data/cache/embeddings/`) do not need to be recomputed
after the driver installation.

---

## How the setup was established

### Step 1 — Install the NVIDIA driver

```bash
sudo apt install nvidia-driver-580-open
sudo reboot
```

After the reboot, `/dev/nvidia*` devices are present and `nvidia-smi` works.

### Step 2 — Install PyTorch with CUDA

The reticulate uv environment installs the CPU variant of PyTorch by default.
The CUDA variant was installed via `uv pip` using the PyTorch CUDA 13.0
index URL:

```bash
~/.cache/R/reticulate/uv/bin/uv pip install torch \
  --index-url https://download.pytorch.org/whl/cu130 \
  --python ~/.cache/R/reticulate/uv/cache/archive-v0/Ki1QH51N4jeG6c5pJ9hTD/bin/python
```

Adjust the index URL if a future driver reports a different CUDA version
(`nvidia-smi` shows the supported version in the top-right corner).

| CUDA version (nvidia-smi) | index URL |
|---|---|
| 12.1 | `https://download.pytorch.org/whl/cu121` |
| 12.4 | `https://download.pytorch.org/whl/cu124` |
| 13.0 | `https://download.pytorch.org/whl/cu130` |

---

## Available drivers (reference)

```
nvidia-driver-580-open        ← installed (recommended)
nvidia-driver-580
nvidia-driver-590-open
nvidia-driver-590
nvidia-driver-570-open
nvidia-driver-570
nvidia-driver-535-open
nvidia-driver-535
nvidia-driver-470
xserver-xorg-video-nouveau    ← open source, no CUDA
```
