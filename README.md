# GridSense

`GridSense` builds kinematic support plots from one or more `PineAPPL` interpolation grids. Grids from other interpolation libraries must be imported via PineAPPL CLI before they can be used.

The current workflow is:

1. read one or more grids from `grids.toml`
2. extract kinematic contributions in Rust
3. project and filter them in Rust
4. build an auxiliary `(x, Q^2)` grid
5. rank cells by relevance and keep the requested coverage
6. pass the resulting support masks to Python
7. save a PDF plot

Rust handles the data pipeline. Python is only used for plotting.

## Requirements

You currently need:

- Rust with `cargo`
- Python 3
- `numpy`
- `matplotlib`
- LHAPDF installed and usable through the Rust `lhapdf` crate
- access to the local `pineappl` fork referenced in `Cargo.toml`

The project currently depends on a local PineAPPL checkout:

```toml
pineappl = { path = "/home/jnieh/uni/projects/pineappl_fork/pineappl/pineappl" }
```

If your PineAPPL fork lives somewhere else, adjust that path in `Cargo.toml`.

## Python setup

The plot is produced by:

```text
python/plot_auxgrid.py
```

`src/main.rs` currently starts it with:

```rust
Command::new("python3")
```

That means:

- `python3` must be available on your `PATH`
- `numpy` and `matplotlib` must be installed in that Python environment

If your Python executable is not called `python3`, you currently have to change the command in [src/main.rs](/home/jnieh/uni/projects/gridsense/src/main.rs).

Rust writes a temporary JSON file and passes its path plus the final PDF path to the plotting script. That JSON file is created with `tempfile::NamedTempFile` and is not intended as a permanent output.

## Running

The project is configured through `grids.toml`.

Run it with:

```bash
cargo run -- --config grids.toml
```

## `grids.toml`

The main configuration file is `grids.toml`.

### Example

```toml
path_to_grids = "/home/user/grids/pineappl"

[grids.wm_atlas]
path = "WZ/matrixhawaii-atlas-wm/.../grid.lz4"
label = 'Wm ATLAS $\frac{d\sigma}{|d\eta_l|}$'

[grids.z0_lhcb]
path = "WZ/matrixhawaii-lhcb-z0/.../grid.lz4"
label = 'Z0 LHCb $\frac{d\sigma}{dy_{\ell\ell}}$'

[output]
process_grid = ["wm_atlas", "z0_lhcb"]
path_plot = "output/plots/fullgrid.pdf"
detail = false

[output.groups]
"Vector bosons" = ["wm_atlas", "z0_lhcb"]

[config]
nx_x = 22
nx_q2 = 15
x_min = 0.0001
x_max = 1.0
q2_min = 100.0
q2_max = 100000000.0
threshold = 0.88

[selection_conf]
pdf = "CT18NNLO"
channels = ["6-17"]
leg_pdgs = [0]
max_as = 1
```

### `path_to_grids`

This is a base directory used together with relative grid paths.

### `[grids.<name>]`

Each grid entry defines one named grid.

Fields:

- `path`: path to the `.lz4` file
- `label`: legend label used in the plot

The section name `wm_atlas` is the internal grid name. You later reference that name in `output.process_grid`.

### `[output]`

This block controls what gets plotted and where the final PDF is written.

- `process_grid`: list of grid names to process
- `path_plot`: output path for the generated PDF
- `detail`: (bool). If a .csv should be generated per grid writing (bin,q2,x,contrib) for each        interpolation grid node.
- `output.groups`: optional named groups. Group members share one color and line style and are plotted as one combined area with a common outer boundary. Every grouped grid must also occur in `process_grid`.


If a name listed in `process_grid` does not exactly match a `[grids.<name>]` section, the run fails.

### `[config]`

This block controls the auxiliary support grid used for plotting.

- `nx_x`: number of bins in `x`
- `nx_q2`: number of bins in `Q^2`
- `x_min`, `x_max`: `x` range
- `q2_min`, `q2_max`: `Q^2` range
- `threshold`: relevance coverage threshold

The auxiliary grid is built in log-space and later shown on logarithmic axes.

### Grid projection selection: `[selection_conf]`

This block selects the PDF, perturbative orders, partonic channels, and incoming legs used for all processed grids.

- `pdf`: LHAPDF set used for the convolution
- `max_as`: optional maximum power of the strong coupling
- `pdgs`: optional initial-state PDG selection
- `channels`: optional list of zero-based PineAPPL channel indices
- `leg_pdgs`: optional list of PDGs whose matching incoming legs are projected

There are three projection modes.

#### Initial-state PDG selection

```toml
[selection_conf]
pdf = "CT18NNLO"
pdgs = [0, 3]
```

Each value in `pdgs` creates a separate projection. A channel entry is retained when the requested PDG occurs anywhere in its initial state, and both incoming legs are projected. Thus this mode describes the support of processes containing the requested parton, rather than only the location of that parton.

#### Channel selection

```toml
[selection_conf]
pdf = "CT18NNLO"
channels = [6, "7-17", 21]
```

`channels` selects complete PineAPPL channels. Integer entries select one channel; quoted ranges are inclusive. Duplicate indices are removed and the resulting indices are sorted. Both incoming legs are projected.

Channel indices are zero-based and refer to the ordering in the individual PineAPPL grid. They can be inspected with:

```bash
pineappl read --channels path/to/grid.lz4
```

#### Leg-specific PDG projection

```toml
[selection_conf]
pdf = "CT18NNLO"
channels = ["6-17"]
leg_pdgs = [0]
```

`leg_pdgs` groups all listed PDGs into one projection and retains only matching incoming legs. For a gluon selection, for example:

```text
[g, q] -> project x1
[q, g] -> project x2
[g, g] -> project x1 and x2
[q, q] -> discard the entry
```

`leg_pdgs` may also be used without `channels`, in which case matching legs are selected across all channels. PDG matching is exact, so a combined quark and antiquark projection must list both signs:

```toml
leg_pdgs = [1, 2, 3, 4, 5, -1, -2, -3, -4, -5]
```

The allowed combinations are:

```text
pdgs                     allowed
channels                 allowed
leg_pdgs                 allowed
channels + leg_pdgs      allowed
pdgs + channels          rejected
pdgs + leg_pdgs          rejected
```

If none of `pdgs`, `channels`, or `leg_pdgs` is given, all channel entries and both incoming legs are used. By convention, PDG `0` denotes the gluon in the TOML; GridSense maps it to `0` or `21` according to the grid convention. Electroweak orders are retained up to the maximum available for the grid.

## Labels and LaTeX-style text

Grid labels are read from `grids.toml` and passed to Matplotlib. You can use Mathtext-style expressions.

## Output

The project may currently create:

- a CSV file with projected contributions:
  - `output/data/grid.csv`
- the final PDF plot:
  - configured through `output.path_plot`

The CSV is currently a fixed debug-style output and may change later.

The temporary JSON sent to Python is not meant to be kept.
