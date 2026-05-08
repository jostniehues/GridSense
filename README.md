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

[config]
nx_x = 22
nx_q2 = 15
x_min = 0.0001
x_max = 1.0
q2_min = 100.0
q2_max = 100000000.0
threshold = 0.88

[selection]
pdg = 0
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


If a name listed in `process_grid` does not exactly match a `[grids.<name>]` section, the run fails.

### `[config]`

This block controls the auxiliary support grid used for plotting.

- `nx_x`: number of bins in `x`
- `nx_q2`: number of bins in `Q^2`
- `x_min`, `x_max`: `x` range
- `q2_min`, `q2_max`: `Q^2` range
- `threshold`: relevance coverage threshold

The auxiliary grid is built in log-space and later shown on logarithmic axes.

### `[selection]`

This block controls global filters applied to all processed grids.

Currently implemented:

- `pdgs`: optional PDG filter. 
- `max_as`: maximum QCD order. 
If wanted the user can choose multiple PDGs in order to create a seperate plot for each pdg.

Example:

```toml
[selection]
pdg = [0,3]
max_as = 1
```

Notes:

- if `pdg` is omitted, no PDG filter is applied
- by convention, `pdg = 0` means gluon selection in the TOML. Internally, `GridSense` maps this to `0` or `21` depending on the grid convention
- electroweak order (if available) is kept up to the maximum available for the grid

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


