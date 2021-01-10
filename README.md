# pyP
`pyP` is a Python GUI tool to pick arrival time of *P* phase

## Usage
### On terminal
- Execute `pyP.py`
```zsh
python pyP.py 7 "work/2020-01-24-mww67-turkey/*.SAC" 38.392 39.085
```

- positional arguments:
  - `displayNum`: Number of traces shown in display (e.g., 7)
  - `sacfiles`: SAC files you want to pick P arrival (e.g., "./*.SAC"). comma-separated list is available. *Do not forget quotation marks!
  - `elat`: Latitude of epicentre (for station azimuth)
  - `elon`: Longitude of epicentre (for station azimuth)

- optional arguments:
  - `-h`, `--help`: show this help message and exit


### On GUI
- Click to pick arrival time

- Control panel by pressing keys:
  - `a`: Zoom in xlim
  - `z`: Zoom out xlim
  - `x`: Reset xlim
  - `.`: Zoom in ylim
  - `,`: Zoom out ylim
  - `↓`: Display next trace
  - `→`: Display further next traces
  - `↑`: Display previous trace
  - `←`: Display further previous traces

- Click `Save` button to save picked time and overwrite `a` marker in SAC header

![](./image/screen.png)
