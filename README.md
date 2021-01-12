# pyP
`pyP` is a Python GUI to pick arrival time of seismic phase

![](./image/screen.gif)

## Usage
- Execute `pyP.py` on termial
```zsh
python pyP.py 7 "work/*.SAC" 38.392 39.085
```

- Positional arguments:
  - `displayNum`: number of traces shown at once in display (e.g., 7)
  - `sacfiles`: SAC files you want to pick (e.g., "./*.SAC")
    - comma-separated list is available (e.g., "./II*.SAC, ./IU.SAC")
    - *Do not forget quotation marks!*
  - `elat`: latitude of epicentre (degree)
  - `elon`: longitude of epicentre (degree)
    - epicentre is used for calculating station azimuth and distance

- Optional arguments:
  - `-h`, `--help`: show this help message and exit

- Control panel by pressing keys:
  - `a`: zoom in xlim
  - `z`: zoom out xlim
  - `x`: reset xlim
  - `.`: zoom in ylim
  - `,`: zoom out ylim
  - `↓`: frame advance
  - `↑`: frame back
  - `→`: page advance
  - `←`: page back

- Click trace to pick arrival time

- Click `Save` button to store arrival time as `a` marker in SAC header

## Dependencies

- [matplotlib](https://matplotlib.org/)
- [ObsPy](https://docs.obspy.org/)
- [geographiclib](https://pypi.org/project/geographiclib/)
- ... and their dependencies

## Note
`pyP` is inspired from `tp` developed by Amato Kasahara
