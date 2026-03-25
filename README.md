# Bouncing Droplets

Direct numerical simulation code infrastructure for **drop impact onto liquid pools in the inertio-capillary regime**, supporting collaborative work with the Harris Lab at Brown and with colleagues at the University of Warwick. The code complements the publication available at [https://doi.org/10.1017/jfm.2023.88](https://doi.org/10.1017/jfm.2023.88) and the [mathematical modelling repository](https://github.com/harrislab-brown/BouncingDroplets), as well as current work in both bouncing and coalescence regimes.

---

## 📌 Features

✅ Axisymmetric Navier-Stokes solver for drop-pool impact scenarios  
✅ Parameter sweep support: drop velocity, resolution level  
✅ Two-phase, non-coalescing VOF implementation for bouncing regime  
✅ Standard (coalescing) VOF implementation also provided  
✅ High-resolution output and animation capabilities  
✅ Full visualization and post-processing ready data  

---

## 🛠️ Installation

### 1. Requirements

- [Basilisk](http://basilisk.fr/) (C-based DNS solver)
  - See [installation page](http://basilisk.fr/src/INSTALL) for instructions
- C compiler
- Visualization tools (optional, can be disabled based on architecture)

### 2. Clone the repository

```bash
git clone https://github.com/rcsc-group/BouncingDroplets
cd BouncingDroplets
```

### 3. Install Basilisk

Follow [Basilisk's installation guide](http://basilisk.fr/src/INSTALL) for your system.

### 4. Core Dependencies

The code includes the **two-phase non-coalescing fluid volume implementation** by V. Sanjay (from [this repository](https://github.com/VatsalSy/Lifting-a-sessile-drop/blob/master/CaseI/two-phaseDOD.h)), which has been successfully employed to limit numerical artifacts during contact time.

---

## ⚙️ Running the Code

After Basilisk is set up, run the driver code using the provided shell script:

```bash
sh run_master_example.sh
```

### Key Parameters

Edit the shell script to control:

- Drop velocity $V_0$ 
- Resolution level
- Domain size and other computational parameters
- Visualization options (inside the driver code)

The script organizes outputs into folders with summaries, VOF data, interface coordinates, simulation slices, and animations for further post-processing.

## 📁 Folder Structure

```bash
.
├── DriverCode/
│   ├── MasterImpact/                           # Two-phase non-coalescing VOF implementation
│   │   ├── DropImpact.c                        # Main simulation source code (bouncing regime)
│   │   └── two-phaseDOD.h                      # Header for two-phase, non-coalescing VOF
│   │
│   ├── MasterImpactSingleVOF/                  # Standard single VOF implementation
│   │   └── DropImpact.c                        # Simulation source code (coalescing regime)
│   │
│   ├── run_master_example.sh                   # Shell script for parameter sweeps (bouncing regime)
│   └── run_SingleVOF.sh                        # Shell script for single VOF simulations (coalescence regime)
│
├── LICENSE                                     # License information
└── README.md                                   # Project documentation
```

---

## 📊 Outputs

The simulation generates:

- **Summary files**: DNS execution information and mass conservation metrics
- **VOF data**: Volume of Fluid field data
- **Interface data**: Droplet and pool interface coordinates
- **Visualizations**: Simulation slices and `.mp4` animations
- **Post-processing data**: Ready for further analysis

Visualization capabilities can be toggled depending on your local architecture and needs.

---

## 📚 Citation

If you use this code or data in your work, please cite the associated publication:

> Alventosa, L. F., Cimpeanu, R., & Harris, D. M. (2023). Inertio-capillary rebound of a droplet impacting a fluid bath. Journal of Fluid Mechanics, 958, A24.

BibTeX:
```bibtex
@article{alventosa2023,
  title={Inertio-capillary rebound of a droplet impacting a fluid bath},
  author={Alventosa, Luke FL and Cimpeanu, Radu and Harris, Daniel M},
  journal={Journal of Fluid Mechanics},
  volume={958},
  pages={A24},
  year={2023},
  publisher={Cambridge University Press}
}
```

---

## 🧑 Contributing

Feel free to:

- Fork this repo
- Open issues for bug reports or feature requests
- Submit pull requests with improvements

---
