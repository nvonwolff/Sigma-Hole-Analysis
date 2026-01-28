# Sigma Hole & Local Extrema Analysis from Cube Files

This repository provides a Python script to generate **PQR files from ORCA electron density and ESP cube files**, detect **sigma holes**, and locate **local extrema** of the electrostatic potential. The generated PQR files can be visualized in **VMD**.

---

## Features

* Convert electron density + ESP cube files to a **PQR surface file** (`sigma_hole.pqr`).
* Detect **local maxima, minima**, and **sigma holes**.
* Generate an **extrema PQR** (`extrema.pqr`) for easy visualization of critical points.
* Fully adjustable **isosurface value** and **neighborhood radius** for extrema search.

---

## Requirements

* **Python 3.10+** recommended
* Create a virtual environment and install dependencies:

```bash
# Create virtual environment
python -m venv .venv

# Activate environment
source .venv/bin/activate  # Linux / Mac
.venv\Scripts\activate     # Windows

# Install dependencies
pip install numpy scipy scikit-image
```

* Optional: **VMD** for visualization of generated PQR files: [VMD Download](https://www.ks.uiuc.edu/Research/vmd/)

---

## Usage

### 0.  Generate Orca cube files

  1.  Generate eldens.cube file from any orca run (tested on Orca 6.1.1)

  -  Make sure that the density is saved and generate the eldens.cube file directly with adjusted grid if desired, e.g.:

     ```
       ! RKS wB97M-V def2-TZVPP def2-tzvpp/c def2/j tightscf rijcosx keepdens ### OR YOUR METHOD OF CHOICE ###
       %plots
        dim1    200
        dim2    200
        dime3   200
        Format  Gaussian_Cube
        Eldens("filename.eldens.cube");
       end
     ```

  This should give you the `filename.eldens.cube` file
  
  2.  Generate esp.cube file via orca_plot

  -  Run orca_plot: 
  `orca_plot filename.gbw -i`

  -  type
  `1 - Enter ypye of plot` and chose `2 - (scf) electron density` or `43 - Electrostatic potential`

  -  define grid
  `4 - Enter number of grid intervals` and enter same grid as in intial orca.inp

  -  select output format `5 - Slect output file format` and chose `7 - 3D  Gaussian cube`

  -  generate the plot `11 - Generate the plot` this will generate the corresponding cube file.
  

### 1. Generate PQR from cube files

Run the script and provide inputs when prompted:

```bash
python sigma_hole_analysis.py
```

You will be asked for:

* Electron density cube file (ELDENS)
* Electrostatic potential cube file (ESP)
* Isosurface value for the PQR
* Neighborhood radius for extrema detection

This produces:

* `sigma_hole.pqr` → isosurface ESP values
* `extrema.pqr` → highlights local maxima, minima, and sigma holes

---

### 2. Workflow Diagram

```
  ┌───────────────────────────┐
  │  ORCA cube files          │
  │  (ELDENS & ESP)           │
  └─────────────┬─────────────┘
                │
                ▼
  ┌───────────────────────────┐
  │  sigma_hole_analysis.py   │
  │  → sigma_hole.pqr         │
  │  → extrema.pqr            │
  └─────────────┬─────────────┘
                │
                ▼
  ┌───────────────────────────┐
  │  VMD visualization        │
  │  - sigma_hole.pqr: ESP    │
  │    isosurface             │
  │  - extrema.pqr: maxima,   │
  │    minima, sigma holes    │
  └───────────────────────────┘
```

---

### 3. Visualize in VMD

#### Isosurface:

1. Open VMD → `File → New Molecule…`
2. Load `sigma_hole.pqr` → click `Load`
3. Adjust representation:

   * Drawing Method → `VDW` or `Licorice`
   * Coloring Method → `Charge` (ESP values)
   * Adjust radius scaling for visibility

#### Extrema and Sigma Holes:

1. Open VMD → `File → New Molecule…`
2. Load `extrema.pqr` → click `Load`
3. Points legend:

   * `SIG` → sigma holes (green)
   * `MAX` → local maxima (red)
   * `MIN` → local minima (blue)
4. Adjust representation → `VDW` for points, `Radius` controls size

#### Combined Visualization:

* Load both `sigma_hole.pqr` and `extrema.pqr` in the same session to visualize **sigma holes and extrema on the ESP isosurface**.

---

### 4. Adjusting Parameters

* **Isosurface value** → controls the density threshold of the surface
* **Radius** → determines neighborhood size for detecting local extrema

Increasing the radius will detect fewer, more robust extrema; decreasing the radius will find more points, including the "negative belt" around sigma holes.

---

### 5. Example Output

* `sigma_hole.pqr` → surface colored by ESP
* `extrema.pqr` → points representing local maxima, minima, and sigma holes

---

### 6. License

MIT License
