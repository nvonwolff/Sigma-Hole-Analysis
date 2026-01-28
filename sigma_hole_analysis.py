import numpy as np
from skimage.measure import marching_cubes
from scipy.interpolate import RegularGridInterpolator
from scipy.spatial import cKDTree

# ---------------- Cube file reader ----------------
def read_cube_file(filename):
    with open(filename, "r") as f:
        lines = f.readlines()

    # First two lines are comments
    header = lines[2].split()
    natoms = int(header[0])
    origin = np.array(header[1:4], dtype=float)

    # Grid vectors (number of points + vector components)
    nx, *vx = lines[3].split()
    ny, *vy = lines[4].split()
    nz, *vz = lines[5].split()
    nx, ny, nz = int(nx), int(ny), int(nz)
    vx, vy, vz = np.array(vx, dtype=float), np.array(vy, dtype=float), np.array(vz, dtype=float)

    # Skip atom lines
    data_start = 6 + natoms

    # Read volumetric data
    values = []
    for line in lines[data_start:]:
        values.extend(map(float, line.split()))
    data = np.array(values)

    if data.size != nx * ny * nz:
        raise ValueError(f"Cube size mismatch: expected {nx*ny*nz}, got {data.size}")

    # ORCA: Z fastest, then Y, then X
    data = data.reshape((nx, ny, nz), order="C")
    grid_vectors = np.vstack([vx, vy, vz])
    return data, origin, grid_vectors

# ---------------- PQR writer ----------------
def write_pqr_file(filename, vertices, esp_values):
    with open(filename, 'w') as f:
        for i, (x, y, z) in enumerate(vertices):
            esp = esp_values[i]
            f.write(
                f"ATOM  {i:5d}  H   UNK X   1  "
                f"{x:8.3f}{y:8.3f}{z:8.3f}"
                f"{esp:8.4f}{1.5:8.4f}\n"
            )

# ---------------- Cube → PQR conversion ----------------
def generate_pqr_from_cube(eldens_file, esp_file, isovalue=0.001, output_pqr="sigma_hole.pqr"):
    eldens_data, origin, grid_vectors = read_cube_file(eldens_file)
    esp_data, esp_origin, esp_vectors = read_cube_file(esp_file)

    if not np.allclose(origin, esp_origin) or not np.allclose(grid_vectors, esp_vectors):
        raise ValueError("ELDENS and ESP cube grids do not match")

    # Bohr → Å
    BOHR_TO_ANG = 0.52917721092
    origin *= BOHR_TO_ANG
    grid_vectors *= BOHR_TO_ANG

    # Marching cubes (Z,Y,X)
    eldens_mc = np.transpose(eldens_data, (2,1,0))
    verts_mc, faces, _, _ = marching_cubes(eldens_mc, level=isovalue)

    # Convert voxel indices → real space
    verts_voxel = verts_mc[:, [2,1,0]]
    verts_real = verts_voxel @ grid_vectors + origin

    # ESP interpolation
    nx, ny, nz = esp_data.shape
    interpolator = RegularGridInterpolator(
        (np.arange(nx), np.arange(ny), np.arange(nz)),
        esp_data,
        bounds_error=False,
        fill_value=0.0
    )
    voxel_coords = np.linalg.solve(grid_vectors.T, (verts_real - origin).T).T
    esp_vals = interpolator(voxel_coords)

    write_pqr_file(output_pqr, verts_real, esp_vals)
    print(f"PQR written: {output_pqr}")

# ---------------- PQR reader ----------------
def read_pqr(filename):
    coords, esp = [], []
    with open(filename) as f:
        for line in f:
            if line.startswith("ATOM"):
                parts = line.split()
                x, y, z = float(parts[6]), float(parts[7]), float(parts[8])
                q = float(parts[9])  # ESP column
                coords.append([x, y, z])
                esp.append(q)
    return np.array(coords), np.array(esp)

# ---------------- Local extrema ----------------
def find_local_extrema(coords, values, radius=0.2, mode="max", tol=1e-4):
    """Find local maxima or minima in a neighborhood."""
    tree = cKDTree(coords)
    extrema = []
    for i, r in enumerate(coords):
        idx = tree.query_ball_point(r, radius)
        if len(idx) < 5:  # skip isolated points
            continue
        neighborhood = values[idx]
        if mode == "max" and values[i] >= neighborhood.max() - tol:
            extrema.append(i)
        elif mode == "min" and values[i] <= neighborhood.min() + tol:
            extrema.append(i)
    return np.array(extrema)

# ---------------- Sigma-hole detection ----------------
def sigma_hole_candidates(coords, esp, center, direction, max_angle_deg=30, min_distance=0.5):
    rvecs = coords - center
    distances = np.linalg.norm(rvecs, axis=1)
    mask = distances > min_distance
    rvecs, esp_sel = rvecs[mask], esp[mask]
    rhat = rvecs / np.linalg.norm(rvecs, axis=1)[:,None]
    cosang = np.dot(rhat, direction)
    angle_mask = cosang > np.cos(np.deg2rad(max_angle_deg))
    return coords[mask][angle_mask], esp_sel[angle_mask]

def detect_sigma_holes(coords, esp, halogen_positions, bonded_atom_positions,
                       max_angle_deg=30, min_distance=1.5):
    sigma_holes = []
    for X, A in zip(halogen_positions, bonded_atom_positions):
        bond_vec = X - A
        bond_vec /= np.linalg.norm(bond_vec)
        sigma_coords, sigma_esp = sigma_hole_candidates(
            coords, esp, center=X, direction=bond_vec,
            max_angle_deg=max_angle_deg, min_distance=min_distance
        )
        if len(sigma_esp) == 0: continue
        max_idx = sigma_esp.argmax()
        sigma_holes.append({
            "halogen_pos": X,
            "bond_vec": bond_vec,
            "sigma_hole_value": sigma_esp[max_idx],
            "sigma_hole_pos": sigma_coords[max_idx]
        })
    return sigma_holes

# ---------------- Write extrema + sigma holes ----------------
def write_extrema_pqr(filename, coords, esp, local_max_idx, local_min_idx, sigma_holes):
    with open(filename, "w") as f:
        atom_id = 0
        # maxima = red
        for i in local_max_idx:
            x, y, z = coords[i]; q = esp[i]
            f.write(f"ATOM  {atom_id:5d}  H   MAX X   1  {x:8.3f}{y:8.3f}{z:8.3f}{q:8.4f}{1.5:8.4f}\n")
            atom_id += 1
        # minima = blue
        for i in local_min_idx:
            x, y, z = coords[i]; q = esp[i]
            f.write(f"ATOM  {atom_id:5d}  H   MIN X   1  {x:8.3f}{y:8.3f}{z:8.3f}{q:8.4f}{1.5:8.4f}\n")
            atom_id += 1
        # sigma holes = green
        for sh in sigma_holes:
            x, y, z = sh["sigma_hole_pos"]; q = sh["sigma_hole_value"]
            f.write(f"ATOM  {atom_id:5d}  H   SIG X   1  {x:8.3f}{y:8.3f}{z:8.3f}{q:8.4f}{1.5:8.4f}\n")
            atom_id += 1
    print(f"Extrema PQR written: {filename}")

if __name__ == "__main__":
    # 0. Ask user for input files and parameters
    eldens_file = input("Enter ELDENS cube file name (e.g., 4Fe.eldens.cube): ").strip()
    esp_file = input("Enter ESP cube file name (e.g., 4Fe.scfp.esp.cube): ").strip()
    
    try:
        isovalue = float(input("Enter isosurface value for electron density (e.g., 0.001): ").strip())
    except ValueError:
        print("Invalid isovalue, using default 0.001")
        isovalue = 0.001

    try:
        extrema_radius = float(input("Enter radius (Å) for local extrema detection (e.g., 0.2): ").strip())
    except ValueError:
        print("Invalid radius, using default 0.2 Å")
        extrema_radius = 0.2

    # 1. Convert cube → PQR
    generate_pqr_from_cube(
        eldens_file,
        esp_file,
        isovalue=isovalue,
        output_pqr="sigma_hole.pqr"
    )

    # 2. Read PQR
    coords, esp = read_pqr("sigma_hole.pqr")

    # 3. Find local extrema
    local_max_idx = find_local_extrema(coords, esp, radius=extrema_radius, mode="max", tol=1e-4)
    local_min_idx = find_local_extrema(coords, esp, radius=extrema_radius, mode="min", tol=1e-4)

    # 4. Define halogen positions / bonded atoms (example, replace with your real coords)
    halogen_positions = np.array([[1.188, -1.930, -6.651]])
    bonded_atom_positions = np.array([[-0.553, -2.266, -7.398]])

    # 5. Detect sigma holes
    sigma_holes = detect_sigma_holes(
        coords, esp,
        halogen_positions,
        bonded_atom_positions,
        max_angle_deg=30,
        min_distance=0.5
    )

    # 6. Print sigma hole strengths
    if len(sigma_holes) == 0:
        print("No sigma holes detected.")
    else:
        print("Detected sigma holes:")
        for i, sh in enumerate(sigma_holes, 1):
            print(f"  σ-hole {i}: ESP = {sh['sigma_hole_value']:.4f} a.u. at position {sh['sigma_hole_pos']}")

    # 7. Write extrema + sigma holes to PQR
    write_extrema_pqr(
        "extrema.pqr",
        coords,
        esp,
        local_max_idx,
        local_min_idx,
        sigma_holes
    )
