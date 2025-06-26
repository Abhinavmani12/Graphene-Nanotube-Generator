import math
import numpy as np
from collections import defaultdict

a = 0.142  # nm, C-C bond length
NM_TO_ANGSTROM = 10.0

def chiral_parameters(n, m):
    Ch_len = math.sqrt(3) * a * math.sqrt(n**2 + n*m + m**2)
    theta = math.atan((math.sqrt(3)*m) / (2*n + m))
    return Ch_len, theta

A1 = np.array([math.sqrt(3)*a, 0.0])
A2 = np.array([math.sqrt(3)*a/2, 3*a/2])

def generate_graphene_in_parallelogram(n, m, basis, T_len):
    Ch_vec = n * A1 + m * A2
    T_dir = np.array([-Ch_vec[1], Ch_vec[0]])
    T_dir /= np.linalg.norm(T_dir)
    T_vec = T_dir * T_len
    M = np.column_stack((A1, A2))
    M_inv = np.linalg.inv(M)
    corners = np.array([[0,0], Ch_vec, T_vec, Ch_vec + T_vec])
    corners_ij = np.dot(M_inv, corners.T).T
    i_min, j_min = corners_ij.min(axis=0).astype(int) - 1
    i_max, j_max = corners_ij.max(axis=0).astype(int) + 2
    atoms = []
    for i in range(i_min, i_max):
        for j in range(j_min, j_max):
            lattice_point = i * A1 + j * A2
            for b in basis:
                pos = lattice_point + b
                try:
                    coeffs = np.linalg.solve(np.column_stack((Ch_vec, T_vec)), pos)
                except np.linalg.LinAlgError:
                    continue
                alpha, beta = coeffs
                if 0 <= alpha <= 1 and 0 <= beta <= 1:
                    atoms.append(pos)
    return np.array(atoms), Ch_vec

def rectangle_to_cylinder(coords_2d, Ch_len, theta):
    radius = Ch_len / (2 * math.pi)
    cos_t = math.cos(-theta)
    sin_t = math.sin(-theta)
    wrapped = np.zeros((len(coords_2d), 3))
    for idx, (ox, oy) in enumerate(coords_2d):
        x_r = ox * cos_t - oy * sin_t
        y_r = ox * sin_t + oy * cos_t
        phi = 2 * math.pi * (x_r / Ch_len)
        wrapped[idx] = [radius * math.cos(phi), radius * math.sin(phi), y_r]
    return wrapped

def deduplicate(coords_3d, tol_nm=0.002):
    unique = []
    idx_map = {}
    for i, p in enumerate(coords_3d):
        for j, q in enumerate(unique):
            if np.linalg.norm(p - q) < tol_nm:
                idx_map[i] = j
                break
        else:
            idx_map[i] = len(unique)
            unique.append(p)
    return np.asarray(unique), idx_map

def find_bonds(coords, min_nm=0.130, max_nm=0.147):
    bonds = set()
    N = len(coords)
    for i in range(N):
        for j in range(i + 1, N):
            dist = np.linalg.norm(coords[i] - coords[j])
            if min_nm <= dist <= max_nm:
                bonds.add((i + 1, j + 1))
    return bonds

def validate_coordination(N_atoms, bonds):
    counts = [0]*N_atoms
    valid = []
    for i, j in sorted(bonds):
        if counts[i-1] < 3 and counts[j-1] < 3:
            valid.append((i, j))
            counts[i-1] += 1
            counts[j-1] += 1
    return valid

def remove_undercoordinated_atoms(coords, bonds, min_bonds=2):
    bond_map = defaultdict(set)
    for i, j in bonds:
        bond_map[i-1].add(j-1)
        bond_map[j-1].add(i-1)
    keep_indices = {i for i, nbrs in bond_map.items() if len(nbrs) >= min_bonds}
    new_idx_map = {old: new for new, old in enumerate(sorted(keep_indices))}
    coords_filtered = np.array([coords[i] for i in sorted(keep_indices)])
    bonds_filtered = []
    for i, j in bonds:
        i0, j0 = i-1, j-1
        if i0 in keep_indices and j0 in keep_indices:
            bonds_filtered.append((new_idx_map[i0]+1, new_idx_map[j0]+1))
    return coords_filtered, bonds_filtered

def add_seaming_bonds(coords, bonds, max_bond_dist=0.147, min_bond_dist=0.130):
    zs = coords[:, 2]
    z_min = zs.min()
    z_max = zs.max()
    edge_tol = 0.05
    start_edge_atoms = [i for i, z in enumerate(zs) if abs(z - z_min) < edge_tol]
    end_edge_atoms = [i for i, z in enumerate(zs) if abs(z - z_max) < edge_tol]
    existing_bonds = set(tuple(sorted((i-1, j-1))) for i, j in bonds)
    new_bonds = []
    for i in start_edge_atoms:
        for j in end_edge_atoms:
            dist = np.linalg.norm(coords[i] - coords[j])
            if min_bond_dist <= dist <= max_bond_dist:
                bond_tuple = tuple(sorted((i, j)))
                if bond_tuple not in existing_bonds:
                    new_bonds.append((i+1, j+1))
    print(f"Adding {len(new_bonds)} seaming bonds.")
    return bonds + new_bonds

def write_pdb(fname, coords, bonds, n, m):
    import datetime
    with open(fname, "w") as f:
        f.write("HEADER    MULTIWALL CARBON NANOTUBE\n")
        f.write(f"REMARK    (n,m)=({n},{m})\n")
        f.write(f"REMARK    {datetime.datetime.now():%Y-%m-%d %H:%M:%S}\n")
        f.write("MODEL        1\n")
        for i, (x,y,z) in enumerate(coords, 1):
            xA, yA, zA = (np.array([x,y,z]) * NM_TO_ANGSTROM)
            f.write(f"ATOM  {i:5d}  C   CNT A   1    "
                    f"{xA:8.3f}{yA:8.3f}{zA:8.3f}  1.00  0.00           C\n")
        cmap = {}
        for i, j in bonds:
            cmap.setdefault(i, []).append(j)
            cmap.setdefault(j, []).append(i)
        for i in sorted(cmap):
            parts = [f"CONECT{i:5d}"] + [f"{j:5d}" for j in cmap[i][:4]]
            f.write("".join(parts).ljust(66) + "\n")
        f.write("ENDMDL\nEND\n")
    print(f"PDB written: {fname}")

def mirror_coords(coords, plane='YZ'):
    mirrored = coords.copy()
    if plane == 'YZ':
        # Mirror across YZ plane = reflect X
        mirrored[:, 0] = -mirrored[:, 0]
    elif plane == 'XZ':
        # Mirror across XZ plane = reflect Y
        mirrored[:, 1] = -mirrored[:, 1]
    elif plane == 'XY':
        # Mirror across XY plane = reflect Z
        mirrored[:, 2] = -mirrored[:, 2]
    else:
        raise ValueError("Plane must be one of 'XY', 'YZ', or 'XZ'")
    return mirrored

# MAIN EXECUTION
print("Multiwalled Carbon Nanotube Generator")
n = int(input("Enter chiral index n (e.g. 10): "))
m = int(input("Enter chiral index m (e.g. 10): "))
orig_n, orig_m = n, m  # store original inputs
n_walls = int(input("Enter number of walls (e.g. 3): "))

if m == 0:
    T_len = 4 * math.sqrt(3) * a  # 4 unit cells ~ 0.984 nm
else:
    T_len = 5.0  # nm

basis = np.array([[0.0, 0.0], [math.sqrt(3)*a/2, a/2]])
all_coords = []
all_bonds = []
atom_offset = 0

for wall in range(n_walls):
    wall_n = n + wall
    wall_m = m + wall
    Ch_len, theta = chiral_parameters(wall_n, wall_m)
    flat, _ = generate_graphene_in_parallelogram(wall_n, wall_m, basis, T_len)
    flat[:, 0] -= (flat[:, 0].max() + flat[:, 0].min()) / 2
    flat[:, 1] -= (flat[:, 1].max() + flat[:, 1].min()) / 2
    coords_3d = rectangle_to_cylinder(flat, Ch_len, theta)
    coords_3d, idx_map = deduplicate(coords_3d)
    bonds = find_bonds(coords_3d)
    bonds = validate_coordination(len(coords_3d), bonds)
    bonds = [(i + atom_offset, j + atom_offset) for i, j in bonds]
    all_coords.extend(coords_3d)
    all_bonds.extend(bonds)
    atom_offset += len(coords_3d)
    print(f"Wall {wall+1}: atoms={len(coords_3d)}, bonds={len(bonds)}")

print(f"Total: {len(all_coords)} atoms, {len(all_bonds)} bonds before filtering")
all_coords_np = np.array(all_coords)
all_coords_np, all_bonds_filtered = remove_undercoordinated_atoms(all_coords_np, all_bonds)
print(f"After filtering: {len(all_coords_np)} atoms, {len(all_bonds_filtered)} bonds")

# Apply final mirroring if needed (only if m > n)
if m > n:
    all_coords_np = mirror_coords(all_coords_np, plane='YZ')

all_bonds_filtered = add_seaming_bonds(all_coords_np, all_bonds_filtered)
output_pdb_name = f"cnt_{orig_n}_{orig_m}.pdb"
write_pdb(output_pdb_name, all_coords_np, all_bonds_filtered, orig_n, orig_m)
