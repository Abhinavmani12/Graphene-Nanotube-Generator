import numpy as np
from math import sqrt

# ───────────── Constants ─────────────
a_cc = 1.42
a = sqrt(3) * a_cc
a1 = np.array([a, 0])
a2 = np.array([a / 2, 3 * a_cc / 2])

# ───────────── Chiral Length ─────────────
def chiral_length(n, m):
    return a * sqrt(n ** 2 + m ** 2 + n * m)

# ───────────── Flat Graphene ─────────────
def generate_chiral_graphene(n, m, sheet_radius):
    atoms = []
    N = int(2 * sheet_radius / a_cc) + 3
    for i in range(-N, N):
        for j in range(-N, N):
            base = i * a1 + j * a2
            for dx, dy in [(0, 0), (0, 1)]:
                pos = base + np.array([0, dy * a_cc])
                if np.linalg.norm(pos[:2]) <= sheet_radius:
                    atoms.append((pos[0], pos[1], 0.0))
    return atoms

# ───────────── Cone Deformation ─────────────
def deform_cone(atoms, height, cone_radius):
    new_atoms = []
    for x, y, _ in atoms:
        r = np.sqrt(x ** 2 + y ** 2)
        z = height * (1 - r / cone_radius) if r <= cone_radius else 0.0
        new_atoms.append((x, y, z))
    return new_atoms

# ───────────── Manual PDB Writer ─────────────
def save_as_pdb(atoms, filename):
    with open(filename, 'w') as f:
        for i, (x, y, z) in enumerate(atoms, 1):
            f.write(
                f"ATOM  {i:5d}  C   GRP A{1:4d}    "
                f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00           C\n"
            )
        f.write("END\n")

# ───────────── Main ─────────────
def generate_cone_graphene(n, m):
    chiral_len = chiral_length(n, m)
    cone_radius = max(10, 0.5 * chiral_len)
    height = 0.6 * cone_radius
    sheet_radius = 2 * cone_radius

    print(f"Using n={n}, m={m}")
    print(f"→ Cone base radius: {cone_radius:.2f} Å")
    print(f"→ Cone height:      {height:.2f} Å")
    print(f"→ Sheet radius:     {sheet_radius:.2f} Å")

    atoms = generate_chiral_graphene(n, m, sheet_radius)
    deformed_atoms = deform_cone(atoms, height, cone_radius)
    save_as_pdb(deformed_atoms, "graphene_cone_fixed.pdb")
    print("✔ Saved as: graphene_cone_fixed.pdb")

# ───────────── Entry ─────────────
if __name__ == "__main__":
    n = int(input("Enter chiral index n: "))
    m = int(input("Enter chiral index m: "))
    generate_cone_graphene(n, m)
