# 🧪 Graphene–Nanotube Generator

A Python-based molecular modeling toolkit to generate **carbon nanotubes (CNTs)** and **graphene-cone junctions** from 2D chiral graphene sheets using vector geometry and atomic transformation.

## 🚀 Features

- Generate **single- or multi-walled carbon nanotubes** based on (n, m) chirality.
- Support for **cone-shaped nanotube ends** attached to flat graphene.
- Outputs in standard **PDB** and **CIF** formats, compatible with:
  - Jmol
  - VMD
  - PyMOL
  - LAMMPS (with conversion)

## 🧬 Sample Output

| Structure              | Visualization |
|------------------------|---------------|
| Multiwall CNT (8, 8)   | ✅ Clean       |
| Graphene Cone (12, 4)  | ✅ Bonded      |
| CNT from flat sheet    | ✅ Seamless    |

## 📦 Installation

```bash
git clone https://github.com/Abhinavmani12/Graphene-Nanotube-Generator.git
cd Graphene-Nanotube-Generator
pip install -r requirements.txt
