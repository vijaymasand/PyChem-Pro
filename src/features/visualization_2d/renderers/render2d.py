import matplotlib.pyplot as plt
import math

class Renderer2D:
    def __init__(self, molecule, coords):
        self.molecule = molecule
        self.coords = coords

    def draw(self):
        fig, ax = plt.subplots(figsize=(10, 8))
        ax.set_aspect('equal')
        ax.axis('off')

        has_label = {}
        for atom in self.molecule.atoms:
            degree = len(list(self.molecule.get_neighbors(atom.index)))
            if atom.symbol == 'C' and degree > 1 and atom.formal_charge == 0:
                has_label[atom.index] = False
            else:
                has_label[atom.index] = True

        # 1. Draw Bonds
        for bond in self.molecule.bonds:
            idx1, idx2 = bond.begin_atom_idx, bond.end_atom_idx
            if idx1 not in self.coords or idx2 not in self.coords: continue

            x1, y1 = self.coords[idx1]
            x2, y2 = self.coords[idx2]

            dx, dy = x2 - x1, y2 - y1
            dist = math.hypot(dx, dy)
            if dist == 0: continue

            # FIX 1: Significantly reduced the shrink amount so bonds get closer to text
            shrink1 = 0.22 if has_label[idx1] else 0.0
            shrink2 = 0.22 if has_label[idx2] else 0.0

            if dist > (shrink1 + shrink2 + 0.05):
                nx, ny = dx/dist, dy/dist
                sx, sy = x1 + nx * shrink1, y1 + ny * shrink1
                ex, ey = x2 - nx * shrink2, y2 - ny * shrink2
            else:
                sx, sy, ex, ey = x1, y1, x2, y2

            if bond.is_double:
                self._draw_double_bond(ax, sx, sy, ex, ey)
            else:
                ax.plot([sx, ex], [sy, ey], color='black', lw=2.5, zorder=1)

        # 2. Draw Atoms
        for atom in self.molecule.atoms:
            if atom.index not in self.coords or not has_label[atom.index]:
                continue

            x, y = self.coords[atom.index]
            label = self._build_atom_label(atom)

            # FIX 2: Removed the "Smart Nudging" block.
            # Text now sits perfectly on the coordinate node.

            ax.text(x, y, label, color='black', fontsize=22,
                    ha='center', va='center', zorder=2, fontweight='bold',
                    clip_on=False,
                    bbox=dict(facecolor='white', edgecolor='none', pad=0.05))

        # Manual Axis Padding
        if self.coords:
            xs = [c[0] for c in self.coords.values()]
            ys = [c[1] for c in self.coords.values()]
            margin = 1.0 # Tighter margin
            ax.set_xlim(min(xs) - margin, max(xs) + margin)
            ax.set_ylim(min(ys) - margin, max(ys) + margin)

        plt.tight_layout()
        plt.show()

    def _build_atom_label(self, atom):
        label = atom.symbol
        h_count = atom.num_implicit_h

        if h_count == 1:
            label += "H"
        elif h_count > 1:
            label += f"H$_{{{h_count}}}$"

        if atom.formal_charge > 0:
            label += f"$^{{+{atom.formal_charge if atom.formal_charge > 1 else ''}}}$"
        elif atom.formal_charge < 0:
            label += f"$^{{-{abs(atom.formal_charge) if abs(atom.formal_charge) > 1 else ''}}}$"

        return label

    def _draw_double_bond(self, ax, x1, y1, x2, y2):
        dx, dy = x2 - x1, y2 - y1
        length = math.hypot(dx, dy)
        if length == 0: return
        nx, ny = -dy / length, dx / length

        # FIX 3: Tighter offset and smaller shrink for cohesive double bonds
        offset = 0.12

        ax.plot([x1, x2], [y1, y2], color='black', lw=2.5, zorder=1)

        shrink_ratio = 0.10 # Reduced from 0.15 so inner lines reach closer to the corners

        ix1, iy1 = x1 + dx * shrink_ratio + nx * offset, y1 + dy * shrink_ratio + ny * offset
        ix2, iy2 = x2 - dx * shrink_ratio + nx * offset, y2 - dy * shrink_ratio + ny * offset

        ax.plot([ix1, ix2], [iy1, iy2], color='black', lw=2.5, zorder=1)