import math
import numpy as np

def generate_sphere_points(n_points=160):
    """Generate uniformly distributed points on a sphere using Golden Section spiral."""
    points = []
    offset = 2.0 / n_points
    increment = math.pi * (3.0 - math.sqrt(5.0))  # golden angle
    for i in range(n_points):
        y = ((i * offset) - 1) + (offset / 2)
        r = math.sqrt(1 - y * y)
        phi = i * increment
        x = math.cos(phi) * r
        z = math.sin(phi) * r
        points.append([x, y, z])
    return np.array(points, dtype=np.float32)

def compute_sasa(molecule, probe_radius=1.4, n_sphere_points=5000):
    """
    Compute Shrake-Rupley Solvent Accessible Surface Area (SASA) for each atom.
    Assigns the result to `atom.sasa` and `atom.sasa_points` (for visualization).
    
    Performance improvements inspired by pymol-rs's SpatialHash algorithm:
    Uses scipy.spatial.cKDTree for spatial neighbor acceleration and fast NumPy
    vectorization for point-cloud collision testing, turning O(N^2) checks into O(N log N).
    """
    if not molecule.atoms:
        return 0.0

    points = generate_sphere_points(n_sphere_points)
    point_area = 4 * math.pi / n_sphere_points
    
    # Pre-extract data for speed
    positions = []
    radii = []
    for atom in molecule.atoms:
        if atom.has_coords:
            positions.append([atom.x, atom.y, atom.z])
            # vdw_radius + probe_radius
            radii.append(atom.element.vdw_radius + probe_radius)
        else:
            positions.append([0.0, 0.0, 0.0])
            radii.append(0.0)
            
    positions = np.array(positions, dtype=np.float32)
    radii = np.array(radii, dtype=np.float32)
    radii_sq = radii ** 2
    
    # Build spatial hash/tree for accelerated lookups (like SpatialHash in pymol-rs)
    from scipy.spatial import cKDTree
    tree = cKDTree(positions)
    max_radius = np.max(radii) if len(radii) > 0 else 0.0
    
    total_sasa = 0.0
    
    for i, atom in enumerate(molecule.atoms):
        atom.sasa = 0.0
        atom.sasa_points = []
        if not atom.has_coords:
            continue
            
        center_i = positions[i]
        radius_i = radii[i]
        
        # Find neighboring atoms that could potentially intersect with atom i's sphere
        # The maximum possible distance for overlap is radius_i + max_radius
        neighbor_indices = tree.query_ball_point(center_i, radius_i + max_radius)
        
        # Filter out self
        neighbor_indices = [idx for idx in neighbor_indices if idx != i]
        
        # Test points forming a sphere around atom i
        test_points = center_i + radius_i * points
        
        if not neighbor_indices:
            # No neighbors, fully accessible
            accessible_count = n_sphere_points
            accessible_dots = test_points.tolist()
        else:
            neighbor_positions = positions[neighbor_indices]
            neighbor_radii_sq = radii_sq[neighbor_indices]
            
            # Fully vectorized point collision check
            # test_points: (5000, 3), neighbor_positions: (M, 3) -> dx: (5000, M, 3)
            dx = test_points[:, np.newaxis, :] - neighbor_positions[np.newaxis, :, :]
            dist_sq = np.sum(dx * dx, axis=2)
            
            # Check if each point intersects with ANY neighbor
            # dist_sq: (5000, M), neighbor_radii_sq: (M,)
            is_intersecting = np.any(dist_sq < neighbor_radii_sq, axis=1)
            is_accessible = ~is_intersecting
            
            accessible_count = np.sum(is_accessible)
            # Convert accessible points directly to list for assignment
            accessible_dots = test_points[is_accessible].tolist()
                
        sasa_i = accessible_count * point_area * (radius_i ** 2)
        atom.sasa = sasa_i
        # Store as tuples for compatibility if needed, or lists. tolist() gives lists.
        # The original code appended tuples, so we map to tuples for exact compatibility.
        atom.sasa_points = [tuple(p) for p in accessible_dots]
        total_sasa += sasa_i
        
    molecule.properties['sasa'] = total_sasa
    return total_sasa

def compute_center_of_mass(molecule):
    """Compute the mass-weighted Center of Mass (COM) for the molecule."""
    total_mass = 0.0
    com = np.zeros(3, dtype=np.float64)
    
    for atom in molecule.atoms:
        if atom.has_coords:
            mass = atom.element.mass
            com[0] += mass * atom.x
            com[1] += mass * atom.y
            com[2] += mass * (atom.z or 0.0)
            total_mass += mass
            
    if total_mass > 0:
        com /= total_mass
        
    molecule.properties['center_of_mass'] = tuple(com)
    return tuple(com)
