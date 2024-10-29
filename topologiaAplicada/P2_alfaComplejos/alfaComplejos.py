"""
-1.Definir una función que calcule la filtración de alfa complejos asociada a un conjunto de puntos del plano.
-2.Definir una función que represente gráficamente dicho alfa complejo.
-3.Definir una función que calcule la filtración de complejos de Vietoris-Rips asociada a un conjunto de puntos.
"""
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
from scipy.spatial import distance, Delaunay, Voronoi, voronoi_plot_2d
import numpy as np
from itertools import combinations

""" Draw points, edges between them if their distance is below the threshold, 
    and circles of radius threshold/2 around each point."""
def draw_points(points, threshold=None):
    plt.figure(figsize=(6, 6))
    ax = plt.gca()
    
    # Plot points and circles
    for i, point in enumerate(points):
        plt.scatter(point[0], point[1], color='blue', zorder=2)
        plt.text(point[0] + 0.05, point[1], f'{i}', fontsize=12, zorder=3)  # Label the points
        
        # Draw a circle of radius threshold/2 around each point
        if threshold is not None:
            circle = Circle((point[0], point[1]), threshold/2, color='r', fill=False, linestyle='--', alpha=0.5, zorder=1)
            ax.add_patch(circle)
    
    # Optionally plot edges if a threshold is given
    if threshold is not None:
        for i, j in combinations(range(len(points)), 2):
            dist = distance.euclidean(points[i], points[j])
            if dist <= threshold:
                plt.plot([points[i][0], points[j][0]], [points[i][1], points[j][1]], 'k-', alpha=0.5, zorder=1)
    
    plt.title(f"Vietoris-Rips Complex (Threshold: {threshold})")
    plt.gca().set_aspect('equal', adjustable='box')
    plt.grid(True)
    plt.xlim(min(points[:,0]) - threshold, max(points[:,0]) + threshold)
    plt.ylim(min(points[:,1]) - threshold, max(points[:,1]) + threshold)
    plt.show()

"""-1.Definir una función que calcule la filtración de alfa complejos asociada a un conjunto de puntos del plano."""
def alpha_complex_filtration(points):
    delaunay = Delaunay(points)
    filtration = []
    
    for simplex in delaunay.simplices:
        vertices = points[simplex]
        circum_radius = np.max([distance.euclidean(vertices[i], vertices[j]) for i, j in combinations(range(len(vertices)), 2)]) / 2
        filtration.append(Simplex(simplex, circum_radius))
    
    return filtration

"""-2.Definir una función que represente gráficamente dicho alfa complejo."""
def draw_alpha_complex(points, filtration):
    plt.figure(figsize=(6, 6))
    ax = plt.gca()
    
    # Plot points
    for i, point in enumerate(points):
        plt.scatter(point[0], point[1], color='blue', zorder=2)
        plt.text(point[0] + 0.05, point[1], f'{i}', fontsize=12, zorder=3)  # Label the points
    
    # Plot edges and triangles
    for simplex in filtration:
        vertices = points[simplex.vertices]
        if len(vertices) == 2:
            plt.plot(vertices[:, 0], vertices[:, 1], 'k-', alpha=0.5, zorder=1)
        elif len(vertices) == 3:
            triangle = plt.Polygon(vertices, edgecolor='k', fill=None, alpha=0.5, zorder=1)
            ax.add_patch(triangle)
    
    plt.title("Alpha Complex")
    plt.gca().set_aspect('equal', adjustable='box')
    plt.grid(True)
    plt.show()

"""-3.Definir una función que calcule la filtración de complejos de Vietoris-Rips asociada a un conjunto de puntos."""
def vietoris_rips_filtration(points, threshold):
    num_points = len(points)
    filtration = []

    # Step 1: Create all 0-simplices (vertices)
    simplices = [Simplex([i], 0) for i in range(num_points)]
    complex = SimplicialComplex(set(simplices))
    filtration.append(complex)
    
    print(f"Step 1: Added {num_points} 0-simplices (vertices)")

    # Get distances between each pair of points
    distancias = [(distance.euclidean(points[i], points[j]), (i, j)) for i, j in combinations(range(num_points), 2)]
    
    # Step 2: Create all higher-dimensional simplices
    for dist, pair in distancias:
        if dist <= threshold:
            complex.add_simplex(Simplex([pair[0], pair[1]], dist))
        else:
            break  # Stop adding new simplices if the distance exceeds the threshold

    return filtration

"""Examples"""
points = np.array([[0, 0], [1, 1], [2, 3], [3, 6]])  # 4 points in a 2D plane
threshold = 2.5

# Draw Vietoris-Rips complex
draw_points(points, threshold)

# Compute and draw Alpha complex
alpha_filtration = alpha_complex_filtration(points)
draw_alpha_complex(points, alpha_filtration)

# Compute Vietoris-Rips filtration
vr_filtration = vietoris_rips_filtration(points, threshold)

##checkeando rama quuique