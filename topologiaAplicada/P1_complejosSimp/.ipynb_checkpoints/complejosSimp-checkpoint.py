"""
Rama 
"""


"""
Clase Complejos simpliciales que permita almacenar complejos simpliciales introduciendo sus símplices maximales y un flotante (filtración) y que contenga métodos que permitan:
1.Calcular la dimensión del complejo simplicial.
2.Calcular el conjunto de todas las caras del complejo simplicial.
3.Calcular el conjunto de todas las caras de una dimensión dada.
4.Calcular la estrella de un símplice.
5.Calcular el link de un símplice.
6.Calcular la característica de Euler.
7.Calcular el número de componentes de conexas.
8.Añadir símplices nuevos (con su flotante).
9.Almacenar la lista de todos los símplices ordenados según el valor del flotante y en caso de que dos tengan el mismo flotante ordenar por dimensión (las caras aparecen primero)
10.Calcular el complejo simplicial formado por todos los símplices cuyo flotante asociado sea menor o igual que un flotante dado.
"""

from itertools import combinations
import numpy as np
from scipy.spatial import distance

class Simplex:
    """Initialization of a simplex with vertices and indice"""
    def __init__(self, vertices, indice=0.0):
        # Vertices should be a tuple or list of hashable items
        self.vertices = tuple(sorted(vertices))  # Sort vertices to maintain consistency
        self.indice = indice

    """Dimension of the simplex"""
    def dimension(self):
        # Dimension of the simplex is the number of vertices - 1
        return len(self.vertices) - 1

    """Faces of the simplex (subsets of vertices)"""
    def faces(self):
        # All subsets of size len(vertices) - 1 form faces
        return [Simplex(vertices, self.indice) for i in range(len(self.vertices)) for vertices in combinations(self.vertices, i)]
    
    """Comparison operators for sorting"""
    def __lt__(self, other):
        # First compare by dimension (number of vertices)
        if self.indice != other.indice:
            return self.indice < other.indice
        if len(self.vertices) != len(other.vertices):
            return len(self.vertices) < len(other.vertices)
        # If dimensions are the same, compare lexicographically by vertices
        return self.vertices < other.vertices

    """Equality check"""
    def __eq__(self, other):
        return self.vertices == other.vertices

    """Hash function for use in sets and dictionaries"""
    def __hash__(self):
        return hash(self.vertices)

    """String representation"""
    def __repr__(self):
        return f"Simplex{self.vertices} {self.indice}"
    
class SimplicialComplex:
    """Initialization of a simplicial complex"""
    def __init__(self, simplices=None):
            # Initialize the complex with a set of simplices
            self.simplices = set()
            self.lista_ordenada = []
            if simplices:
                for simplex in simplices:
                    self.add_simplex(simplex)

    """Add a simplex to the complex"""
    def add_simplex(self, simplex):
        if isinstance(simplex, Simplex):
            if len(self.lista_ordenada) == 0:
                self.simplices.add(simplex)
                self.lista_ordenada = sorted(self.simplices)
                return 
            
            is_contained = False
            new_bigger_simplex = False
            end_of_searching = False
            
            remove_simplex = set()

            index=0
            #As there are only simplex of the same index, we can clean comparing dimensions
            while index<=len(self.lista_ordenada)-1 and not end_of_searching:
                if(self.lista_ordenada[index].indice == simplex.indice):
                    new_bigger_simplex = False
                    if simplex > self.lista_ordenada[index]:
                        big_simplex = simplex
                        small_simplex = self.lista_ordenada[index]
                        new_bigger_simplex = True
                    else:
                        big_simplex = self.lista_ordenada[index]
                        small_simplex = simplex
                    is_contained = all(vertice in big_simplex.vertices for vertice in small_simplex.vertices)
                    if is_contained and new_bigger_simplex:
                        remove_simplex.add(self.lista_ordenada[index])
                    elif is_contained:
                        end_of_searching = True
                index = index+1
                for not_useful_simplex in remove_simplex:
                    self.simplices.discard(not_useful_simplex)
                    if not_useful_simplex in self.lista_ordenada:
                        self.lista_ordenada.remove(not_useful_simplex)
                if not end_of_searching:
                    self.simplices.add(simplex)
                    self.lista_ordenada = sorted(self.simplices)
        else:
            raise TypeError("Expected a Simplex object")
        

    def get_lista_ord(self):
        return self.lista_ordenada

    def get_simplices_indice_menor_n(self, n):
        lista_simplices = []
        i = 0
        while i<len(self.lista_ordenada) and n>=self.lista_ordenada[i].indice:
            lista_simplices.append(self.lista_ordenada[i])
            i = i+1
        return lista_simplices

"""1.Calcular la dimensión del complejo simplicial."""
def dimension(self):
    """Return the dimension of the simplicial complex (largest simplex dimension)"""
    return max(s.dimension() for s in self.simplices)

"""2.Calcular el conjunto de todas las caras del complejo simplicial."""
"""3.Calcular el conjunto de todas las caras de dimensión k del complejo simplicial."""
    # If no dimension value is given, it will return all posible faces
def get_faces(self, dim=-1):
        all_faces = set()
        for simplex in self.simplices:
            if dim < 0 or (dim >= 0 and simplex.dimension() == dim): 
                    all_faces.add(simplex)
            for face in simplex.faces():
                if dim < 0 or (dim >= 0 and face.dimension() == dim): 
                    all_faces.add(face)
        return sorted(all_faces)

"""4.Calcular la estrella de un símplice."""
def get_estrella(self, simplice):
    if not isinstance(simplice, Simplex):
        return "This was not instance of Simplex"
        
    actual_dim = simplice.dimension()
    all_faces = sorted(self.get_faces(), reverse=True)
    number_faces = len(all_faces)
    
    if simplice not in all_faces:
        return "there is not such simplice in the complex"
    estrella = set()
    actual_vertices = simplice.vertices
    index = 0
    while index < number_faces and all_faces[index].dimension() >= actual_dim:
        is_contained = all(vertice in all_faces[index].vertices for vertice in actual_vertices)
        if is_contained:
            estrella.add(all_faces[index])
        index += 1
    return sorted(estrella)

"""5.Calcular el link de un símplice."""
def get_link(self, simplice):
    estrella = self.get_estrella(simplice)
    vertices = simplice.vertices
    estrella_cerrada = set()
    link = set()
    for cara in estrella:
        if cara not in estrella_cerrada:
            estrella_cerrada.update(cara.faces())
    for cara in estrella_cerrada:
        intersection = any(vertice in cara.vertices for vertice in vertices)
        if not intersection:
            link.add(cara)
    return sorted(link)

"""6.Calcular la característica de Euler."""
def get_euler(self):
    """Calculate the Euler characteristic of the simplicial complex"""
    faces = self.get_faces()
    caracteristic = 0
    for face in faces:
        if face.dimension() != -1:
            caracteristic += (1 - 2 * (face.dimension() % 2))
    return caracteristic

"""7.Calcular el número de componentes conexas."""
def get_numero_componentes_conexas(self):
    """Calculate the number of connected components in the simplicial complex"""
    numero_componentes_conexas = 1
    vertices = set()
    complejo_simplicial = set(self.simplices)
    aux = complejo_simplicial.pop().vertices
    for vertice in aux:
        vertices.add(vertice)

    simplice_actual = None
    intersection = False
    while complejo_simplicial:
        for simplice in complejo_simplicial:
            simplice_actual = simplice
            intersection = any(vertice in simplice.vertices for vertice in vertices)
            if intersection:
                break
            
        if not intersection:
            numero_componentes_conexas += 1

        vertices.update(simplice_actual.vertices)
        complejo_simplicial.remove(simplice_actual)
    return numero_componentes_conexas

"""8.Añadir simplices nuevos (con su flotante)."""
def add_simplex(self, simplex):
    if not isinstance(simplex, Simplex):
        raise TypeError("Expected a Simplex object")

    if len(self.lista_ordenada) == 0:
        self.simplices.add(simplex)
        self.lista_ordenada = sorted(self.simplices)
        return 

    remove_simplex = self._find_and_remove_contained_simplices(simplex)
    self._update_simplices_and_list(simplex, remove_simplex)

def _find_and_remove_contained_simplices(self, simplex):
    remove_simplex = set()
    index = 0
    end_of_searching = False

    while index <= len(self.lista_ordenada) - 1 and not end_of_searching:
        if self.lista_ordenada[index].indice == simplex.indice:
            new_bigger_simplex, big_simplex, small_simplex = self._compare_simplices(simplex, self.lista_ordenada[index])
            is_contained = all(vertice in big_simplex.vertices for vertice in small_simplex.vertices)
            if is_contained and new_bigger_simplex:
                remove_simplex.add(self.lista_ordenada[index])
            elif is_contained:
                end_of_searching = True
        index += 1

    return remove_simplex

def _compare_simplices(self, simplex1, simplex2):
    if simplex1 > simplex2:
        return True, simplex1, simplex2
    else:
        return False, simplex2, simplex1

def _update_simplices_and_list(self, simplex, remove_simplex):
    for not_useful_simplex in remove_simplex:
        self.simplices.discard(not_useful_simplex)
        if not_useful_simplex in self.lista_ordenada:
            self.lista_ordenada.remove(not_useful_simplex)

    if not any(simplex in remove_simplex for simplex in self.lista_ordenada):
        self.simplices.add(simplex)
        self.lista_ordenada = sorted(self.simplices)
    
"""9.Almacenar la lista de todos los símplices ordenados según el valor del flotante y en caso de que dos tengan el mismo flotante ordenar por dimensión (las caras aparecen primero)"""
def get_lista_ord(self):
    return self.lista_ordenada

"""10.Calcular el complejo simplicial formado por todos los símplices cuyo flotante asociado sea menor o igual que un flotante dado:"""
def get_simplices_indice_menor_n(self, n):
    """Return the simplicial complex formed by all simplices with filtration value less than or equal to n"""
    lista_simplices = []
    i = 0
    while i < len(self.lista_ordenada) and n >= self.lista_ordenada[i].indice:
        lista_simplices.append(self.lista_ordenada[i])
        i += 1
    return lista_simplices

def __repr__(self):
    return f"SimplicialComplex({list(self.simplices)})"

def thresholdvalues(list):
    return [distance.euclidean(a,b) for a, b in combinations(list, 2)]

def vietoris_rips_filtration(points, threshold):
    """
    Construye la filtración de Vietoris-Rips para un conjunto de puntos.
    
    Parámetros:
    - points: una lista de puntos (coordenadas 2D).
    - max_dimension: dimensión máxima de los simplíces a considerar.
    - threshold: umbral de distancia para agregar aristas.
    
    Retorna:
    - filtration: una lista de objetos de SimplicialComplex, uno para cada paso en la filtración.
    """
    num_points = len(points)
    filtration = []

    # Paso 1: Crear todos los 0-simplices (vértices)
    simplices = [Simplex([i], 0) for i in range(num_points)]
    complex = SimplicialComplex(set(simplices))
    filtration.append(complex)
    
    print(f"Paso 1: Añadidos {num_points} 0-simplices (vértices)")

    # Obtener las distancias entre cada par de puntos
    distancias = [(distance.euclidean(points[i], points[j]), (i, j)) for i, j in combinations(range(num_points), 2)]
    
    
    # Paso 2: Crear todos los simplíces de mayor dimensión
    new_simplices = []
    for dist, pair in distancias:
        if dist <= threshold:
            complex.add_simplex(Simplex([pair[0],pair[1]],dist)) #FINALIZAR AQUIIIIII AÑADIR LOS SIMPLICES Y YA
            draw_points(points, dist);
        else:
            break  # Si la distancia excede el umbral, dejamos de añadir nuevas simplices

        # Si hemos añadido nuevos simplíces, los incorporamos al complejo simplicial
    if new_simplices:
        new_complex = SimplicialComplex(complex.simplices.union(new_simplices))
        filtration.append(new_complex)
        complex = new_complex

    return filtration

def my_v_r(points, threshold):
    n_points = len(points)
    print(n_points)
        
    '''
num_points = len(points)
        filtration = []
        
        # Step 1: Create all 0-simplices (vertices)
        simplices = [Simplex([i], 0) for i in range(num_points)]
        complex = SimplicialComplex(simplices)
        filtration.append(complex)
        
        # Step 2: Create all higher-dimensional simplices
        for k in range(1, max_dimension + 1):
            new_simplices = SimplicialComplex()
            for comb in combinations(range(num_points), k + 1):
                if all(distance.euclidean(points[i], points[j]) <= threshold for i, j in combinations(comb, 2)):
                    new_simplices.add_simplex(Simplex(comb, k))
            
            if new_simplices:
                for simplice in new_simplices.simplices:
                    new_complex = SimplicialComplex(complex.add_simplex(simplice))
                filtration.append(new_complex)
                '''
    
"""Examples"""
points = np.array([[0, 0], [1, 0], [0, 1], [1, 1]])  # 4 points in a 2D plane
threshold = 2

my_v_r(points, threshold)

# Create some simplices
#s1 = Simplex([0, 1])  # A 1-simplex (edge)
s2 = Simplex([1, 2, 3], 0.0)  # A 2-simplex (triangle)
s3 = Simplex([1, 2, 3, 4, 5], 0.0)
# Initialize a simplicial complex
complex1 = SimplicialComplex([s2])
print(complex1)
complex1.add_simplex(Simplex([5], 1.0))
complex1.add_simplex(Simplex([6], 1.0))
complex1.add_simplex(Simplex([2], 1.0))
print(complex1)
complex1.add_simplex(Simplex([2,4], 1.0))
print(complex1)
# 

#complex1.add_simplex(Simplex([1,2,3,4]))
#print(complex1)


# Add another simplex
#complex1.add_simplex(Simplex([2, 3]))
points = np.array([[0, 0], [1, 0], [0, 1], [1, 1]])  # 4 points in a 2D plane
max_dimension = 2  # We want simplices up to dimension 2 (i.e., triangles)
threshold = 1.5  # Distance threshold

filtration = vietoris_rips_filtration(points, threshold)

# Display the filtration
for i, f in enumerate(filtration):
    print(f"Filtration step {i}: {f}")
