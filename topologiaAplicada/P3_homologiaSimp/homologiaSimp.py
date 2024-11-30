"""
-Definir una función que calcule la forma normal de Smith de una matriz con coeficientes en Z2
-Definir dentro de la clase Complejo simplicial métodos que permitan:
    1.Calcular la matriz borde para cada dimensión.
    2.Calcular los números de Betti.
-Calcular los números de Betti de los siguientes complejos simpliciales:
    1.El tetraedro.
    2.El borde del tetraedro.
    3.El toro con las dos triangulaciones vistas en clase.
    4.El plano proyectivo.
    5.La botella de Klein.
    6.El anillo.
    7.El sombrero del asno.
    8.Del complejo simplicial de la transparencia 4 del documento Homología Simplicial II.
    9.Del doble toro.
    10.De algunos alfa complejos.
-Crear una función que calcule los números de Betti b_0 y b_1 de un complejo simplicial contenido en el plano utilizando el algoritmo incremental. 
-Calcular los números de Betti de algunos alfa complejos del plano utilizando el algoritmo incremental.
"""

# Función que calcula la forma normal de Smith de una matriz con coeficientes en Z2 (fuera del objeto Complejo simplicial)
def smith_normal_form_z2(matrix):
    """
    Calcula la forma normal de Smith de una matriz con coeficientes en Z2.
    Args:
        matrix (list of list of int): Matriz con valores en {0, 1}.
    Returns:
        np.ndarray: Matriz en forma normal de Smith.
    """
    mat = np.array(matrix, dtype=int) % 2
    rows, cols = mat.shape
    row, col = 0, 0

    while row < rows and col < cols:
        # Buscar pivote en la columna actual
        pivot_row = next((r for r in range(row, rows) if mat[r, col] == 1), None)

        if pivot_row is not None:
            # Intercambiar filas para mover el pivote a la posición (row, col)
            mat[[row, pivot_row]] = mat[[pivot_row, row]]
        else:
            # Si no hay pivote en la columna, buscar en la fila
            pivot_col = next((c for c in range(col, cols) if mat[row, c] == 1), None)
            if pivot_col is not None:
                # Intercambiar columnas para mover el pivote a la posición (row, col)
                mat[:, [col, pivot_col]] = mat[:, [pivot_col, col]]
            else:
                # Avanzar si no se encuentra pivote en columna ni en fila
                col += 1
                continue

        # Reducir otras filas usando el pivote
        for r in range(rows):
            if r != row and mat[r, col] == 1:
                mat[r] = (mat[r] + mat[row]) % 2

        # Reducir otras columnas usando el pivote
        for c in range(cols):
            if c != col and mat[row, c] == 1:
                mat[:, c] = (mat[:, c] + mat[:, col]) % 2

        row += 1
        col += 1

    return mat

'''
Ejemplo de smith: 
matrix = [
    [1, 0, 1, 1, 0],
    [0, 1, 1, 1, 1],
    [1, 1, 0, 1, 0],
    [1, 0, 1, 0, 1],
    [0, 1, 1, 0, 0],
    [1, 1, 0, 1, 1]
]
smith_form = smith_normal_form_z2(matrix)
print("Forma normal de Smith:")
print(smith_form)
'''

# Meter en el complejo simplicial estos métodos
    def calcular_matriz_borde(self, dim):
        """
        Calcula la matriz borde para los símplices de dimensión dim en un complejo simplicial.
        
        Args:
            complejo (SimplicialComplex): El complejo simplicial.
            dim (int): La dimensión para la cual calcular la matriz borde.
            
        Returns:
            np.ndarray: La matriz borde.
        """

        if(dim == 0):
            num_vertices = len(self.get_faces(0))
            return np.zeros((1, num_vertices), dtype=int)
        if(dim == self.dimension):
            num_hiperplanos = len(self.get_faces(self.dimension-1))
            return np.ones((num_vertices, 1), dtype=int)
        
        # Obtener los símplices de dimensión dim y dim-1
        simplices_dim = [s for s in self.get_faces(dim) if s.dimension() == dim]
        simplices_dim_menos_1 = [s for s in self.get_faces(dim - 1) if s.dimension() == dim - 1]
    
        # Crear un diccionario para mapear cada simplex a un índice
        indice_simplices_dim = {s: i for i, s in enumerate(simplices_dim)}
        indice_simplices_dim_menos_1 = {s: i for i, s in enumerate(simplices_dim_menos_1)}
    
        # Inicializar la matriz borde con ceros
        matriz_borde = np.zeros((len(simplices_dim_menos_1), len(simplices_dim)), dtype=int)
    
        # Llenar la matriz borde
        for j, simplex in enumerate(simplices_dim):
            for face in simplex.faces():
                if face in indice_simplices_dim_menos_1:
                    i = indice_simplices_dim_menos_1[face]
                    matriz_borde[i, j] = 1  # En Z2, sólo consideramos 0 o 1
    
        return matriz_borde

    def calcular_numero_betti_smith(self):
        """
        Calcula los números de Betti de un complejo simplicial usando la forma normal de Smith.
        
        Args:
            complejo (SimplicialComplex): El complejo simplicial.
        
        Returns:
            list: Los números de Betti para cada dimensión.
        """
        max_dim = self.dimension()  # Dimensión máxima del complejo
        #print(max_dim)
        betti_numbers = []
    
        for k in range(max_dim + 1):
            # Matriz borde para la dimensión k y k+1
            matriz_borde_k = self.calcular_matriz_borde(k)
            matriz_borde_kp1 = self.calcular_matriz_borde(k + 1) if k + 1 <= max_dim else None
    
            # Forma normal de Smith
            matriz_smith_k = smith_normal_form_z2(matriz_borde_k)
            matriz_smith_kp1 = smith_normal_form_z2(matriz_borde_kp1) if matriz_borde_kp1 is not None else None
            print("Matriz ",k)
            print(matriz_smith_k)
            print("Matriz ",k+1)
            print(matriz_smith_kp1)
            # Dimensiones
            dim_kernel_k = matriz_smith_k.shape[1] - np.sum(np.diag(matriz_smith_k) == 1) if matriz_smith_k.size > 0 else 0
            dim_image_kp1 = np.sum(np.diag(matriz_smith_kp1) == 1) if matriz_smith_kp1 is not None and matriz_smith_kp1.size > 0 else 0
            print("Dim Z",k,":")
            print(dim_kernel_k)
            print("Dim B",k,":")
            print(dim_image_kp1)
            # Número de Betti
            betti_k = dim_kernel_k - dim_image_kp1
            betti_numbers.append(betti_k)
    
        return betti_numbers
    

# Tetraedro
simplices = [
    Simplex([0, 1, 2, 3])
]
complejo = SimplicialComplex(simplices)
numeros_de_betti = complejo.calcular_numero_betti_smith()
print(numeros_de_betti)

# Borde del Tetraedro
simplices = [
    Simplex([0, 1, 2]), Simplex([0, 1, 3]), Simplex([0, 2, 3]), Simplex([1, 2, 3])
]
complejo = SimplicialComplex(simplices)
numeros_de_betti = complejo.calcular_numero_betti_smith()
print(numeros_de_betti)


# Anillo
simplices = [
    Simplex([0, 1, 3]), Simplex([0, 2, 5]), Simplex([0, 3, 5]), Simplex([1, 2, 4]), 
    Simplex([1, 3, 4]), Simplex([2, 4, 5])
]
complejo = SimplicialComplex(simplices)
print(complejo)
numeros_de_betti = complejo.calcular_numero_betti_smith()
print(numeros_de_betti)


# Toro 1
toro1 = [
    Simplex([0, 1, 2]), Simplex([0, 2, 3]), Simplex([0, 3, 4]),
    Simplex([0, 4, 5]), Simplex([0, 5, 1]), Simplex([1, 2, 5]),
    Simplex([2, 3, 5]), Simplex([3, 4, 5])
]
complejoToro1 = SimplicialComplex(toro1)
numeros_de_betti = complejoToro1.calcular_numero_betti_smith()
print(numeros_de_betti)


# Toro 2
toro2 = [
    Simplex([0, 1, 3]), Simplex([1, 2, 3]), Simplex([2, 3, 4]),
    Simplex([2, 4, 5]), Simplex([3, 4, 5]), Simplex([0, 3, 5]),
    Simplex([0, 1, 5]), Simplex([1, 2, 5])
]
complejoToro2 = SimplicialComplex(toro2)
numeros_de_betti = complejoToro2.calcular_numero_betti_smith()
print(numeros_de_betti)