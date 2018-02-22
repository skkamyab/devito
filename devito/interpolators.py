from collections import OrderedDict
from math import floor

import sympy

from devito.symbolics import indexify, retrieve_indexed


class LinearInterpolator(object):
    def __init__(self, grid, coordinates):
        self.grid = grid
        self.coordinates = coordinates
        self.npoint = 2**grid.dim
        self.r = 1

    def calculate_coefficients(self, array, args):
        print("Calculating coefficients")
        array[:] = 1

    def calculate_grid_points(self, array, args):
        assert(array.shape == self.coordinates.shape)
        print(args)
        for i, point in enumerate(self.coordinates.data):
            for j in range(self.grid.dim):
                origin = args[self.grid.origin[j].name]
                spacing = self.grid.spacing[j]
                
                array[i, j] = (point[j] - origin)/spacing

    def coefficients(self):
        """Symbolic expression for the coefficients for sparse point
        interpolation according to:
        https://en.wikipedia.org/wiki/Bilinear_interpolation.

        :returns: List of coefficients, eg. [b_11, b_12, b_21, b_22]
        """
        # Grid indices corresponding to the corners of the cell
        x1, y1, z1, x2, y2, z2 = sympy.symbols('x1, y1, z1, x2, y2, z2')
        # Coordinate values of the sparse point
        px, py, pz = self.point_symbols
        if self.grid.dim == 2:
            A = sympy.Matrix([[1, x1, y1, x1*y1],
                              [1, x1, y2, x1*y2],
                              [1, x2, y1, x2*y1],
                              [1, x2, y2, x2*y2]])

            p = sympy.Matrix([[1],
                              [px],
                              [py],
                              [px*py]])

            # Map to reference cell
            x, y = self.grid.dimensions
            reference_cell = {x1: 0, y1: 0, x2: x.spacing, y2: y.spacing}

        elif self.grid.dim == 3:
            A = sympy.Matrix([[1, x1, y1, z1, x1*y1, x1*z1, y1*z1, x1*y1*z1],
                              [1, x1, y2, z1, x1*y2, x1*z1, y2*z1, x1*y2*z1],
                              [1, x2, y1, z1, x2*y1, x2*z1, y2*z1, x2*y1*z1],
                              [1, x1, y1, z2, x1*y1, x1*z2, y1*z2, x1*y1*z2],
                              [1, x2, y2, z1, x2*y2, x2*z1, y2*z1, x2*y2*z1],
                              [1, x1, y2, z2, x1*y2, x1*z2, y2*z2, x1*y2*z2],
                              [1, x2, y1, z2, x2*y1, x2*z2, y1*z2, x2*y1*z2],
                              [1, x2, y2, z2, x2*y2, x2*z2, y2*z2, x2*y2*z2]])

            p = sympy.Matrix([[1],
                              [px],
                              [py],
                              [pz],
                              [px*py],
                              [px*pz],
                              [py*pz],
                              [px*py*pz]])

            # Map to reference cell
            x, y, z = self.grid.dimensions
            reference_cell = {x1: 0, y1: 0, z1: 0, x2: x.spacing,
                              y2: y.spacing, z2: z.spacing}
        else:
            raise NotImplementedError('Interpolation coefficients not implemented '
                                      'for %d dimensions.' % self.grid.dim)

        A = A.subs(reference_cell)
        return A.inv().T.dot(p)
    
    def gridpoints(self, point, grid):
        # List of indirection indices for all adjacent grid points
        index_matrix = [tuple(idx + ii + offset for ii, idx
                              in zip(inc, self.coordinate_indices))
                        for inc in self.point_increments]

        # Generate index substituions for all grid variables
        idx_subs = []
        for i, idx in enumerate(index_matrix):
            v_subs = [(v, v.base[v.indices[:-self.grid.dim] + idx])
                      for v in variables]
            idx_subs += [OrderedDict(v_subs)]
        return idx_subs
    
