import numpy as np
from scipy.linalg import expm


class HHLAlgorithm:
    def __init__(self, A, b, t, M):
        """
        Initializes the HHLAlgorithm object with the input matrix A,
        the right-hand side vector b, the time evolution parameter t,
        and the number of qubits M used to represent the solution.
        """
        self.A = A
        self.b = b
        self.t = t
        self.M = M

    def solve(self):
        """
        Solves the linear system of equations Ax = b using the HHL algorithm,
        and returns the solution x.
        """
        n = len(self.b)
        I = np.eye(n)
        H = np.array([[1, 1], [1, -1]]) / np.sqrt(2)

        # Construct the matrix U
        A_t = expm(-1j * self.t * self.A)
        U = np.kron(I, H) @ np.kron(A_t, I) @ np.kron(I, H)

        # Construct the matrix b_tilde
        b_norm = np.linalg.norm(self.b)
        b_tilde = np.concatenate((self.b, np.zeros(n)), axis=None)
        b_tilde = b_tilde / b_norm

        # Solve the linear system of equations
        psi = np.linalg.solve(U, b_tilde)
        x = psi[:n] / psi[n]

        # calculate expected value
        y = np.linalg.solve(self.A, self.b)

        return x, y


if __name__ == '__main__':
    A = np.array([[3, 1], [1, 2]])
    b = np.array([1, 0])
    t = 0.5
    M = 2

    hhl = HHLAlgorithm(A, b, t, M)
    x = hhl.solve()

    print('The solution to Ax = b is:', x)

    A = np.array([[1, 1], [2, 3]])
    b = np.array([2, 5])
    t = 0.5
    M = 2

    hhl = HHLAlgorithm(A, b, t, M)
    x = hhl.solve()
    print('The solution to Ax = b is:', x)
