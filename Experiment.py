import numpy as np
import GateGeneration as gg
import MathematicaScripts as scripts

class Experiment:
    def __init__(self, num_qubits):
        """Initialize an experiment with a given number of qubits, all starting in the |0⟩ state."""
        self.num_qubits = num_qubits
        self.state = np.zeros(2**num_qubits, dtype=complex)
        self.state[0] = 1  # Start in |00...0⟩ state

    def set_state(self, state):
        self.state = state

    def apply_gate(self, gate, qubit):
        """Apply a single-qubit gate to a specific qubit."""
        full_gate = self._expand_gate(gate, qubit)
        self.state = full_gate @ self.state

    def apply_global_gate(self, gate):
        """Apply a global state to the whole system"""
        self.state = gate @ self.state

    def measure(self):
        """Measure the quantum state in the computational basis."""
        probabilities = np.abs(self.state)**2
        outcome = np.random.choice(len(probabilities), p=probabilities)
        return bin(outcome)[2:].zfill(self.num_qubits)

    def print_state(self):
        """Print the quantum state as a superposition of computational basis states."""
        print("Quantum State:")
        for i, amplitude in enumerate(self.state):
            if np.abs(amplitude) > 1e-6:  # Only show significant amplitudes
                basis_state = bin(i)[2:].zfill(self.num_qubits)
                print(f"|{basis_state}⟩: {amplitude.real:.4f} + {amplitude.imag:.4f}i")


    def global_rotation(self, theta, rotation, NIons):
        """Apply a global rotation to all qubits around the bloch sphere"""
        I = np.eye(2)
        sx = np.array([[0, 1], [1, 0]])
        sy = np.array([[0, -1j], [1j, 0]])
        sz = np.array([[1, 0], [0, -1]]) 

        #different types of rotations
        def rotation_x(self, theta):
            """Rotation about X-axis"""
            return np.cos(theta/2) * I - 1j * np.sin(theta/2) * sx
    
        def rotation_y(self, theta):
            """Rotation about Y-axis"""
            return np.cos(theta/2) * I - 1j * np.sin(theta/2) * sy
    
        def rotation_z(self, theta):
            """Rotation about Z-axis"""
            return np.array([[np.exp(-1j*theta/2), 0], 
                            [0, np.exp(1j*theta/2)]])
        
        #generate matrix operator
        matrix = np.eye(2**NIons)
        if rotation == "x":
            operator = rotation_x(self, theta)
        elif rotation == "y":
            operator = rotation_y(self, theta)
        elif rotation == "z":
            operator = rotation_z(self, theta)
        else:
            print("invalid rotation direction: use x, y, or z")
            return
        
        for i in range(NIons):
            local_rotation = np.array([1])
            for j in range(NIons):
                if i == j:
                    local_rotation = np.kron(local_rotation, operator)
                else:
                    local_rotation = np.kron(local_rotation, I)
            matrix = matrix @ local_rotation

        self.apply_global_gate(matrix)
    

    def global_entanglement(self, theta, Jij, rotation, NIons):
        I = np.eye(2)
        sx = np.array([[0, 1], [1, 0]])
        sy = np.array([[0, -1j], [1j, 0]])
        sz = np.array([[1, 0], [0, -1]]) 
        def entanglement_x(self, Jij, NIons):
            return gg.MS_Hamiltonian(Jij, N=NIons, phase=0)

        def entanglement_y(self, Jij, NIons):
            return gg.MS_Hamiltonian(Jij, NIons, np.pi)

        def entanglement_z(self, Jij, NIons):
            self.global_rotation(-np.pi/2, 'y', NIons)
            self.global_entanglement(theta, Jij, 'x', NIons)
            self.global_rotation(np.pi/2, 'y', NIons)
        if rotation == "x":
            matrix = entanglement_x(self, Jij, NIons)
        elif rotation == "y":
            matrix = entanglement_y(self, Jij, NIons)
        elif rotation == "z":
            entanglement_z(self, Jij, NIons)
            return
        else:
            print("invalid rotation direction: use x, y, or z")
            return
        matrix = matrix / np.linalg.norm(matrix)
        effective_matrix = np.cos(theta/2) * np.eye(2**NIons) - 1j * np.sin(theta/2) * matrix 
        # print("matrix:")
        # gg.print_matrix(matrix)
        # print("effective matrix:")
        # gg.print_matrix(effective_matrix)
        self.apply_global_gate(effective_matrix)
        return

