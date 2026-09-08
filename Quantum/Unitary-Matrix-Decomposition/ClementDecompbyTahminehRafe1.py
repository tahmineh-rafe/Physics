"""
Unitary Matrix Decomposition for Photonic Quantum Computing

Author:
    Tahmineh Rafe
    2026

Description:
    This code implements the decomposition of arbitrary unitary matrices
    into a sequence of two-mode transformations using the methods proposed
    by Reck et al. and Clements et al.

    The extracted parameters (theta and phi) define the configuration of
    Mach-Zehnder interferometers (MZIs) required to reconstruct the target
    unitary transformation in a programmable photonic circuit.

    The implementation, matrix decomposition procedure, reconstruction,
    and numerical verification were developed by Tahmineh Rafe.

References:
    [1] Reck, M., Zeilinger, A., Bernstein, H. J., & Bertani, P. (1994).
        Experimental realization of any discrete unitary operator.
        Physical Review Letters, 73(1), 58.

    [2] Clements, W. R., Humphreys, P. C., Metcalf, B. J.,
        Kolthammer, W. S., & Walmsley, I. A. (2016).
        Optimal design for universal multiport interferometers.
        Optica, 3(12), 1460-1465.

    [3] Capmany, J., & Pérez, D. (2020).
        Programmable integrated photonics.
        Oxford University Press.

    [4] Arrazola, J. M., Bergholm, V., Brádler, K., Bromley, T. R.,
        Collins, M. J., Dhand, I., ... & Zhang, Y. (2021).
        Quantum circuits with many photons on a programmable nanophotonic chip.
        Nature, 591(7848), 54-60.

Citation:
    If this implementation is used in academic work, please cite:

    Rafe, Tahmineh. (2026).
    Unitary Matrix Decomposition for Photonic Quantum Computing.
    GitHub repository: Physics.
"""

# Unit block:
# T_mn(theta, phi) = i exp(i theta/2) [[exp(i phi) sin(theta), cos(theta)],
#                                      [exp(i phi) cos(theta), -sin(theta)]]

import os
import numpy as np
import matplotlib.pyplot as plt

if os.getcwd()[-9:] == "notebooks":
    os.chdir("..")


# ==========================
# 1. Utility functions
# ==========================

def unitary_error(U):
    """
    Compute unitarity error:
        ||U†U - I||F
    """
    U = np.asarray(U, dtype=complex)
    N = U.shape[0]
    return np.linalg.norm(U.conj().T @ U - np.eye(N), ord="fro")


def project_to_unitary(A):
    """
    Project an almost-unitary matrix to the nearest unitary matrix
    using polar decomposition through SVD.

    This is useful when the input matrix is rounded.
    """
    W, S, Vh = np.linalg.svd(A)
    return W @ Vh


def random_unitary(N, seed=None):
    """
    Generate a random N x N unitary matrix using QR decomposition.
    This replaces itf.random_unitary(N).
    """
    rng = np.random.default_rng(seed)

    Z = rng.normal(size=(N, N)) + 1j * rng.normal(size=(N, N))
    Q, R = np.linalg.qr(Z)

    phases = np.diag(R) / np.abs(np.diag(R))
    U = Q * phases

    return U




def chapter4_T_block(N, m, n, theta, phi):
    """
    Build the full N x N T_mn(theta, phi) matrix based on Chapter 4.

    The 2 x 2 block is:

        T(theta, phi) =
        i exp(i theta/2)
        [[exp(i phi) sin(theta),  cos(theta)],
         [exp(i phi) cos(theta), -sin(theta)]]

    This block acts only on modes m and n.
    Python indexing starts from 0.
    """

    T = np.eye(N, dtype=complex)

    s = np.sin(theta)
    c = np.cos(theta)

    global_phase = 1j * np.exp(1j * theta / 2)
    exp_i_phi = np.exp(1j * phi)

    T[m, m] = global_phase * exp_i_phi * s
    T[m, n] = global_phase * c
    T[n, m] = global_phase * exp_i_phi * c
    T[n, n] = global_phase * (-s)

    return T


def chapter4_T_inverse_block(N, m, n, theta, phi):
    """
    Since T is unitary:

        T^{-1} = T†

    This is also consistent with Chapter 4.
    """

    T = chapter4_T_block(N, m, n, theta, phi)
    return T.conj().T


# ==========================
# 3. Clements zeroing rules
# ==========================

def get_clements_right_zeroing_angles(a, b, tol=1e-12):
    """
    Compute theta and phi for right-side zeroing.

    For right multiplication:

        U_new = U @ T^{-1}

    To zero the target element using Chapter 4 condition:

        a exp(-i phi) sin(theta) + b cos(theta) = 0

    Therefore:

        tan(theta) = |b| / |a|

        exp(-i phi) = - (b / a) cot(theta)

    Special cases are handled to avoid division by zero.
    """

    if abs(a) < tol:
        theta = np.pi / 2
        phi = 0.0
        return theta, phi

    if abs(b) < tol:
        theta = 0.0
        phi = 0.0
        return theta, phi

    theta = np.arctan2(abs(b), abs(a))

    exp_minus_i_phi = -(b / a) / np.tan(theta)
    phi = -np.angle(exp_minus_i_phi)

    return theta, phi


def get_clements_left_zeroing_angles(a, b, tol=1e-12):
    """
    Compute theta and phi for left-side zeroing.

    For left multiplication:

        U_new = T @ U

    The affected rows are mixed. Using the Chapter 4 block, the lower
    transformed element is zeroed by:

        exp(i phi) cos(theta) a - sin(theta) b = 0

    Therefore:

        tan(theta) = |a| / |b|

        exp(i phi) = (b / a) tan(theta)

    Special cases are handled to avoid division by zero.
    """

    if abs(a) < tol:
        theta = 0.0
        phi = 0.0
        return theta, phi

    if abs(b) < tol:
        theta = np.pi / 2
        phi = 0.0
        return theta, phi

    theta = np.arctan2(abs(a), abs(b))

    exp_i_phi = (b / a) * np.tan(theta)
    phi = np.angle(exp_i_phi)

    return theta, phi


# ==========================
# 4. Clements decomposition
# ==========================

def chapter4_clements_decomposition(U, tol=1e-12):
    """
    Clements rectangular decomposition using the Chapter 4 MZI block.

    The algorithm alternates between right-side and left-side multiplications.

    For odd steps:

        U_new = U @ T^{-1}

    For even steps:

        U_new = T @ U

    At the end:

        D = L_total @ U @ R_total

    Therefore:

        U = L_total† @ D @ R_total†

    Outputs:
        D                 final diagonal phase matrix
        operations        list of MZI blocks and parameters
        U_reconstructed   reconstructed matrix
        A_final           final diagonalized matrix
        error             Frobenius reconstruction error
        fidelity          normalized matrix fidelity
    """

    U = np.asarray(U, dtype=complex)

    if U.shape[0] != U.shape[1]:
        raise ValueError("Input matrix must be square.")

    N = U.shape[0]
    A = U.copy()

    L_total = np.eye(N, dtype=complex)
    R_total = np.eye(N, dtype=complex)

    operations = []

    for i in range(1, N):

        if i % 2 == 1:
            # Odd step: right-side multiplication
            for j in range(0, i):

                row = N - 1 - j
                m = i - j - 1
                n = m + 1

                a = A[row, m]
                b = A[row, n]

                theta, phi = get_clements_right_zeroing_angles(a, b, tol=tol)

                T = chapter4_T_block(N, m, n, theta, phi)
                T_inv = chapter4_T_inverse_block(N, m, n, theta, phi)

                A = A @ T_inv
                R_total = R_total @ T_inv

                if abs(A[row, m]) < 1e-10:
                    A[row, m] = 0.0 + 0.0j

                operations.append({
                    "side": "right",
                    "m": m,
                    "n": n,
                    "theta": theta,
                    "phi": phi,
                    "T": T,
                    "T_inverse": T_inv
                })

        else:
            # Even step: left-side multiplication
            for j in range(1, i + 1):

                n = N + j - i - 1
                m = n - 1
                col = j - 1

                a = A[m, col]
                b = A[n, col]

                theta, phi = get_clements_left_zeroing_angles(a, b, tol=tol)

                T = chapter4_T_block(N, m, n, theta, phi)
                T_inv = chapter4_T_inverse_block(N, m, n, theta, phi)

                A = T @ A
                L_total = T @ L_total

                if abs(A[n, col]) < 1e-10:
                    A[n, col] = 0.0 + 0.0j

                operations.append({
                    "side": "left",
                    "m": m,
                    "n": n,
                    "theta": theta,
                    "phi": phi,
                    "T": T,
                    "T_inverse": T_inv
                })

    D = np.diag(np.diag(A))

    U_reconstructed = L_total.conj().T @ D @ R_total.conj().T

    error = np.linalg.norm(U - U_reconstructed, ord="fro")
    fidelity = abs(np.trace(U.conj().T @ U_reconstructed)) / N

    return {
        "D": D,
        "operations": operations,
        "U_reconstructed": U_reconstructed,
        "A_final": A,
        "L_total": L_total,
        "R_total": R_total,
        "error": error,
        "fidelity": fidelity
    }


# ==========================
# 5. Simple Clements circuit drawing
# ==========================

def draw_clements_network(operations, N, title="Clements rectangular decomposition"):
    """
    Draw a simple schematic of the Clements network.
    This is not a physical layout generator; it is a clear visual map
    of which modes are coupled at each step.
    """

    fig, ax = plt.subplots(figsize=(1.2 * len(operations) + 2, 0.8 * N + 2))

    # Draw mode lines
    for mode in range(N):
        ax.hlines(y=mode, xmin=0, xmax=len(operations) + 1, linewidth=1)
        ax.text(-0.4, mode, f"Mode {mode}", va="center", ha="right")

    # Draw MZI blocks
    for k, op in enumerate(operations, start=1):
        m = op["m"]
        n = op["n"]

        ax.plot([k, k], [m, n], linewidth=2)
        ax.scatter([k, k], [m, n], s=120)

        label = f"{op['side']}\nθ={op['theta']:.2f}\nφ={op['phi']:.2f}"
        ax.text(k, (m + n) / 2, label, ha="center", va="center", fontsize=8)

    ax.set_title(title)
    ax.set_xlim(-1, len(operations) + 1)
    ax.set_ylim(N - 0.5, -0.5)
    ax.axis("off")

    plt.show()


# ==========================
# 6. Example: use your own matrix
# ==========================

U_raw = np.array([
    [-0.3250-0.4190j,  -0.3362-0.1772j,  -0.3803+0.1376j,   0.3149+0.5582j],
    [-0.6316-0.03660j,  0.4966-0.4024j,  -0.2632-0.0127j,  -0.1553-0.31230j],
    [-0.1518+0.0883j,  -0.5153+0.0818j,  -0.3063-0.7007j,  -0.0142-0.3346j],
    [ 0.3766-0.3819j,  -0.2557-0.3306j,  -0.1885+0.3830j,   0.0749-0.5915j]
], dtype=complex)

N = U_raw.shape[0]

print("Original unitarity error:")
print(unitary_error(U_raw))

# Because the matrix is rounded, project it to the nearest unitary matrix.
U = project_to_unitary(U_raw)

print("\nProjected unitarity error:")
print(unitary_error(U))

result = chapter4_clements_decomposition(U)

D = result["D"]
operations = result["operations"]
U_rec = result["U_reconstructed"]
A_final = result["A_final"]
error = result["error"]
fidelity = result["fidelity"]

print("\nFinal matrix after Clements zeroing:")
print(np.round(A_final, 6))

print("\nFinal phase matrix D:")
print(np.round(D, 6))

print("\nReconstructed matrix:")
print(np.round(U_rec, 6))

print("\nReconstruction error ||U - U_rec||F:")
print(error)

print("\nMatrix fidelity:")
print(fidelity)

print("\nNumber of Clements MZI blocks:")
print(len(operations))

print("\nExpected number N(N-1)/2:")
print(N * (N - 1) // 2)

print("\nClements parameters based on Chapter 4 block:")
for i, op in enumerate(operations, start=1):
    print(
        f"{i}: side = {op['side']}, "
        f"modes ({op['m']}, {op['n']}), "
        f"theta = {op['theta']:.10f}, "
        f"phi = {op['phi']:.10f}"
    )

draw_clements_network(operations, N, title="Clements decomposition based on Chapter 4 MZI block")
