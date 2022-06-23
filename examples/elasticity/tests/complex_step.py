"""
This script generate results given random inputs for expressions in a2dtmp.h and
a2dtmp2d.h
"""

import numpy as np

def print_flatten_8(mat, end='\n'):
    s = ''
    for v in mat.flatten():
        s +=f'{v:.8f}, '
    print(s[:-2], end=end)

def print_flatten_16(mat, end='\n'):
    s = ''
    for v in mat.flatten():
        s +=f'{v:.16f}, '
    print(s[:-2], end=end)

def MatMat_forward(func, A, B, Ab, Bb):
    """
    Compute dC given dA and dB for the expression C = f(A, B), where
    C can be either matrix or scalar
    """
    h = 1e-30
    Cb = func(A + 1j*h*Ab, B + 1j*h*Bb).imag / h
    return Cb

def Mat_forward(func, A, Ab):
    """
    Compute dC given dA for the expression C = f(A), where
    C can be either matrix or scalar
    """
    h = 1e-30
    Cb = func(A + 1j*h*Ab).imag / h
    return Cb

def MatMat_reverse(func, A, B, Cb):
    """
    Compute dA and dB given dC for expression C = f(A, B), where C is matrix
    """
    h = 1e-30

    if not isinstance(Cb, np.ndarray):
        Cb = np.array(Cb).reshape(1, 1)

    dCdA = np.zeros((*Cb.shape, *A.shape))
    dCdB = np.zeros((*Cb.shape, *B.shape))

    Ab = np.zeros(A.shape)
    Bb = np.zeros(B.shape)

    # Compute derivative of function w.r.t. A and B, note that results
    # here are 4-th order tensors
    for k in range(A.shape[0]):
        for l in range(A.shape[1]):
            A_ = np.copy(A).astype(complex)
            B_ = np.copy(B).astype(complex)
            A_[k, l] += 1j * h
            B_[k, l] += 1j * h
            dCdA[:, :, k, l] = func(A_, B).imag / h
            dCdB[:, :, k, l] = func(A, B_).imag / h

    # Compute bar(A) and bar(B)
    for k in range(A.shape[0]):
        for l in range(A.shape[1]):
            Ab[k, l] = np.sum(Cb*dCdA[:, :, k, l])
            Bb[k, l] = np.sum(Cb*dCdB[:, :, k, l])

    return Ab, Bb

def Mat_reverse(func, A, Cb):
    """
    Compute dA given dC for expression C = f(A)
    """
    h = 1e-30

    if not isinstance(Cb, np.ndarray):
        Cb = np.array(Cb).reshape(1, 1)

    dCdA = np.zeros((*Cb.shape, *A.shape))
    Ab = np.zeros(A.shape)

    # Compute derivative of function w.r.t. A
    for k in range(A.shape[0]):
        for l in range(A.shape[1]):
            A_ = np.copy(A).astype(complex)
            A_[k, l] += 1j * h
            dCdA[:, :, k, l] = func(A_).imag / h

    # Compute bar(A)
    for k in range(A.shape[0]):
        for l in range(A.shape[1]):
            Ab[k, l] = np.sum(Cb*dCdA[:, :, k, l])

    return Ab

def _compute_hess(func, x):
    """
    Compute Hessian matrix given function and input vector x
    """
    hc = 1e-15
    h = 1e-8
    N = len(x)

    test = func(x)
    if isinstance(test, np.ndarray):
        H = np.zeros((*test.shape, N, N))
    else:
        H = np.zeros((N, N))

    for i in range(N):
        x1 = x.copy().astype(complex)
        x1[i] += 1j * hc
        for j in range(N):
            x2 = x1.copy().astype(complex)
            x2[j] += h
            x3 = x1.copy().astype(complex)
            x3[j] -= h
            u2 = func(x2)
            u3 = func(x3)
            H[..., i, j] = (u2 - u3).imag / (2 * h * hc)
    return H


def MatMat_hreverse(func, A, Ap, B, Bp, Cb, Ch):
    """
    Compute the Hessian-vector product for the following binary expression:
        C = func(A, B)
    """

    if not isinstance(Cb, np.ndarray):
        Cb = np.array(Cb).reshape(1, 1)
        Ch = np.array(Ch).reshape(1, 1)
        _func = lambda A, B: np.array(func(A, B)).reshape(1, 1)
    else:
        _func = func

    # Compute dCdA and dCdB
    def _compute_grad(A, B):
        hc = 1e-15
        dCdA = np.zeros((*Cb.shape, *A.shape))
        dCdB = np.zeros((*Cb.shape, *B.shape))
        for k in range(A.shape[0]):
            for l in range(A.shape[1]):
                A_ = np.copy(A).astype(complex)
                B_ = np.copy(B).astype(complex)
                A_[k, l] += 1j * hc
                B_[k, l] += 1j * hc
                dCdA[:, :, k, l] = _func(A_, B).imag / hc
                dCdB[:, :, k, l] = _func(A, B_).imag / hc
        return dCdA, dCdB

    dCdA, dCdB = _compute_grad(A, B)

    def fun_hess(x):
        _A = x[:np.prod(A.shape)].reshape(A.shape)
        _B = x[np.prod(A.shape):].reshape(B.shape)
        return _func(_A, _B)

    x = np.concatenate([A.flatten(), B.flatten()])
    H = _compute_hess(fun_hess, x)
    nA = np.prod(A.shape)
    d2CdA2 = H[:, :, :nA, :nA]
    d2CdAB = H[:, :, :nA, nA:]
    d2CdBA = H[:, :, nA:, :nA]
    d2CdB2 = H[:, :, nA:, nA:]

    # Compute Ah and Bh
    Ah = np.einsum("ij,ijkl -> kl", Ch, dCdA)
    At = np.einsum("ijkl,l -> ijk", d2CdA2, Ap.flatten())
    At += np.einsum("ijkl,l -> ijk", d2CdAB, Bp.flatten())
    Ah += np.einsum("ij,ijk -> k", Cb, At).reshape(A.shape)

    Bh = np.einsum("ij, ijkl -> kl", Ch, dCdB)
    Bt = np.einsum("ijkl,l -> ijk", d2CdB2, Bp.flatten())
    Bt += np.einsum("ijkl,l -> ijk", d2CdBA, Ap.flatten())
    Bh += np.einsum("ij,ijk -> k", Cb, Bt).reshape(B.shape)

    # return Ah, Bh, dCdA, dCdB, d2CdA2, d2CdB2, d2CdAB, d2CdBA
    return Ah, Bh

def test_grad_hess():
    np.random.seed(0)
    H1 = np.random.rand(8, 8)
    H2 = np.random.rand(8, 8)
    H3 = np.random.rand(8, 8)
    H4 = np.random.rand(8, 8)
    b = np.random.rand(8)
    x = np.random.rand(8)
    H1 = H1 + H1.T
    H2 = H2 + H2.T
    H3 = H3 + H3.T
    H4 = H4 + H4.T
    poly = lambda x: np.array([
        0.5 * x.dot(H1.dot(x)),
        0.5 * x.dot(H2.dot(x)),
        0.5 * x.dot(H3.dot(x)),
        0.5 * x.dot(H4.dot(x)),
        ]).reshape(2,2)
    func = lambda A, B: poly(np.concatenate([A.flatten(), B.flatten()]))
    A = x[:4].reshape(2,2)
    B = x[4:].reshape(2,2)
    Ah, Bh, dCdA, dCdB, d2CdA2, d2CdB2, d2CdAB, d2CdBA = MatMat_hreverse(func, A=A, Ap=A, B=B, Bp=B, Cb=np.zeros((2,2)), Ch=np.zeros((2,2)))


    print(np.abs(H1.dot(x)[:4] - dCdA[0, 0, ...].flatten()).max())
    print(np.abs(H1.dot(x)[4:] - dCdB[0, 0, ...].flatten()).max())
    print(np.abs(H2.dot(x)[:4] - dCdA[0, 1, ...].flatten()).max())
    print(np.abs(H2.dot(x)[4:] - dCdB[0, 1, ...].flatten()).max())
    print(np.abs(H3.dot(x)[:4] - dCdA[1, 0, ...].flatten()).max())
    print(np.abs(H3.dot(x)[4:] - dCdB[1, 0, ...].flatten()).max())
    print(np.abs(H4.dot(x)[:4] - dCdA[1, 1, ...].flatten()).max())
    print(np.abs(H4.dot(x)[4:] - dCdB[1, 1, ...].flatten()).max())

    print(np.abs(H1[:4,:4] - d2CdA2[0, 0, ...]).max())
    print(np.abs(H1[:4,4:] - d2CdAB[0, 0, ...]).max())
    print(np.abs(H1[4:,:4] - d2CdBA[0, 0, ...]).max())
    print(np.abs(H1[4:,4:] - d2CdB2[0, 0, ...]).max())

    print(np.abs(H2[:4,:4] - d2CdA2[0, 1, ...]).max())
    print(np.abs(H2[:4,4:] - d2CdAB[0, 1, ...]).max())
    print(np.abs(H2[4:,:4] - d2CdBA[0, 1, ...]).max())
    print(np.abs(H2[4:,4:] - d2CdB2[0, 1, ...]).max())

    print(np.abs(H3[:4,:4] - d2CdA2[1, 0, ...]).max())
    print(np.abs(H3[:4,4:] - d2CdAB[1, 0, ...]).max())
    print(np.abs(H3[4:,:4] - d2CdBA[1, 0, ...]).max())
    print(np.abs(H3[4:,4:] - d2CdB2[1, 0, ...]).max())

    print(np.abs(H4[:4,:4] - d2CdA2[1, 1, ...]).max())
    print(np.abs(H4[:4,4:] - d2CdAB[1, 1, ...]).max())
    print(np.abs(H4[4:,:4] - d2CdBA[1, 1, ...]).max())
    return

def Mat_hreverse(func, A, Ap, Cb, Ch):
    """
    Compute the Hessian-vector product for the following unitary expression:
        C = func(A)
    """
    hc = 1e-30
    h = 1e-8

    if not isinstance(Cb, np.ndarray):
        Cb = np.array(Cb).reshape(1, 1)
        Ch = np.array(Ch).reshape(1, 1)

    # Compute dCdA
    dCdA = np.zeros((*Cb.shape, *A.shape))
    for k in range(A.shape[0]):
        for l in range(A.shape[1]):
            A_ = np.copy(A).astype(complex)
            A_[k, l] += 1j * hc
            dCdA[:, :, k, l] = func(A_).imag / hc

    # Compute d2C/dA2
    d2CdA2 = np.zeros((*Cb.shape, *A.shape, *A.shape))
    for i in range(A.shape[0]):
        for j in range(A.shape[1]):
            A1 = A.copy().astype(complex)
            A1[i, j] += 1j * hc
            for k in range(A.shape[0]):
                for l in range(A.shape[1]):
                    A2 = A1.copy().astype(complex)
                    A2[k, l] += h
                    A3 = A1.copy().astype(complex)
                    A3[k, l] -= h
                    u2 = func(A2)
                    u3 = func(A3)
                    d2CdA2[:, :, i, j, k, l] = (u2 - u3).imag / (2 * h * hc)

    # Compute Ah
    Ah = np.einsum("ij, ijkl -> kl", Ch, dCdA)
    At = np.einsum("ijklmn, mn->ijkl", d2CdA2, Ap)
    Ah += np.einsum("ij, ijkl -> kl", Cb, At)

    return Ah


def det_func(A):
    ndim = A.shape[0]
    if ndim == 3:
        det = A[0, 0] * A[1, 1] * A[2, 2] + A[0, 1] * A[1, 2] * A[2, 0] + A[0, 2] * A[1, 0] * A[2, 1] \
            - A[0, 2] * A[1, 1] * A[2, 0] - A[0, 1] * A[1, 0] * A[2, 2] - A[0, 0] * A[1, 2] * A[2, 1]
    else:
        det = A[0, 0] * A[1, 1] - A[1, 0] * A[0, 1]
    return det

def inv_func(A):
    ndim = A.shape[0]
    if ndim == 3:
        invA = np.zeros((3, 3), dtype=A.dtype)
        det = det_func(A)
        invA[0, 0] =  (A[1, 1] * A[2, 2] - A[1, 2] * A[2, 1]) / det
        invA[0, 1] = -(A[0, 1] * A[2, 2] - A[0, 2] * A[2, 1]) / det
        invA[0, 2] =  (A[0, 1] * A[1, 2] - A[0, 2] * A[1, 1]) / det
        invA[1, 0] = -(A[1, 0] * A[2, 2] - A[1, 2] * A[2, 0]) / det
        invA[1, 1] =  (A[0, 0] * A[2, 2] - A[0, 2] * A[2, 0]) / det
        invA[1, 2] = -(A[0, 0] * A[1, 2] - A[0, 2] * A[1, 0]) / det
        invA[2, 0] =  (A[1, 0] * A[2, 1] - A[1, 1] * A[2, 0]) / det
        invA[2, 1] = -(A[0, 0] * A[2, 1] - A[0, 1] * A[2, 0]) / det
        invA[2, 2] =  (A[0, 0] * A[1, 1] - A[0, 1] * A[1, 0]) / det
    else:
        invA = np.zeros((2, 2), dtype=A.dtype)
        det = det_func(A)
        invA[0, 0] = A[1, 1] / det
        invA[1, 0] = -A[1, 0] / det
        invA[0, 1] = -A[0, 1] / det
        invA[1, 1] = A[0, 0] / det
    return invA

def symm2mat(S):
    if len(S.flatten()) == 6:
        ndim = 3
    elif len(S.flatten()) == 3:
        ndim = 2
    else:
        raise ValueError(f"Unknown ndim({ndim}) given S({S})")
    Svals = S.flatten()
    M = np.zeros((ndim, ndim), dtype=S.dtype)
    for i in range(ndim):
        for j in range(ndim):
            if i >= j:
                M[i, j] = Svals[j + i * (i + 1) // 2]
            else:
                M[i, j] = Svals[i + j * (j + 1) // 2]
    return M

def mat2symm(Mat):
    N = Mat.shape[0]
    temp = np.tril(Mat).flatten()
    S = temp[temp != 0]
    S = S.reshape(-1, 1)
    return S

def SymmTrace(Ssymm):
    return symm2mat(Ssymm).trace()

def SymmMultTrace(Asymm, Bsymm):
    A = symm2mat(Asymm)
    B = symm2mat(Bsymm)
    return A.dot(B).diagonal().sum()

def SymmIsoConstitutive(mu, lam, Esymm):
    E = symm2mat(Esymm)
    S = 2 * mu * E + lam * E.trace() * np.eye(E.shape[0])
    return mat2symm(S)

def SymmIsoEnergy(mu, lam, Esymm):
    S = SymmIsoConstitutive(mu, lam, Esymm)
    S = symm2mat(S)
    E = symm2mat(Esymm)
    return 0.5 * np.dot(E, S).trace()

def AugSymmIsoEnergy(AugE):
    Esymm = AugE[0:-2]
    mu = AugE[-2, 0]
    lam = AugE[-1, 0]
    return SymmIsoEnergy(mu, lam, Esymm)

def GreenStrain(Ux):
    E = 0.5 * (Ux + Ux.T + Ux.T.dot(Ux))
    return mat2symm(E)

def LinearGreenStrain(Ux):
    E = 0.5 * (Ux + Ux.T)
    return mat2symm(E)

def generate_reference_values(ndim, A, B, Ab, Bb, Cb, Ch, S, E, Sb, Eb, Sh, sb,
                              sh, mu, lam, mub, lamb, AugE, AugEb):
    # ADMatxADMat
    print("\nADMatxADMat\n")
    cases = ["AB", "ATB", "ABT", "ATBT"]
    apply_to_A = [lambda x:x.copy(), lambda x:x.T.copy(), lambda x:x.copy(), lambda x:x.T.copy()]
    apply_to_B = [lambda x:x.copy(), lambda x:x.copy(), lambda x:x.T.copy(), lambda x:x.T.copy()]
    for case, to_A, to_B in zip(cases, apply_to_A, apply_to_B):
        _A = to_A(A)
        _B = to_B(B)
        _Ab = to_A(Ab)
        _Bb = to_B(Bb)
        _Cb = MatMat_forward(np.dot, _A, _B, _Ab, _Bb)
        _Ab, _Bb = MatMat_reverse(np.dot, _A, _B, _Cb)
        print(f"\n//{case}")
        print(f"const T {case}_dC_out[{ndim * ndim}] = {{", end='')
        print_flatten_16(_Cb, end='')
        print("};")
        print(f"const T {case}_dA_out[{ndim * ndim}] = {{", end='')
        print_flatten_16(to_A(_Ab), end='')
        print("};")
        print(f"const T {case}_dB_out[{ndim * ndim}] = {{", end='')
        print_flatten_16(to_B(_Bb), end='')
        print("};")

    # A2DMatxA2DMat
    print("\nA2DMatxA2DMat\n")
    for case, to_A, to_B in zip(cases, apply_to_A, apply_to_B):
        _A = to_A(A)
        _B = to_B(B)
        _Ap = to_A(Ab)
        _Bp = to_B(Bb)
        _Cb = Cb.copy()
        _Ch = Ch.copy()
        _Cp = MatMat_forward(np.dot, _A, _B, _Ap, _Bp)
        _Ab, _Bb = MatMat_reverse(np.dot, _A, _B, _Cb)
        _Ah, _Bh = MatMat_hreverse(np.dot, _A, _Ap, _B, _Bp, _Cb, _Ch)
        print(f"\n//{case}")
        print(f"const T {case}_Cp_out[{ndim * ndim}] = {{", end='')
        print_flatten_16(_Cp, end='')
        print("};")
        print(f"const T {case}_Ab_out[{ndim * ndim}] = {{", end='')
        print_flatten_16(to_A(_Ab), end='')
        print("};")
        print(f"const T {case}_Bb_out[{ndim * ndim}] = {{", end='')
        print_flatten_16(to_B(_Bb), end='')
        print("};")
        print(f"const T {case}_Ah_out[{ndim * ndim}] = {{", end='')
        print_flatten_16(to_A(_Ah), end='')
        print("};")
        print(f"const T {case}_Bh_out[{ndim * ndim}] = {{", end='')
        print_flatten_16(to_B(_Bh), end='')
        print("};")


    # ADMatxMat
    print("\nADMatxMat\n")
    for case, to_A, to_B in zip(cases, apply_to_A, apply_to_B):
        _A = to_A(A)
        _B = to_B(B)
        _Ab = to_A(Ab)
        _Cb = Mat_forward(lambda X: X.dot(_B), _A, _Ab)
        _Ab = Mat_reverse(lambda X: X.dot(_B), _A, _Cb)
        print(f"\n//{case}")
        print(f"const T {case}_dC_out[{ndim * ndim}] = {{", end='')
        print_flatten_16(_Cb, end='')
        print("};")
        print(f"const T {case}_dA_out[{ndim * ndim}] = {{", end='')
        print_flatten_16(to_A(_Ab), end='')
        print("};")

    # MatxADMat
    print("\nMatxADMat\n")
    for case, to_A, to_B in zip(cases, apply_to_A, apply_to_B):
        _A = to_A(A)
        _B = to_B(B)
        _Bb = to_B(Bb)
        _Cb = Mat_forward(lambda X: _A.dot(X), _B, _Bb)
        _Bb = Mat_reverse(lambda X: _A.dot(X), _B, _Cb)
        print(f"\n//{case}")
        print(f"const T {case}_dC_out[{ndim * ndim}] = {{", end='')
        print_flatten_16(_Cb, end='')
        print("};")
        print(f"const T {case}_dB_out[{ndim * ndim}] = {{", end='')
        print_flatten_16(to_B(_Bb), end='')
        print("};")

    # A2DMatxMat
    print("\nA2DMatxMat\n")
    for case, to_A, to_B in zip(cases, apply_to_A, apply_to_B):
        _A = to_A(A)
        _B = to_B(B)
        _Ap = to_A(Ab)
        _Cb = Cb.copy()
        _Ch = Ch.copy()
        _Cp = Mat_forward(lambda A: A.dot(_B), _A, _Ap)
        _Ab = Mat_reverse(lambda A: A.dot(_B), _A, _Cb)
        _Ah = Mat_hreverse(lambda A: A.dot(_B), _A, _Ap, _Cb, _Ch)
        print(f"\n//{case}")
        print(f"const T {case}_Cp_out[{ndim * ndim}] = {{", end='')
        print_flatten_16(_Cp, end='')
        print("};")
        print(f"const T {case}_Ab_out[{ndim * ndim}] = {{", end='')
        print_flatten_16(to_A(_Ab), end='')
        print("};")
        print(f"const T {case}_Ah_out[{ndim * ndim}] = {{", end='')
        print_flatten_16(to_A(_Ah), end='')
        print("};")

    # MatxA2DMat
    print("\nMatxA2DMat\n")
    for case, to_A, to_B in zip(cases, apply_to_A, apply_to_B):
        _A = to_A(A)
        _B = to_B(B)
        _Bp = to_B(Bb)
        _Cb = Cb.copy()
        _Ch = Ch.copy()
        _Cp = Mat_forward(lambda B: _A.dot(B), _B, _Bp)
        _Bb = Mat_reverse(lambda B: _A.dot(B), _B, _Cb)
        _Bh = Mat_hreverse(lambda B: _A.dot(B), _B, _Bp, _Cb, _Ch)
        print(f"\n//{case}")
        print(f"const T {case}_Cp_out[{ndim * ndim}] = {{", end='')
        print_flatten_16(_Cp, end='')
        print("};")
        print(f"const T {case}_Bb_out[{ndim * ndim}] = {{", end='')
        print_flatten_16(to_B(_Bb), end='')
        print("};")
        print(f"const T {case}_Bh_out[{ndim * ndim}] = {{", end='')
        print_flatten_16(to_B(_Bh), end='')
        print("};")

    # DetA
    print("\nDetA\n")

    # DetA: AD
    case = "AD"
    _A = A.copy()
    _Ab = Ab.copy()
    _sb = Mat_forward(lambda x: det_func(x), _A, _Ab)
    _Ab = Mat_reverse(lambda x: det_func(x), _A, _sb)
    print(f"\n//{case}")
    print(f"const T {case}_sb_out = {_sb:.16f};")
    print(f"const T {case}_Ab_out[{ndim * ndim}] = {{", end='')
    print_flatten_16(_Ab, end='')
    print("};")

    # DetA: A2D
    case = "A2D"
    _A = A.copy()
    _Ap = Ab.copy()
    _sb = sb
    _sh = sh
    _sp = Mat_forward(lambda x: det_func(x), _A, _Ap)
    _Ab = Mat_reverse(lambda x: det_func(x), _A, _sb)
    _Ah = Mat_hreverse(lambda x: det_func(x), _A, _Ap, _sb, _sh)
    print(f"\n//{case}")
    print(f"const T {case}_sp_out = {_sp:.16f};")
    print(f"const T {case}_Ab_out[{ndim * ndim}] = {{", end='')
    print_flatten_16(_Ab, end='')
    print("};")
    print(f"const T {case}_Ah_out[{ndim * ndim}] = {{", end='')
    print_flatten_16(_Ah, end='')
    print("};")

    # InvA
    print("\nInvA\n")

    # DetA: AD
    case = "AD"
    _A = A.copy()
    _Ab = Ab.copy()
    _Cb = Mat_forward(lambda x: inv_func(x), _A, _Ab)
    _Ab = Mat_reverse(lambda x: inv_func(x), _A, _Cb)
    print(f"\n//{case}")
    print(f"const T {case}_invAb_out[{ndim * ndim}] = {{", end='')
    print_flatten_16(_Cb, end='')
    print("};")
    print(f"const T {case}_Ab_out[{ndim * ndim}] = {{", end='')
    print_flatten_16(_Ab, end='')
    print("};")

    # DetA: A2D
    case = "A2D"
    _A = A.copy()
    _Ap = Ab.copy()
    _Cb = Cb.copy()
    _Ch = Ch.copy()
    _Cp = Mat_forward(lambda x: inv_func(x), _A, _Ap)
    _Ab = Mat_reverse(lambda x: inv_func(x), _A, _Cb)
    _Ah = Mat_hreverse(lambda x: inv_func(x), _A, _Ap, _Cb, _Ch)
    print(f"\n//{case}")
    print(f"const T {case}_Ab_out[{ndim * ndim}] = {{", end='')
    print_flatten_16(_Ab, end='')
    print("};")
    print(f"const T {case}_Ah_out[{ndim * ndim}] = {{", end='')
    print_flatten_16(_Ah, end='')
    print("};")
    print(f"const T {case}_Cp_out[{ndim * ndim}] = {{", end='')
    print_flatten_16(_Cp, end='')
    print("};")

    # SymmTrace
    print("\nSymmTrace\n")

    # SymmTrace: A
    _S = S.copy()
    _tr = SymmTrace(_S)
    print(f"const T tr_out = {_tr:.16f};")

    # SymmTrace: AD
    case = "AD"
    _S = S.copy()
    _Sb = Sb.copy()
    _trb = Mat_forward(SymmTrace, _S, _Sb)
    _Sb = Mat_reverse(SymmTrace, _S, _trb)
    print(f"\n//{case}")
    print(f"const T {case}_trb_out = {_trb:.16f};")
    print(f"const T {case}_Sb_out[{ndim * ndim}] = {{", end='')
    print_flatten_16(symm2mat(_Sb), end='')
    print("};")

    # SymmTrace: A2D
    case = "A2D"
    _S = S.copy()
    _Sp = Sb.copy()
    _sb = sb
    _sh = sh
    _sp = Mat_forward(SymmTrace, _S, _Sp)
    _Sb = Mat_reverse(SymmTrace, _S, _sb)
    _Sh = Mat_hreverse(SymmTrace, _S, _Sp, _sb, _sh)
    print(f"\n//{case}")
    print(f"const T {case}_Sb_out[{ndim * ndim}] = {{", end='')
    print_flatten_16(symm2mat(_Sb), end='')
    print("};")
    print(f"const T {case}_Sh_out[{ndim * ndim}] = {{", end='')
    print_flatten_16(symm2mat(_Sh), end='')
    print("};")
    print(f"const T {case}_trp_out = {_sp:.16f};")

    # SymmMultTrace
    print("\nSymmMultTrace\n")

    # SymmMultTrace: AD
    case = "AD"
    _S = S.copy()
    _E = E.copy()
    _Sb = Sb.copy()
    _Eb = Eb.copy()
    _Cb = MatMat_forward(SymmMultTrace, _S, _E, _Sb, _Eb)
    _Sb, _Eb = MatMat_reverse(SymmMultTrace, _S, _E, _Cb)
    print(f"\n//{case}")
    print(f"const T {case}_Cb_out = {_Cb:.16f};")
    print(f"const T {case}_Sb_out[{ndim * ndim}] = {{", end='')
    print_flatten_16(symm2mat(_Sb), end='')
    print("};")
    print(f"const T {case}_Eb_out[{ndim * ndim}] = {{", end='')
    print_flatten_16(symm2mat(_Eb), end='')
    print("};")

    # SymmMultTrace: A2D
    case = "A2D"
    _S = S.copy()
    _E = E.copy()
    _Sp = Sb.copy()
    _Ep = Eb.copy()
    _sb = sb
    _sh = sh
    _sp = MatMat_forward(SymmMultTrace, _S, _E, _Sp, _Ep)
    _Sb, _Eb = MatMat_reverse(SymmMultTrace, _S, _E, _sb)
    _Sh, _Eh = MatMat_hreverse(SymmMultTrace, _S, _Sp, _E, _Ep, _sb, _sh)
    print(f"\n//{case}")
    print(f"const T {case}_Sb_out[{ndim * ndim}] = {{", end='')
    print_flatten_16(symm2mat(_Sb), end='')
    print("};")
    print(f"const T {case}_Sh_out[{ndim * ndim}] = {{", end='')
    print_flatten_16(symm2mat(_Sh), end='')
    print("};")
    print(f"const T {case}_Eb_out[{ndim * ndim}] = {{", end='')
    print_flatten_16(symm2mat(_Eb), end='')
    print("};")
    print(f"const T {case}_Eh_out[{ndim * ndim}] = {{", end='')
    print_flatten_16(symm2mat(_Eh), end='')
    print("};")
    print(f"const T {case}_sp_out = {_sp:.16f};")

    # SymmIsoConstitutive
    print("\nSymmIsoConstitutive\n")

    # SymmIsoConstitutive: A
    _E = E
    _S = SymmIsoConstitutive(mu, lam, _E)
    print(f"const T S_out[{ndim * ndim}] = {{", end='')
    print_flatten_16(symm2mat(_S), end='')
    print("};")

    # SymmIsoConstitutive: AD
    case = "AD"
    _E = E.copy()
    _Eb = Eb.copy()
    _Sb = Mat_forward(lambda x: SymmIsoConstitutive(mu, lam, x), _E, _Eb)
    _Eb = Mat_reverse(lambda x: SymmIsoConstitutive(mu, lam, x), _E, _Sb)
    print(f"\n//{case}")
    print(f"const T {case}_Sb_out[{ndim * ndim}] = {{", end='')
    print_flatten_16(symm2mat(_Sb), end='')
    print("};")
    print(f"const T {case}_Eb_out[{ndim * ndim}] = {{", end='')
    print_flatten_16(symm2mat(_Eb), end='')
    print("};")

    # SymmIsoConstitutive: A2D
    case = "A2D"
    _E = E.copy()
    _Ep = Eb.copy()
    _Sb = Sb.copy()
    _Sh = Sh.copy()
    _Sp = Mat_forward(lambda x: SymmIsoConstitutive(mu, lam, x), _E, _Ep)
    _Eb = Mat_reverse(lambda x: SymmIsoConstitutive(mu, lam, x), _E, _Sb)
    _Eh = Mat_hreverse(lambda x: SymmIsoConstitutive(mu, lam, x), _E, _Ep, _Sb, _Sh)
    print(f"\n//{case}")
    print(f"const T {case}_Sp_out[{ndim * ndim}] = {{", end='')
    print_flatten_16(symm2mat(_Sp), end='')
    print("};")
    print(f"const T {case}_Eb_out[{ndim * ndim}] = {{", end='')
    print_flatten_16(symm2mat(_Eb), end='')
    print("};")
    print(f"const T {case}_Eh_out[{ndim * ndim}] = {{", end='')
    print_flatten_16(symm2mat(_Eh), end='')
    print("};")

    # SymmIsoEnergy
    print("\nSymmIsoEnergy\n")
    _E = E.copy()
    _w = SymmIsoEnergy(mu, lam, _E)
    print(f"const T w_out = {_w:.16f};")

    # SymmIsoEnergy: AD_E_only
    case = "AD_E_only"
    _E = E.copy()
    _Eb = Eb.copy()
    _wb = Mat_forward(lambda x: SymmIsoEnergy(mu, lam, x), _E, _Eb)
    _Eb = Mat_reverse(lambda x: SymmIsoEnergy(mu, lam, x), _E, _wb)
    print(f"\n//{case}")
    print(f"const T {case}_wb_out = {_wb:.16f};")
    print(f"const T {case}_Eb_out[{ndim * ndim}] = {{", end='')
    print_flatten_16(symm2mat(_Eb), end='')
    print("};")

    # SymmIsoEnergy: AD_E_mu_lam
    case = "AD_E_mu_lam"
    _E = AugE.copy()
    _Eb = AugEb.copy()
    _wb = Mat_forward(AugSymmIsoEnergy, _E, _Eb)
    _Eb = Mat_reverse(AugSymmIsoEnergy, _E, _wb)
    print(f"\n//{case}")
    print(f"const T {case}_wb_out = {_wb:.16f};")
    print(f"const T {case}_mub_out = {_Eb[-2,0]:.16f};")
    print(f"const T {case}_lamb_out = {_Eb[-1,0]:.16f};")
    print(f"const T {case}_Eb_out[{ndim * ndim}] = {{", end='')
    print_flatten_16(symm2mat(_Eb[0:-2]), end='')
    print("};")

    # SymmIsoEnergy: A2D_E_only
    case = "A2D_E_only"
    _E = E.copy()
    _Ep = Eb.copy()
    _wb = sb
    _wh = sh
    _wp = Mat_forward(lambda x: SymmIsoEnergy(mu, lam, x), _E, _Ep)
    _Eb = Mat_reverse(lambda x: SymmIsoEnergy(mu, lam, x), _E, _wb)
    _Eh = Mat_hreverse(lambda x: SymmIsoEnergy(mu, lam, x), _E, _Ep, _wb, _wh)
    print(f"\n//{case}")
    print(f"const T {case}_wp_out = {_wp:.16f};")
    print(f"const T {case}_Eb_out[{ndim * ndim}] = {{", end='')
    print_flatten_16(symm2mat(_Eb), end='')
    print("};")
    print(f"const T {case}_Eh_out[{ndim * ndim}] = {{", end='')
    print_flatten_16(symm2mat(_Eh), end='')
    print("};")

    # SymmIsoEnergy: A2D_E_mu_lam
    case = "A2D_E_mu_lam"
    _E = AugE.copy()
    _Ep = AugEb.copy()
    _wb = sb
    _wh = sh
    _wp = Mat_forward(AugSymmIsoEnergy, _E, _Ep)
    _Eb = Mat_reverse(AugSymmIsoEnergy, _E, _wb)
    _Eh = Mat_hreverse(AugSymmIsoEnergy, _E, _Ep, _wb, _wh)
    print(f"\n//{case}")
    print(f"const T {case}_wp_out = {_wp:.16f};")
    print(f"const T {case}_Eb_out[{ndim * ndim}] = {{", end='')
    print_flatten_16(symm2mat(_Eb[0:-2]), end='')
    print("};")
    print(f"const T {case}_Eh_out[{ndim * ndim}] = {{", end='')
    print_flatten_16(symm2mat(_Eh[0:-2]), end='')
    print("};")
    print(f"const T {case}_mub_out = {_Eb[-2,0]:.16f};")
    print(f"const T {case}_muh_out = {_Eh[-2,0]:.16f};")
    print(f"const T {case}_lamb_out = {_Eb[-1,0]:.16f};")
    print(f"const T {case}_lamh_out = {_Eh[-1,0]:.16f};")

    # GreenStrain
    print("\nGreenStrain\n")
    _Ux = A
    _E = GreenStrain(_Ux)
    print(f"const T E_out[{ndim * ndim}] = {{", end='')
    print_flatten_16(symm2mat(_E), end='')
    print("};")

    # GreenStrain: AD
    case = "AD"
    _Ux = A.copy()
    _Uxb = Ab.copy()
    _Eb = Mat_forward(GreenStrain, _Ux, _Uxb)
    _Uxb = Mat_reverse(GreenStrain, _Ux, _Eb)
    print(f"\n//{case}")
    print(f"const T {case}_Eb_out[{ndim * ndim}] = {{", end='')
    print_flatten_16(symm2mat(_Eb), end='')
    print("};")
    print(f"const T {case}_Uxb_out[{ndim * ndim}] = {{", end='')
    print_flatten_16(_Uxb, end='')
    print("};")

    # GreenStrain: A2D
    case = "A2D"
    _Ux = A.copy()
    _Uxp = Ab.copy()
    _Eb = Sb.copy()
    _Eh = Sh.copy()
    _Ep = Mat_forward(GreenStrain, _Ux, _Uxp)
    _Uxb = Mat_reverse(GreenStrain, _Ux, _Eb)
    _Uxh = Mat_hreverse(GreenStrain, _Ux, _Uxp, _Eb, _Eh)
    print(f"\n//{case}")
    print(f"const T {case}_Uxb_out[{ndim * ndim}] = {{", end='')
    print_flatten_16(_Uxb, end='')
    print("};")
    print(f"const T {case}_Uxh_out[{ndim * ndim}] = {{", end='')
    print_flatten_16(_Uxh, end='')
    print("};")
    print(f"const T {case}_Ep_out[{ndim * ndim}] = {{", end='')
    print_flatten_16(symm2mat(_Ep), end='')
    print("};")

    # LinearGreenStrain
    print("\nLinerGreenStrain\n")
    _Ux = A
    _E = LinearGreenStrain(_Ux)
    print(f"const T E_out[{ndim * ndim}] = {{", end='')
    print_flatten_16(symm2mat(_E), end='')
    print("};")

    # LinearGreenStrain: AD
    case = "AD"
    _Ux = A.copy()
    _Uxb = Ab.copy()
    _Eb = Mat_forward(LinearGreenStrain, _Ux, _Uxb)
    _Uxb = Mat_reverse(LinearGreenStrain, _Ux, _Eb)
    print(f"\n//{case}")
    print(f"const T {case}_Eb_out[{ndim * ndim}] = {{", end='')
    print_flatten_16(symm2mat(_Eb), end='')
    print("};")
    print(f"const T {case}_Uxb_out[{ndim * ndim}] = {{", end='')
    print_flatten_16(_Uxb, end='')
    print("};")

    # LinearGreenStrain: A2D
    case = "A2D"
    _Ux = A.copy()
    _Uxp = Ab.copy()
    _Eb = Sb.copy()
    _Eh = Sh.copy()
    _Ep = Mat_forward(LinearGreenStrain, _Ux, _Uxp)
    _Uxb = Mat_reverse(LinearGreenStrain, _Ux, _Eb)
    _Uxh = Mat_hreverse(LinearGreenStrain, _Ux, _Uxp, _Eb, _Eh)
    print(f"\n//{case}")
    print(f"const T {case}_Uxb_out[{ndim * ndim}] = {{", end='')
    print_flatten_16(_Uxb, end='')
    print("};")
    print(f"const T {case}_Uxh_out[{ndim * ndim}] = {{", end='')
    print_flatten_16(_Uxh, end='')
    print("};")
    print(f"const T {case}_Ep_out[{ndim * ndim}] = {{", end='')
    print_flatten_16(symm2mat(_Ep), end='')
    print("};")

    return

def run_3d_cases():
    # Randomly generated inputs
    A = np.array([0.54881350, 0.71518937, 0.60276338,
                    0.54488318, 0.42365480, 0.64589411,
                    0.43758721, 0.89177300, 0.96366276]).reshape(3,3)

    B = np.array([0.63263687, 0.19519233, 0.13896359,
                    0.83668126, 0.50298017, 0.20476354,
                    0.57524785, 0.83518588, 0.58878908]).reshape(3,3)

    Ab = np.array([0.69626791, 0.53987667, 0.65374594,
                    0.45510150, 0.98931708, 0.74223721,
                    0.40945653, 0.74347344, 0.88483858]).reshape(3,3)

    Bb = np.array([0.14453098, 0.50253050, 0.44495442,
                    0.81206852, 0.59942161, 0.72829853,
                    0.99976353, 0.17018835, 0.75282306]).reshape(3,3)

    Cb = np.array([0.95793759, 0.25352624, 0.24823463,
                    0.26332010, 0.32297026, 0.40684991,
                    0.86658516, 0.82340961, 0.60516882]).reshape(3,3)

    Ch = np.array([0.66773557, 0.09141460, 0.95390316,
                    0.44180188, 0.53497833, 0.96218661,
                    0.45728781, 0.18113042, 0.98689272]).reshape(3,3)

    # Symmetric matrices
    S = np.array([0.16665530, 0.09586054, 0.29580571,
                 0.16117561, 0.23525668, 0.79122117]).reshape(6, 1)

    E = np.array([0.28147380, 0.47227010, 0.39857781,
                  0.92539657, 0.00320169, 0.15774955]).reshape(6, 1)

    Sb = np.array([0.14503781, 0.66755052, 0.83809867,
                   0.61659135, 0.99207311, 0.10221423]).reshape(6, 1)

    Eb = np.array([0.87502649, 0.41741302, 0.66047610,
                   0.84442668, 0.65391596, 0.25608671]).reshape(6, 1)

    Sh = np.array([0.52028303, 0.73053919, 0.44468893,
                   0.68831122, 0.54873182, 0.19163654]).reshape(6, 1)

    sb = 0.33622324
    sh = 0.73157930

    mu = 0.84010356
    lam = 0.76764533

    mub = 0.45061325
    lamb = 0.30576147

    AugE = np.append(E, [mu, lam]).reshape(len(E.flatten()) + 2, 1)
    AugEb = np.append(Eb, [mub, lamb]).reshape(len(E.flatten()) + 2, 1)

    ndim = 3
    generate_reference_values(ndim, A, B, Ab, Bb, Cb, Ch, S, E, Sb, Eb, Sh, sb,
                              sh, mu, lam, mub, lamb, AugE, AugEb)

def run_2d_cases():
    # Randomly generated inputs
    A  = np.array([0.96155402, 0.02855176, 0.95787560, 0.45439794]).reshape(2, 2)
    B  = np.array([0.80766462, 0.60212270, 0.86418474, 0.65304149]).reshape(2, 2)
    Ab = np.array([0.69626791, 0.53987667, 0.65374594, 0.45510150]).reshape(2, 2)
    Bb = np.array([0.14453098, 0.50253050, 0.44495442, 0.81206852]).reshape(2, 2)
    Cb = np.array([0.95793759, 0.25352624, 0.24823463, 0.26332010]).reshape(2, 2)
    Ch = np.array([0.66773557, 0.09141460, 0.95390316, 0.44180188]).reshape(2, 2)

    S  = np.array([0.16665530, 0.09586054, 0.29580571]).reshape(3, 1)
    E  = np.array([0.28147380, 0.47227010, 0.39857781]).reshape(3, 1)
    Sb = np.array([0.14503781, 0.66755052, 0.83809867]).reshape(3, 1)
    Eb = np.array([0.87502649, 0.41741302, 0.66047610]).reshape(3, 1)
    Sh = np.array([0.52028303, 0.73053919, 0.44468893]).reshape(3, 1)

    sb   = 0.33622324
    sh   = 0.73157930
    mu   = 0.84010356
    lam  = 0.76764533
    mub  = 0.45061325
    lamb = 0.30576147

    AugE = np.append(E, [mu, lam]).reshape(5, 1)
    AugEb = np.append(Eb, [mub, lamb]).reshape(5, 1)

    ndim = 2

    generate_reference_values(ndim, A, B, Ab, Bb, Cb, Ch, S, E, Sb, Eb, Sh, sb,
                              sh, mu, lam, mub, lamb, AugE, AugEb)

if __name__ == '__main__':
    # run_3d_cases()
    run_2d_cases()







