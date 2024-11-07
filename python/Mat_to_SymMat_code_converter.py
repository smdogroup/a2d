def convert(s, mat):
    s_mod = ""

    # general mat index -> SymMat index
    mapper = {
        "0": "0",
        "1": "1",
        "2": "3",
        "3": "1",
        "4": "2",
        "5": "4",
        "6": "3",
        "7": "4",
        "8": "5",
    }

    for i, c in enumerate(s):
        if c in mapper and s[i - 2] == mat:
            c = mapper[c]
        s_mod += c

    print("========== Before ==========")
    print(s)

    print("========== After ==========")
    print(s_mod)


if __name__ == "__main__":
    s = """
    C[0] = A[0] * S[0] + A[1] * S[3] + A[2] * S[6];
    C[1] = A[0] * S[1] + A[1] * S[4] + A[2] * S[7];
    C[2] = A[0] * S[2] + A[1] * S[5] + A[2] * S[8];
    C[3] = A[3] * S[0] + A[4] * S[3] + A[5] * S[6];
    C[4] = A[3] * S[1] + A[4] * S[4] + A[5] * S[7];
    C[5] = A[3] * S[2] + A[4] * S[5] + A[5] * S[8];
    C[6] = A[6] * S[0] + A[7] * S[3] + A[8] * S[6];
    C[7] = A[6] * S[1] + A[7] * S[4] + A[8] * S[7];
    C[8] = A[6] * S[2] + A[7] * S[5] + A[8] * S[8];
  } else if constexpr (opA == MatOp::TRANSPOSE) {
    C[0] = A[0] * S[0] + A[3] * S[3] + A[6] * S[6];
    C[1] = A[0] * S[1] + A[3] * S[4] + A[6] * S[7];
    C[2] = A[0] * S[2] + A[3] * S[5] + A[6] * S[8];
    C[3] = A[1] * S[0] + A[4] * S[3] + A[7] * S[6];
    C[4] = A[1] * S[1] + A[4] * S[4] + A[7] * S[7];
    C[5] = A[1] * S[2] + A[4] * S[5] + A[7] * S[8];
    C[6] = A[2] * S[0] + A[5] * S[3] + A[8] * S[6];
    C[7] = A[2] * S[1] + A[5] * S[4] + A[8] * S[7];
    C[8] = A[2] * S[2] + A[5] * S[5] + A[8] * S[8];
    """

    mat = "S"

    convert(s, mat)
