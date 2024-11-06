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
    C[0] = S[0] * B[0] + S[1] * B[3] + S[2] * B[6];
    C[1] = S[0] * B[1] + S[1] * B[4] + S[2] * B[7];
    C[2] = S[0] * B[2] + S[1] * B[5] + S[2] * B[8];
    C[4] = S[3] * B[1] + S[4] * B[4] + S[5] * B[7];
    C[5] = S[3] * B[2] + S[4] * B[5] + S[5] * B[8];
    C[8] = S[6] * B[2] + S[7] * B[5] + S[8] * B[8];
  } else if constexpr (opB == MatOp::TRANSPOSE) {
    C[0] = S[0] * B[0] + S[1] * B[1] + S[2] * B[2];
    C[1] = S[0] * B[3] + S[1] * B[4] + S[2] * B[5];
    C[2] = S[0] * B[6] + S[1] * B[7] + S[2] * B[8];
    C[4] = S[3] * B[3] + S[4] * B[4] + S[5] * B[5];
    C[5] = S[3] * B[6] + S[4] * B[7] + S[5] * B[8];
    C[8] = S[6] * B[6] + S[7] * B[7] + S[8] * B[8];
    """

    mat = "S"

    convert(s, mat)
