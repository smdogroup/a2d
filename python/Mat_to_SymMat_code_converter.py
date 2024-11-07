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
        if c in mapper and s[i - 2] in mat:
            c = mapper[c]
        s_mod += c

    print("========== Before ==========")
    print(s)

    print("========== After ==========")
    print(s_mod)


if __name__ == "__main__":
    s = """
    Sh[0] += (S[8] * Sp[4] - S[7] * Sp[5] + Sp[8] * S[4] - Sp[7] * S[5]) * bdet;
    Sh[1] += (S[6] * Sp[5] - S[8] * Sp[3] + Sp[6] * S[5] - Sp[8] * S[3]) * bdet;
    Sh[2] += (S[7] * Sp[3] - S[6] * Sp[4] + Sp[7] * S[3] - Sp[6] * S[4]) * bdet;
    Sh[4] += (S[8] * Sp[0] - S[6] * Sp[2] + Sp[8] * S[0] - Sp[6] * S[2]) * bdet;
    Sh[5] += (S[6] * Sp[1] - S[7] * Sp[0] + Sp[6] * S[1] - Sp[7] * S[0]) * bdet;
    Sh[8] += (S[0] * Sp[4] - S[3] * Sp[1] + Sp[0] * S[4] - Sp[3] * S[1]) * bdet;

    Sh[0] += (S[8] * S[4] - S[7] * S[5]) * hdet;
    Sh[1] += (S[6] * S[5] - S[8] * S[3]) * hdet;
    Sh[2] += (S[7] * S[3] - S[6] * S[4]) * hdet;
    Sh[4] += (S[8] * S[0] - S[6] * S[2]) * hdet;
    Sh[5] += (S[6] * S[1] - S[7] * S[0]) * hdet;
    Sh[8] += (S[0] * S[4] - S[3] * S[1]) * hdet;
    """

    mat = {"S", "p", "h"}

    convert(s, mat)
