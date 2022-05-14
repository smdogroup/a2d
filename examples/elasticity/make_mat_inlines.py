
def get_mat(s, mator='normal'):
    mat = []
    if mator == 'normal':
        for k in [
            (0, 0), (0, 1), (0, 2),
            (1, 0), (1, 1), (1, 2),
            (2, 0), (2, 1), (2, 2)]:
            mat.append('%s(%d, %d)'%(s, k[0], k[1]))
    elif mator == 'transpose':
        for k in [
            (0, 0), (1, 0), (2, 0),
            (0, 1), (1, 1), (2, 1),
            (0, 2), (1, 2), (2, 2)]:
            mat.append('%s(%d, %d)'%(s, k[0], k[1]))
    elif mator == 'symmetric':
        for k in [
            (0, 0), (0, 1), (0, 2),
            (1, 0), (1, 1), (1, 2),
            (2, 0), (2, 1), (2, 2)]:
            mat.append('%s(%d, %d)'%(s, k[0], k[1]))
    return mat

def get_subroutine(Aor='normal', Bor='normal', op='=', scale=False):
    A = get_mat('A', Aor)
    B = get_mat('B', Bor)
    C = get_mat('C', 'normal')

    product = []
    for i in range(3):
        for j in range(3):
            s = ''
            for k in range(3):
                s += '%s * %s'%(A[3*i + k], B[3*k + j])
                if k < 2:
                    s += ' + '
            product.append(s)

    name = ''
    if Aor == 'normal':
        name = 'Mat3x3'
    elif Aor == 'transpose':
        name = 'MatTrans3x3'
    elif Aor == 'symmetric':
        name = 'Symm3x3'

    if Bor == 'normal':
        name += 'Mat'
    elif Bor == 'transpose':
        name += 'MatTrans'
    elif Bor == 'symmetric':
        name += 'Symm'

    name += 'Mult'

    if op == '+=':
        name += 'Add'
    elif op == '-=':
        name += 'Sub'

    args = 'const AMatType& A, const BMatType& B, CMatType& C'
    if scale:
        args = 'T scale, ' + args
        name += 'Scale'

    name += 'Core'

    if scale:
        result = 'template<typename T, class AMatType, class BMatType, class CMatType>\n'
    else:
        result = 'template<class AMatType, class BMatType, class CMatType>\n'

    result += 'inline void %s( %s ){\n'%(name, args)
    for k in range(9):
        if scale:
            result += '  %s %s scale * (%s);\n'%(C[k], op, product[k])
        else:
            result += '  %s %s (%s);\n'%(C[k], op, product[k])

    result += '}\n\n'

    return result

or_pairs = [
    ['symmetric', 'symmetric'],
    ['symmetric', 'normal'],
    ['symmetric', 'transpose'],
    ['normal', 'symmetric'],
    ['transpose', 'symmetric'],

    ['normal', 'normal'],
    ['normal', 'transpose'],
    ['transpose', 'normal'],
    ['transpose', 'transpose']]

kwargs_list = [
    {'op':'=', 'scale':False},
    {'op':'=', 'scale':True},
    {'op':'+=', 'scale':False},
    {'op':'-=', 'scale':False},
    {'op':'+=', 'scale':True},]

s = ''
for pair in or_pairs:
    for kwargs in kwargs_list:
        s += get_subroutine(Aor=pair[0], Bor=pair[1], **kwargs)

print(s)
