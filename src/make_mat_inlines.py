
def get_mat(s, mator='normal'):
    mat = []
    if mator == 'normal':
        for k in range(9):
            mat.append('%s[%d]'%(s, k))
    elif mator == 'transpose':
        # 0 3 6
        # 1 4 7
        # 2 5 8
        for k in [0, 3, 6, 1, 4, 7, 2, 5, 8]:
            mat.append('%s[%d]'%(s, k))
    elif mator == 'symmetric':
        # 0 1 2
        # 1 3 4
        # 2 4 5
        for k in [0, 1, 2, 1, 3, 4, 2, 4, 5]:
            mat.append('%s[%d]'%(s, k))
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

    args = 'const TacsScalar A[], const TacsScalar B[], TacsScalar C[]'
    if scale:
        args = 'TacsScalar scale, ' + args
        name += 'Scale'

    name += 'Core'

    result = 'inline void %s( %s ){\n'%(name, args)
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
