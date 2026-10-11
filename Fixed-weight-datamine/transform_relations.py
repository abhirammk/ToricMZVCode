"""Substitute MZV coordinates into toric polynomials, then reduce over ZZ."""
from sage.all import QQ, PolynomialRing, matrix, prod
from basis_coordinates import coordinate_system
from reduce_and_export_relations import reduce_export


def transport_relations(data, name, weight, descriptors, C):
    """Substitution only; do not evaluate the equations in the MZV algebra."""
    names = [f'B{i}' for i in range(C.ncols())]
    R = PolynomialRing(QQ, names=names)
    forms = [sum((c*x for c,x in zip(row,R.gens())),R.zero()) for row in C.rows()]
    relations = []
    for rel in data['binoms']:
        p = prod((f**e for f,e in zip(forms,rel['mon1'])),R.one()) - prod((f**e for f,e in zip(forms,rel['mon2'])),R.one())
        terms = [{'coefficient':str(c),'exponents':list(map(int,e))} for e,c in sorted(p.dict().items(),reverse=True)]
        relations.append({'terms':terms})
    result = {'weight':weight,'basis_name':name,
              'variables':[{'name':n,**d} for n,d in zip(names,descriptors)],
              'change_of_basis':{'matrix':[[str(c) for c in r] for r in C.rows()]},
              'relations':relations}
    return result


def convert(data, name='hoffmann', custom=None, progress=None):
    weight, descriptors, C = coordinate_system(data,name,custom,progress)
    result = transport_relations(data,name,weight,descriptors,C)
    result = reduce_export(result,progress=progress)
    verify(data,result)
    return result


def verify(data, result):
    assert 'source' not in result
    C = matrix(QQ, result['change_of_basis']['matrix'])
    degrees = {sum(r['mon1']) for r in data['binoms']}
    if len(degrees) > 1:
        for degree in sorted(degrees):
            subset = dict(data, binoms=[r for r in data['binoms'] if sum(r['mon1']) == degree])
            selected = dict(result, relations=[r for r in result['relations'] if r['terms'] and sum(r['terms'][0]['exponents']) == degree])
            verify(subset, selected)
        return
    # Quadratic inverse substitution is D^T Q D. Matrix arithmetic avoids
    # repeatedly expanding thousands of products of dense linear forms.
    if (all(sum(r['mon1']) == sum(r['mon2']) == 2 for r in data['binoms'])
            and all(sum(t['exponents']) == 2 for r in result['relations'] for t in r['terms'])):
        n = C.nrows()
        D = C.inverse()
        pairs = [(i,j) for i in range(n) for j in range(i,n)]
        originals = []
        for relation in data['binoms']:
            coefficients = {}
            for sign,side in ((1,'mon1'),(-1,'mon2')):
                indices = [i for i,e in enumerate(relation[side]) for _ in range(e)]
                key = tuple(indices)
                coefficients[key] = coefficients.get(key,0)+sign
            originals.append([coefficients.get(pair,0) for pair in pairs])
        recovered = []
        for relation in result['relations']:
            Q = matrix(QQ,n,n)
            for term in relation['terms']:
                c = QQ(term['coefficient'])
                assert c.denominator() == 1
                i,j = [i for i,e in enumerate(term['exponents']) for _ in range(e)]
                if i == j:
                    Q[i,i] += c
                else:
                    Q[i,j] += c/2
                    Q[j,i] += c/2
            original_Q = D.transpose()*Q*D
            recovered.append([original_Q[i,j]*(1 if i == j else 2) for i,j in pairs])
        A = matrix(QQ,originals)
        B = matrix(QQ,recovered)
        assert B.row_space() == A.row_space()
        assert B.nrows() == A.rank()
        return
    R = PolynomialRing(QQ, names=[v['name'] for v in result['variables']])
    # For the few higher-degree generators, expand forward rather than
    # substituting dense inverse forms into tens of thousands of terms.
    # Invertibility of C makes equality of target spans equivalent to
    # equality of source spans.
    assert C.rank() == C.nrows() == C.ncols()
    forms = [sum((c*x for c,x in zip(row,R.gens())),R.zero()) for row in C.rows()]
    recovered_relations = []
    original_relations = []
    for converted in result['relations']:
        assert all(QQ(t['coefficient']).denominator() == 1 for t in converted['terms'])
        recovered_relations.append({tuple(t['exponents']): QQ(t['coefficient']) for t in converted['terms']})
    for original in data['binoms']:
        expected = prod((x**e for x,e in zip(forms,original['mon1'])),R.one())-prod((x**e for x,e in zip(forms,original['mon2'])),R.one())
        original_relations.append({tuple(e): c for e,c in expected.dict().items()})
    columns = sorted({e for row in recovered_relations+original_relations for e in row})
    def rows(dictionaries):
        return matrix(QQ, len(dictionaries), len(columns),
                      [row.get(e,0) for row in dictionaries for e in columns])
    A, B = rows(original_relations), rows(recovered_relations)
    assert B.row_space() == A.row_space()
    assert B.nrows() == A.rank()


