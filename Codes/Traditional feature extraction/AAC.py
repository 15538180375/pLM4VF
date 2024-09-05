import re
from collections import Counter

def AAC(fastas, **kw):         
    if kw['order'] != None:
        AA = kw['order']
    else:
        'ACDEFGHIKLMNPQRSTVWY'

    encodings = []
    header = ['#']
    for i in AA:
        header.append(i)
    encodings.append(header)

    for i in fastas:
        name, sequence = i[0], re.sub('-','',i[1])
        count = Counter(sequence)
        for key in count:
            count[key] = count[key]/len(sequence)
        code = [name]
        for aa in AA:
            code.append(count[aa])
        encodings.append(code)
    return encodings



