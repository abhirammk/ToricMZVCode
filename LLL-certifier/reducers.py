"""Untrusted external reducers; certificates are always checked independently."""
import json
import os
import re
import subprocess


SIZE_POWERS = {'': 0, 'K': 1, 'M': 2, 'G': 3}
STACK_NOTICE = re.compile(
    r'\s*\*\*\*.*Warning: (?:new maximum stack size = [0-9]+ '
    r'\([^\n]+\)|increasing stack size to [0-9]+)\.\s*'
)

def gp_stack_bytes(value):
    match = re.fullmatch(r'([1-9][0-9]*)([KMG]?)', value.upper())
    if not match:
        raise ValueError('GP_STACK and GP_STACK_MAX must be positive sizes such as 256M or 4G')
    return int(match[1]) * 1024 ** SIZE_POWERS[match[2]]


def run_gp(script, gp):
    initial = gp_stack_bytes(os.environ.get('GP_STACK', '256M'))
    maximum = max(initial, gp_stack_bytes(os.environ.get('GP_STACK_MAX', '4G')))
    # Changing parisizemax can reset GP's stack, so use a separate input line.
    source = f'default(parisizemax,{maximum});\n'+script
    process = subprocess.run([gp, '-fq', '-s', str(initial)],
                             input=source, capture_output=True, text=True, check=True)
    # GP reports successful automatic growth on stderr; retain all other
    # diagnostics as errors (GP can exit zero even after a computation fails).
    diagnostics = [
        line for line in process.stderr.splitlines()
        if line.strip() and not STACK_NOTICE.fullmatch(line)
    ]
    if diagnostics:
        message = '\n'.join(diagnostics)
        if 'stack overflows' in message or 'stack overflow' in message:
            message += '\nIncrease GP_STACK_MAX (for example GP_STACK_MAX=8G).'
        raise ValueError(f'GP reduction failed: {message}')
    return process.stdout


def transform(rounded, multiplier, gp='gp'):
    """Return an integer column transform; the verifier checks unimodularity."""
    n = len(rounded)
    source = (f'r={rounded}; K={multiplier}; n=#r; '
              'M=matrix(n+1,n,i,j,if(i<=n,if(i==j,K,0),r[j])); '
              'U=qflll(M,1); print(vector(n,i,Vec(U[i,]))); quit\n')
    U = json.loads(run_gp(source, gp))
    if (type(U) is not list or len(U) != n
            or any(type(row) is not list or len(row) != n for row in U)):
        raise ValueError('reducer returned a matrix of the wrong size')
    if any(type(v) is not int for row in U for v in row):
        raise ValueError('reducer returned noninteger transform entries')
    return U


