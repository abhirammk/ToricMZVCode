"""Generate integer LLL relation files and reusable MZV coordinates.

Run: sage -python generate_relation_files.py [--weight 10] [--basis hoffmann]
Default: weights 10–18, both bases; keep complete outputs and reuse coordinates.
"""
import argparse
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import tempfile
from sage.all import pari
from basis_coordinates import compute_coordinates, read_json
from transform_relations import convert, verify
from reduce_and_export_relations import write_export, basis_definitions
from product_basis_cache import with_cached_product_bases

ROOT = Path(__file__).resolve().parent


def toric_source(weight):
    source_dir = ROOT/'toric-source'
    path = source_dir/f'toric-relations-weight-{weight}.json'
    if not path.exists():
        executable = os.environ.get('M2_EXECUTABLE') or shutil.which('M2') or '/Applications/Macaulay2/bin/M2'
        with tempfile.TemporaryDirectory(prefix='mzv-m2-') as temporary:
            script = Path(temporary)/'compute_toric_relations.m2'
            template = (source_dir/'compute_toric_relations.m2').read_text()
            script.write_text(re.sub(r'(?m)^j[ \t]*=[ \t]*[0-9]+[ \t]*$',f'j = {weight}',template,count=1))
            subprocess.run([executable,'--script',str(script)],cwd=source_dir,check=True)
    return read_json(path)


def save_coordinates(custom, path):
    path.parent.mkdir(parents=True,exist_ok=True)
    temporary = path.with_suffix('.json.tmp')
    temporary.write_text(json.dumps(custom,indent=2)+'\n',encoding='utf-8')
    temporary.replace(path)


@with_cached_product_bases
def generate(weight, name, output_dir=ROOT, custom=None, force=False, progress=print):
    """Produce one coordinate file and one verified JSON/text relation pair."""
    output_dir = Path(output_dir)
    prefix = output_dir/f'{name}-basis'/f'relations-weight-{weight}-{name}'
    coordinate_path = output_dir/'mzv-linear-bases'/f'change-of-basis-weight-{weight}-{name}.json'
    complete = all(prefix.with_suffix(s).exists() for s in ('.json','.txt'))
    if complete and coordinate_path.exists() and not force and custom is None:
        progress(f'Weight {weight}, {name}: keeping completed files')
        return
    source = toric_source(weight)
    coordinates = custom
    if coordinates is None:
        cached = coordinate_path if coordinate_path.exists() else ROOT/'mzv-linear-bases'/coordinate_path.name
        if cached.exists():
            coordinates = read_json(cached)
            if coordinates.get('name') != name:
                raise ValueError(f'basis name does not match {cached}')
    if coordinates is None or 'coordinates' not in coordinates:
        coordinates = compute_coordinates(source,name,coordinates,progress=progress)
    save_coordinates(coordinates,coordinate_path)
    # Missing coordinates may have been restored for an existing export. Check
    # that it is the same target basis and relation span before keeping it.
    if complete and not force:
        saved = read_json(prefix.with_suffix('.json'))
        variables = [{'name':f'B{i}','indices':w} for i,w in enumerate(coordinates['basis'])]
        if saved['basis'] == basis_definitions({'variables':variables}):
            verify(source,dict(saved,variables=variables,change_of_basis={'matrix':coordinates['coordinates']}))
            progress(f'Weight {weight}, {name}: coordinates restored; existing exports verified')
            return
    progress(f'Weight {weight}, {name}: converting and reducing relations')
    result = convert(source,name,coordinates,progress=progress)
    write_export(result,prefix)
    progress(f'Weight {weight}, {name}: exported {len(result["relations"])} verified integer relations')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--weight',type=int,help='generate a single weight')
    parser.add_argument('--min-weight',type=int)
    parser.add_argument('--max-weight',type=int)
    parser.add_argument('--basis',choices=['hoffmann','depth-minimal','both'],default='both')
    parser.add_argument('--basis-file',help='custom JSON basis with optional coordinates')
    parser.add_argument('--force',action='store_true',help='regenerate relation files using saved coordinates')
    parser.add_argument('--output-dir',type=Path,default=ROOT)
    args = parser.parse_args()
    if args.weight is not None and (args.min_weight is not None or args.max_weight is not None):
        parser.error('--weight cannot be combined with a weight range')
    minimum = args.weight if args.weight is not None else (args.min_weight if args.min_weight is not None else 10)
    maximum = args.weight if args.weight is not None else (args.max_weight if args.max_weight is not None else 18)
    if not 2 <= minimum <= maximum <= 18:
        parser.error('require 2 <= min-weight <= max-weight <= 18')
    custom = read_json(args.basis_file) if args.basis_file else None
    if custom:
        if minimum != maximum:
            parser.error('a custom basis requires a single weight')
        names = [custom.get('name','custom')]
    else:
        names = ['hoffmann','depth-minimal'] if args.basis == 'both' else [args.basis]
    if any(not re.fullmatch(r'[A-Za-z0-9_-]+',name) for name in names):
        parser.error('basis names must contain only letters, digits, underscores and hyphens')
    pari.allocatemem(268435456,8589934592)
    for weight in range(minimum,maximum+1):
        for name in names:
            generate(weight,name,args.output_dir,custom,args.force,progress=lambda s:print(s,flush=True))


if __name__ == '__main__':
    main()
